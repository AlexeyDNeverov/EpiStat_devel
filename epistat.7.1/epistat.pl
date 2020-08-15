#!/usr/bin/env perl
#This script calculates epistatic statistics [Kryazhimsky11]
#options:
#	-x <FN> - input file for EpiStat application in XPAR format [Kryazhimsky11]. Contains tree in the adjacent file format
#	[-m <"0|1|2|3">]  PrintMatrixMode - mode for printing of mutation matrix  (see below): 0- don't print;1 - background;2 - target; 3 - both
#			Default="0"
#	[-a <"dist(ances)"+"epistat"+"figtree">] Separator='+', Run analysis:
#			Default="epistat"
#		epistat - calculate epistatic statistics
#		distances - average distances for site pairs: root-site_pair, bgr-fgr, site_pair-closest terminal, site_pair-farest terminal, site_pair-average terminal
#		figtree_site_pairs - draw mutations on branches for specified site pairs in FigTree format. This option requires site pair file specified by '-s' option
#		figtree_subst_dens - colour branches on the tree according to differences of observed and expected (assuming uniformity) numbers of substitutions in their subtrees.
#	[-s <site_pair_fn>] - Specifies a file (CSV) with a list of site pairs those mutations will be drawn on the tree
#	-p - do not print site pairs
#Params:
#<parameters> - The file with script params.
#		Contents:
#		[TAU="<ufloat>"] - timescale parameter (>0)
#		PairsType=("intra[gene]"|"inter[gene]")
#		[StatType=("all_pairs"|"independent_pairs"|"branch_pairs")] - switches the way of convolution into a site pair epistatic statistics of epistatic statistics of corresponding mutation pairs
#			"all_pairs" - all consequtive mutation pairs are treated independently of each other
#			"independent_pairs" - average epistatic statistics for mutation pairs having common background mutation
#			"branch_pairs" - only events on same branches are accounted. This option is adjusted for searching for genotype-phenotype associations (GWAS)
#			"all_time_determined_pairs" - the same as the "all_pairs" mode excluding the mutation pairs on the same branch
#			Default="all_pairs"
#		[StatFunc=("id(entity)*"|"exp")] - specifies a transformation function for distances between mutations
#			Default="exp" - exponential with weighting parameter tau
#		[IgnoreTerminals="0|1"] - specifies whether to ignore terminal branches
#			Default="1"
#		[OutSitePairs="Suffix"] - list of site pairs for which epistatic statistics is calculated
#			Default=".site_pairs"
#		[OutEpiStatistics="Suffix"] - file with epistatic statistics for pairs
#			Default=".stat.exp".$TAU
#		[OutMutMatrix="Suffix"] - binary matrix: raws - branches; columns - sites
#			Default=".mmtx"
#		[OutIdx2Site="Suffix"] - table to convert columns' indexes of OutMutMatrix into site positions
#			Default=".idx2site"
#		[OutIdx2Branch="Suffix"] - table to convert raws' indexes of OutMutMatrix into branch names
#			Default=".idx2branch"
#		[PrintAllPairs="0|1"] - specifies whether to print ALL site pairs (even not consequtive)
#		[BranchDistUnits="(DEFAULT||SYN|NSYN|TOTAL)"] - Units for branch length distances
#			Default="SYN"
#		[BranchConfidences="<FN>"] - A file containing confidences for tree branches
#		[BranchConfidenceThreshold="<ufloat>"] - A threshold to accept of a branch support

#server inst
use lib "$ENV{EPISTAT_HOME}";
use lib "$ENV{EPISTAT_LIB}";

use strict;
use Bio::Phylo::IO;
use Class::Struct;
use Getopt::Std;
use File::Basename;
use EpistatAnalysis::PhyloDistances qw(mk_distance_analysis);
use EpistatAnalysis::PhyloDrawSites qw(site_pairs2figtree_str sites2figtree_str subst_distrib2figtree_str);
use List::BinarySearch qw(binsearch);
use SitePairInfo;
use Pairs::ConsequtivePairs;
use Pairs::BranchPairs;
use ParallelSubstWeights;
use BranchSiteMutationsMatrix qw(print_branch_site_matrix);

my %args;
if(!getopts('x:m:a:ps:',\%args)){
	die "\nError in option string!";
}

my $xpar_fn=$args{x};
die "\nThe input XPAR file is not specified! Use -x option!" unless defined $xpar_fn;
my $branch_conf_fn;
my %branch_confidences;
my $branch_conf_threshold;
my $stat_fn;
my $pairs_fn=".site_pairs";
my $idx2site_fn=".idx2site";
my $idx2branch_fn=".idx2branch";
my $mmtx_fn=".mmtx";
my $tau;
my $f_intragene;
my $f_ignore_terminals=1;
my $f_mtx_mode=0;
my $f_print_all_pairs=0;
my $f_print_pairs=1;
my $f_stat_type=0; #all_pairs
my $f_stat_func=1; #exp
my %f_analysis;
my $dist_fn=".analysis.dist";
my $input_site_pair_fn;
my $npheno;
my @pheno_labels;

if(defined($args{s})){
	$input_site_pair_fn=$args{s};
}

if(defined($args{a})){
	my @splitter=split '\+',$args{a};
	foreach my $str(@splitter){
		if($str=~m/^dist(ances)?/i){
			$f_analysis{dist}=1; 
		}elsif($str=~m/^epistat/i){
			$f_analysis{epistat}=1; 
		}elsif($str=~m/^figtree_site_pairs/i){
			$f_analysis{figtree}=1;
			die "\nList of site pairs for -a figtree_site_pairs analysis was not accounted!" unless defined $input_site_pair_fn;
		}elsif($str=~m/^figtree_subst_dens/i){
			$f_analysis{subst_dens}=1;
			die "\nList of site pairs for -a figtree_subst_dens analysis was not accounted!" unless defined $input_site_pair_fn;
		
		}else{
			die "\nUnknown analysis: -a $str!";
		}
	}
}else{
	$f_analysis{epistat}=1;
}

$f_print_pairs=0 if(defined $args{p});

if(defined $args{m}){
	$f_mtx_mode=$args{m};
	die "\nWrong output matrix mode ($f_mtx_mode) is accounted in the -m option!\n\tUINT in range is [0..3] expected!" unless $f_mtx_mode=~/^[0123]\s*$/;
}		
#0 - default
#1 - number of synonimous substitutions
#2 - number of non synonimous substitutions
#3 - total number of substitutions
my $distance_mode=1;	
open INFILE, "$ARGV[0]" or die "\nUnable to open input file: $ARGV[0]!";
while(<INFILE>){
	$_=$` if(/#/);
	if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
		my ($key,$value)=($1,$2);
		if($key eq "OutSitePairs"){
			$pairs_fn=$value;
		}elsif($key eq "OutEpiStatistics"){
			$stat_fn=$value;
		}elsif($key eq "IgnoreTerminals"){
			$f_ignore_terminals=$value;
			die "\nUnknown value ($value) for the IgnoreTerminals parameter in file $ARGV[0]!\n\tUINT in range is [0..1] expected!" unless $value=~/^[01]\s*$/;
		}elsif($key eq "OutMutMatrix"){
			$mmtx_fn=$value;
		}elsif($key eq "OutIdx2Site"){
			$idx2site_fn=$value;
		}elsif($key eq "OutIdx2Branch"){
			$idx2branch_fn=$value;
		}elsif($key eq "PairsType"){
			if($value=~m/^intra(gene)*/i){
				$f_intragene=1;
			}elsif($value=~m/^inter(gene)*/i){
				$f_intragene=0;
			}else{
				die "\nUnknown value ($value) for the PairsType parameter in file $ARGV[0]!"; 
			};
		}elsif($key eq "PrintAllPairs"){
			$f_print_all_pairs=$value;
			die "\nUnknown value ($value) for the PrintAllPairs parameter in file $ARGV[0]!\n\tUINT in range is [0..1] expected!" unless $value=~/^[01]\s*$/;
		}elsif($key eq "TAU"){
			if($value=~/\d+(\.\d+)*/){
				$tau=$value;
			}else{
				die "\nWrong value for the TAU parameter in file $ARGV[0]!\n\tPositive number expected!";
			}
		}elsif($key eq "BranchDistUnits"){
			if($value=~m/^default/i){
				$distance_mode=0;
			}elsif($value=~m/^syn/i){
				$distance_mode=1;
			}elsif($value=~m/^nsyn/i){
				$distance_mode=2;
			}elsif($value=~m/^total/i){
				$distance_mode=3;
			}else{
				die "\nUndefined distance mode -m $args{u}";
			}
		}elsif($key eq "StatType"){
			if($value=~m/^all_pairs/i){
				$f_stat_type=0;
			}elsif($value=~m/^independent_pairs/i){
				$f_stat_type=1;
			}elsif($value=~m/^branch_pairs/i){
				$f_stat_type=2;
			}elsif($value=~m/^all_time_determined_pairs/i){
				$f_stat_type=3;
			}else{
				die "\nUndefined type of statistics: StatType=$value!";
			}
		}elsif($key eq"StatFunc"){
			if($value=~m/^id(entity)?$/i){
				$f_stat_func=0;
			}elsif($value=~m/^exp/i){
				$f_stat_func=1;
			}else{
				die "\nUndefined type of statistics: StatFunc=$value!";
			}
		}elsif($key eq "BranchConfidences"){
			$branch_conf_fn=$value;
		}elsif($key eq "BranchConfidenceThreshold"){
			$branch_conf_threshold=$value;
		}elsif($key eq "Phenotypes"){
			if($value==1){
				$npheno=0;
			}else{
				die "\nUnpropper value for the 'Phenotype' parameter! Valid values: 0 or 1" unless $value==0;
			}
		}else{
			die "\nUnknown parameter: $key in the input file $ARGV[0]!";
		}
	}
}
close INFILE;

if(defined($branch_conf_fn)){
	die "\nUndefined the minimal accepted value for branches supports!" unless defined $branch_conf_threshold;
	open INFILE, "<$branch_conf_fn" or die "\nUnable to open input file: $branch_conf_fn!";
	while(<INFILE>){
		chomp;
		if(/^(\S+)\t([-+]?\d*\.?\d+([eE][-+]?\d+)?)/){
			$branch_confidences{$1}=$2;
		}
	}
	close INFILE;
}

die "\nError: Parameter PairsType is undefined in the file $ARGV[0]!" unless defined $f_intragene;

unless(defined $stat_fn){
	$stat_fn=".stat";
	if($f_stat_func==1){
		$stat_fn.=".exp";
		die "\nError: The TAU parameter has no defaul value!" unless defined($tau)||($f_stat_type==2);
		$stat_fn.=".$tau";
	}elsif($f_stat_func==0){
		$stat_fn.=".ident";
	}
}

my ($basename,$dir,$ext) = fileparse($xpar_fn,'\.[^\.]*$');
$pairs_fn=$dir.$basename.$pairs_fn;
$stat_fn=$dir.$basename.$stat_fn;
$dist_fn=$dir.$basename.$dist_fn if $f_analysis{dist};
my $tree=Bio::Phylo::IO->parse(
    -file => $xpar_fn,
    -format => 'adjacency');
#!!!The return value of Adjacency parser not corresponded to interface of the Bio::Phylo::IO->parse
$tree=$tree->[0]->first;
my %terminals;
if($f_ignore_terminals){
	foreach my $node(@{$tree->get_terminals}){
		my $name=$node->get_name;
		$terminals{$name}=1;
	};
};

#my $str=Bio::Phylo::IO->unparse(-phylo => $tree,-format => 'newick');
#print $str;

my %site_pairs;
my %bg_site_idx; #converts background site coordinate into array index
my @bg_site2cord;
my %tg_site_idx; #converts target site coordinate into array index
my @tg_site2cord;
#substitution counts
my %bg_nsubst;
my %tg_nsubst;
 
#substitutions
my %tg_subst_map;
my %bg_subst_map;
my $rh_bg_subst_branches;
my $rh_tg_subst_branches;
#phenotypes
my %phenotypes;
	
struct SubstInfo => {
	site => '$',
	weight => '$',
	bases => '@'
};

sub parse_subst_abbr{
	my $str=shift;
	my $ra_bases=shift;
	my $spos;
	if($str=~/([A-Za-z-])(\d+)([[A-Za-z-])/){
		@{$ra_bases}=($1,$3) if(defined $ra_bases);
		$spos=$2;
	}else{
		$str=~s/\D//g;
		$spos=$str;
	};
	return $spos;
}

sub calc_subst_numbers{
	my $sp_info=shift;
	my %bg_branches;
	my %tg_branches;
	foreach my $sb_info(@{$sp_info->subst_info}){
		$bg_branches{$sb_info->bg_branch()->get_name()}=1;
		$tg_branches{$sb_info->tg_branch()->get_name()}=1;
	};
	return (scalar(keys %bg_branches),scalar(keys %tg_branches));
}

sub init_subst_map{
	my $xpar_fn=shift;
	my $tree=shift;
	my ($rh_bg_subst_branches,$rh_tg_subst_branches)=@_;
	open INPF, "<$xpar_fn" or die "\nUnable to open input file $xpar_fn!";
	$_=<INPF>; #skip header line
	chomp;
	s/\s+$//;
	my @line=split "\t";
	if(defined $npheno){
		if($line[6]=~/^phenotypes:/){
			@pheno_labels=split ",", $';
			$npheno=@pheno_labels;
		}else{
			"\nNo phenotypes defined on branches of a tree in the file: $xpar_fn!";
		}
	}
	die "\nHeader line was omitted in file: $xpar_fn!" unless /^child\tparent\tlength/;
	$_=<INPF>; #skip root node statement
	die "\nThe root node statement required in file: $xpar_fn!" unless /^[\w.]+\s*$/;
	my %syn_counts;
	my %nsyn_counts;
	my $I=0;
	my $J=0;
	while(<INPF>){
		chomp;
		s/\s+$//;
		@line=split "\t";
		next if $f_ignore_terminals&&$terminals{$line[0]};
		my $n=0;
		if($line[4] ne ""){
			my @sites=split(/;/, $line[4]);
			$n=@sites;
			$tg_subst_map{$line[0]}={};
			foreach my $str(@sites){
				my @bases;
				my $spos=parse_subst_abbr($str,\@bases);
				if(!defined $tg_site_idx{$spos}){
					$tg_site_idx{$spos}=$J;
					$tg_site2cord[$J++]=$spos;
				};
				$tg_nsubst{$spos}++;
				my $si=SubstInfo->new();
				$si->site($spos);
				$si->bases([@bases]) if @bases;
				$tg_subst_map{$line[0]}->{$spos}=$si;
				if(defined $rh_tg_subst_branches){
					$rh_tg_subst_branches->{$spos}=[] unless defined $rh_tg_subst_branches->{$spos};
					push @{$rh_tg_subst_branches->{$spos}},$line[0];
				}
			}
		}
		$nsyn_counts{$line[0]}=$n;
		if($line[5] ne ""){
			my @sites=split(/;/, $line[5]);
			$bg_subst_map{$line[0]}={};
			foreach my $str(@sites){
				my @bases;
				my $spos=parse_subst_abbr($str,\@bases);
				if(!defined $bg_site_idx{$spos}){
					$bg_site_idx{$spos}=$I;
					$bg_site2cord[$I++]=$spos;
				};
				$bg_nsubst{$spos}++;
				my $si=SubstInfo->new();
				$si->site($spos);
				$si->bases([@bases]) if @bases;
				$bg_subst_map{$line[0]}->{$spos}=$si;
				if(defined $rh_bg_subst_branches){
					$rh_bg_subst_branches->{$spos}=[] unless defined $rh_bg_subst_branches->{$spos};
					push @{$rh_bg_subst_branches->{$spos}},$line[0];
				}
			};			
		}
		if($npheno){
			if($line[6]=~/^[10]+$/){
				my @pheno=split "", $line[6];
				if(@pheno==$npheno){
					$phenotypes{$line[0]}=[@pheno];
				}else{
					die "\nNumber of labels in the $ARGV[0] doesn't equal to the number of phenotypes on the branch $line[0]!";
				}
			}else{
				die "\nUndefined phenotypes on the branch: $line[0]!";
			}
		}
		$n=0;
		if($line[3] ne ""){
			while($line[3]=~m/;/g){$n++};
			$n++;
		};
		$syn_counts{$line[0]}=$n;
	};
	close INPF;
	#reweight tree
	if($distance_mode){
		$tree->visit_breadth_first(
			-in => sub{
				my $node=shift;
				my $w=0;
				unless($node->is_root||($f_ignore_terminals&&$node->is_terminal)){
					if($distance_mode==1){
						$w=$syn_counts{$node->get_name};
					}elsif($distance_mode==2){
						$w=$nsyn_counts{$node->get_name};
					}elsif($distance_mode==3){
						$w=$syn_counts{$node->get_name}+$nsyn_counts{$node->get_name};
					}
				}
				$node->set_branch_length($w);
			}
		);
	}
	#calculates weights for homoplastic mutations according for branches' confidences
	if(defined $branch_conf_fn){
		my %mweights=ParallelSubstWeights::calculate_mutations_confidences($tree,\%tg_subst_map,\%branch_confidences,$branch_conf_threshold);
		foreach my $name(keys %tg_subst_map){
			foreach my $site(%{$tg_subst_map{$name}}){
				$tg_subst_map{$name}->{$site}->weight($mweights{$name}->{$site}) if defined $mweights{$name}->{$site};
			}
		}
		%mweights=ParallelSubstWeights::calculate_mutations_confidences($tree,\%bg_subst_map,\%branch_confidences,$branch_conf_threshold);
		foreach my $name(keys %bg_subst_map){
			foreach my $site(%{$bg_subst_map{$name}}){
				$bg_subst_map{$name}->{$site}->weight($mweights{$name}->{$site}) if defined $mweights{$name}->{$site};
			}
		}
	}
}

#begin script
if($f_analysis{dist}||$f_analysis{figtree}||$f_analysis{subst_dens}){
	$rh_bg_subst_branches={};
	$rh_tg_subst_branches={};
}
init_subst_map($xpar_fn,$tree,$rh_bg_subst_branches,$rh_tg_subst_branches);
#calculating site pairs
my @site_pair_info;
if($f_stat_type==0||$f_stat_type==1||$f_stat_type==3){
	%site_pairs=Pairs::ConsequtivePairs::find($f_intragene,$tree,$f_ignore_terminals,\%bg_subst_map,\%tg_subst_map,\%bg_site_idx,\%tg_site_idx,\@bg_site2cord);
}elsif($f_stat_type==2){
	%site_pairs=Pairs::BranchPairs::find($f_intragene,$tree,$f_ignore_terminals,\%bg_subst_map,\%tg_subst_map);
}
foreach my $str(keys %site_pairs){
	my $epistat;
	if($f_stat_type==0||$f_stat_type==1||$f_stat_type==3){
		if(defined $branch_conf_fn){
			$epistat=Pairs::ConsequtivePairs::calc_statistics($site_pairs{$str},$f_stat_type,$f_stat_func,$tau,\%bg_subst_map,\%tg_subst_map);
		}else{
			$epistat=Pairs::ConsequtivePairs::calc_statistics($site_pairs{$str},$f_stat_type,$f_stat_func,$tau);
		}
	}elsif($f_stat_type==2){
		if(defined $branch_conf_fn){
			$epistat=Pairs::BranchPairs::calc_statistics($site_pairs{$str},$tau,\%bg_subst_map,\%tg_subst_map);
		}else{
			$epistat=Pairs::BranchPairs::calc_statistics($site_pairs{$str},$tau);
		}
	}
	$site_pairs{$str}->epistat($epistat);
	push @site_pair_info, $site_pairs{$str};
};
%site_pairs=();
#all site pairs initialized
#Printing results
@site_pair_info=sort {
	$a->bg_site<=>$b->bg_site||$a->tg_site<=>$b->tg_site} @site_pair_info;
my @tg_sites_srt=sort {$a<=>$b} @tg_site2cord;
my @bg_sites_srt=sort {$a<=>$b} @bg_site2cord;
my $k=0;
my $spi=$site_pair_info[$k];
if($f_analysis{epistat}){
	if($f_print_pairs){
		open OPF_PAIRS, ">$pairs_fn" or die "\nUnable to open output file $pairs_fn!";
		print OPF_PAIRS "bg site\ttarget site\tbg site subst\ttarget site subst\n";
	};
	open OPF_STAT, ">$stat_fn" or die "\nUnable to open output file $stat_fn!";
	if($f_print_all_pairs){
		for(my $i=0;$i<@bg_sites_srt;$i++){
			my $bgs=$bg_sites_srt[$i];
			for(my $j=0;$j<@tg_sites_srt;$j++){
				my $tgs=$tg_sites_srt[$j];
				next if($f_intragene&&$bgs==$tgs);
				my $bgn=$bg_nsubst{$bgs};
				my $tgn=$tg_nsubst{$tgs};
				print OPF_PAIRS "$bgs\t$tgs\t$bgn\t$tgn\n" if($f_print_pairs);
				if($bgs==$spi->bg_site&&$tgs==$spi->tg_site){
					my $val=$spi->epistat;
					$val=sprintf("%.6f",$val);
					print OPF_STAT "$val\t";
					$val=$spi->npair_subst;
					$val=sprintf("%.2f",$val);
					print OPF_STAT "$val\t";
					($bgn,$tgn)=calc_subst_numbers($spi);
					print OPF_STAT "$bgn\t$tgn\n";
					$k++;
					$spi=$site_pair_info[$k] if $k<@site_pair_info;
				}else{
					print OPF_STAT "0\t0.00\t0\t0\n";
				}
			}
		};
	}else{
		while($k<@site_pair_info){
			my $bgs=$spi->bg_site;
			my $tgs=$spi->tg_site;
			my $bgn=$bg_nsubst{$bgs};
			my $tgn=$tg_nsubst{$tgs};
			print OPF_PAIRS "$bgs\t$tgs\t$bgn\t$tgn\n" if($f_print_pairs);
			my $val=$spi->epistat;
			$val=sprintf("%.6f",$val);
			print OPF_STAT "$val\t";
			$val=$spi->npair_subst;
			$val=sprintf("%.2f",$val);
			print OPF_STAT "$val\t";
			($bgn,$tgn)=calc_subst_numbers($spi);
			print OPF_STAT "$bgn\t$tgn\n";
			$k++;
			$spi=$site_pair_info[$k] if $k<@site_pair_info;
		}
	};
	close OPF_PAIRS if($f_print_pairs);
	close OPF_STAT;
}

sub mk_phenotype_subst_map{
	my ($rh_phenotypes,$nphen,$rh_subst_map,$rao_pheno_subst_map,$rao_pheno_site_sets)=@_;
	@{$rao_pheno_subst_map}=();
	@{$rao_pheno_site_sets}=();
	my @tmp_pheno_site_set;
	for(my $i=0;$i<$nphen;$i++){
		push @{$rao_pheno_subst_map},[({},{})];
		push @tmp_pheno_site_set,[({},{})];
		push @{$rao_pheno_site_sets},[([],[])];
	}
	foreach my $brname(keys %{$rh_subst_map}){
		for(my $i=0;$i<$nphen;$i++){
			my $t=$rh_phenotypes->{$brname}->[$i];
			$rao_pheno_subst_map->[$i]->[$t]->{$brname}=$rh_subst_map->{$brname};
			foreach my $site(keys %{$rh_subst_map->{$brname}}){
				$tmp_pheno_site_set[$i]->[$t]->{$site}++;
			}
		}
	}
	for(my $i=0;$i<$nphen;$i++){
		for(my $j=0;$j<2;$j++){
			@{$rao_pheno_site_sets->[$i]->[$j]}=keys %{$tmp_pheno_site_set[$i]->[$j]};
		}
	}
}

if($f_mtx_mode){
	unless($npheno){
		my @bg_branches=keys %bg_subst_map;
		my @tg_branches=keys %tg_subst_map;
		print_branch_site_matrix($dir.$basename,$f_mtx_mode,$f_intragene,\@bg_site2cord,\@tg_site2cord,
			\@bg_branches,\@tg_branches,\%bg_subst_map,\%tg_subst_map);
	}else{
		my @bg_pheno_subst_map;
		my @bg_pheno_site_set;
		mk_phenotype_subst_map(\%phenotypes,$npheno,\%bg_subst_map,\@bg_pheno_subst_map,\@bg_pheno_site_set);
		my @tg_pheno_subst_map;
		my @tg_pheno_site_set;
		mk_phenotype_subst_map(\%phenotypes,$npheno,\%tg_subst_map,\@tg_pheno_subst_map,\@tg_pheno_site_set);
		for(my $i=0;$i<$npheno;$i++){
			my $ofn_pref=$dir.$basename.".".$pheno_labels[$i];
			for(my $j=0;$j<2;$j++){
				my @bg_branches=keys %{$bg_pheno_subst_map[$i]->[$j]};
				my @tg_branches=keys %{$tg_pheno_subst_map[$i]->[$j]};
				print_branch_site_matrix($ofn_pref."#$j",$f_mtx_mode,$f_intragene,
					$bg_pheno_site_set[$i]->[$j],$tg_pheno_site_set[$i]->[$j],
					\@bg_branches,\@tg_branches,
					$bg_pheno_subst_map[$i]->[$j],$tg_pheno_subst_map[$i]->[$j]);
			}
		}
	}
}

#additional analysis
my @in_site_pairs;
if($input_site_pair_fn){
	open INPF, "<$input_site_pair_fn" or die "\nUnable to open input file $input_site_pair_fn!";
	while(<INPF>){
		chomp;
		if(/\S+/){
			my @line=split ',';
			if($line[0]=~/\d+/ || $line[1]=~/\d+/){
				my $p=[];
				$p->[0]=$line[0] if $line[0]=~/\d+/;
				$p->[1]=$line[1] if $line[1]=~/\d+/;
				push @in_site_pairs,$p;
			}else{
				die "\nError format of input file $input_site_pair_fn!";
			}
		}
	}
	close INPF;
}
#distance analysis
if($f_analysis{dist}){
	my @sp_distance_info=mk_distance_analysis($tree,\@site_pair_info,$rh_bg_subst_branches,$rh_tg_subst_branches);
	#print results
	open OPF_STAT, ">$dist_fn" or die "\nUnable to open output file: $dist_fn!";
	print OPF_STAT "bgr_fgr\tbgr_root\tfgr_root_indicator\tbgr_free_root\tbgr_free_root_indicator\tfgr_free_root\tfgr_free_root_indicator\n";
	$k=0;
	$spi=$site_pair_info[$k] if $k<@site_pair_info;
	my $sp_di=$sp_distance_info[$k] if $k<@site_pair_info;
	foreach my $bgs (@bg_sites_srt){
		foreach my $tgs (@tg_sites_srt){
			next if($f_intragene&&$bgs==$tgs);
			if($bgs==$spi->bg_site&&$tgs==$spi->tg_site){
				my $d=sprintf("%.6f",$sp_di->[0]);
				my $root_dist=sprintf("%.6f",$sp_di->[1]->root());
				my $root_ind=sprintf("%.6f",$sp_di->[1]->root_indicator());
				print OPF_STAT "$d\t$root_dist\t$root_ind\t";
				$root_dist=sprintf("%.6f",$sp_di->[2]->root());
				$root_ind="";
				$root_ind=sprintf("%.6f",$sp_di->[2]->root_indicator()) if defined $sp_di->[2]->root_indicator();
				print OPF_STAT "$root_dist\t$root_ind\t";
				$root_dist=sprintf("%.6f",$sp_di->[3]->root());
				$root_ind="";
				$root_ind=sprintf("%.6f",$sp_di->[3]->root_indicator()) if defined $sp_di->[3]->root_indicator();
				print OPF_STAT "$root_dist\t$root_ind\n";
				$k++;
				$spi=$site_pair_info[$k] if $k<@site_pair_info;
				$sp_di=$sp_distance_info[$k] if $k<@site_pair_info;
			}else{
				print OPF_STAT "0\t0\t0\t0\t0\n";
			}
		}
	}
	close OPF_STAT;
}

if($f_analysis{figtree}){
	my @fig_site_pairs;
	my @site_pos;
	for(my $i=0;$i<@in_site_pairs;$i++){
		if(defined($in_site_pairs[$i]->[0])&&defined($in_site_pairs[$i]->[1])){
			my $key=SitePairInfo->new();
			$key->bg_site($in_site_pairs[$i]->[0]);
			$key->tg_site($in_site_pairs[$i]->[1]);
			my $j=binsearch {$a->bg_site<=>$b->bg_site||$a->tg_site<=>$b->tg_site} $key, @site_pair_info;
			push @fig_site_pairs, $site_pair_info[$j] if(defined $j);
		}else{
			push @site_pos,$in_site_pairs[$i];
		}
	}
	my ($basename,$dir,$ext) = fileparse($input_site_pair_fn,'\.[^\.]*$');
	my $out_figtree_fn=$dir.$basename.".site_pairs.tre";
	open OPF, ">$out_figtree_fn" or die "\nUnable to open output file $out_figtree_fn!";
	print OPF "#NEXUS\nbegin trees;";
	if(@site_pos){
		print OPF "\n\ttree tree_sites = [&R] ";
		my $tree_str;
		if(defined $branch_conf_fn){
			$tree_str=sites2figtree_str(\@site_pos,$rh_bg_subst_branches,$rh_tg_subst_branches,\%bg_subst_map,\%tg_subst_map,$tree,\%branch_confidences);
		}else{
			$tree_str=sites2figtree_str(\@site_pos,$rh_bg_subst_branches,$rh_tg_subst_branches,\%bg_subst_map,\%tg_subst_map,$tree);
		}
		print OPF $tree_str;
	}
	if(@fig_site_pairs){
		print OPF "\n\ttree tree_pairs = [&R] ";
		my $tree_str;
		if(defined $branch_conf_fn){
			$tree_str=site_pairs2figtree_str(\@fig_site_pairs,$rh_bg_subst_branches,$rh_tg_subst_branches,\%bg_subst_map,\%tg_subst_map,$tree,\%branch_confidences);
		}else{
			$tree_str=site_pairs2figtree_str(\@fig_site_pairs,$rh_bg_subst_branches,$rh_tg_subst_branches,\%bg_subst_map,\%tg_subst_map,$tree);
		}
		print OPF $tree_str;
	}
	print OPF "\nend;";
	close OPF;
}

if($f_analysis{subst_dens}){
	my @bg_site_pos;
	my @fg_site_pos;
	for(my $i=0;$i<@in_site_pairs;$i++){
		if(defined $in_site_pairs[$i]->[0]){
			push @bg_site_pos,$in_site_pairs[$i]->[0];
		}
		if(defined $in_site_pairs[$i]->[1]){
			push @fg_site_pos,$in_site_pairs[$i]->[1];
		}
	}
	my ($basename,$dir,$ext) = fileparse($input_site_pair_fn,'\.[^\.]*$');
	my $out_figtree_fn=$dir.$basename.".subst_dens.tre";
	open OPF, ">$out_figtree_fn" or die "\nUnable to open output file $out_figtree_fn!";
	print OPF "#NEXUS\nbegin trees;";
	if(@bg_site_pos){
		print OPF "\n\ttree tree_bg_sites = [&R] ";
		my $tree_str=subst_distrib2figtree_str(\@bg_site_pos,\%bg_subst_map,$tree);
		print OPF $tree_str;
	}
	if(@fg_site_pos){
		print OPF "\n\ttree tree_fg_sites = [&R] ";
		my $tree_str=subst_distrib2figtree_str(\@fg_site_pos,\%tg_subst_map,$tree);
		print OPF $tree_str;
	}
	print OPF "\nend;";
	close OPF;
}
