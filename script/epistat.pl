#!/usr/bin/env perl
#This script calculates epistatic statistics [Kryazhimsky11]
#options:
#	-x <FN> - input file for EpiStat application in XPAR format [Kryazhimsky11]. Contains tree in the adjacent file format
#	[-m <"0|1|2|3">]  PrintMatrixMode - mode for printing of mutation matrix  (see below): 0- don't print;1 - background;2 - target; 3 - both
#			Default="0"
#	[-a <"epistat"+"figtree">] Separator='+', Run analysis:
#			Default="epistat"
#		epistat - calculate epistatic statistics
#		figtree_site_pairs - draw mutations on branches for specified site pairs in FigTree format. This option requires site pair file specified by '-s' option
#		figtree_subst_dens - colour branches on the tree according to differences of observed and expected (assuming uniformity) numbers of substitutions in their subtrees.
#	[-s <site_pair_fn>] - Specifies a file (CSV) with a list of site pairs those mutations will be drawn on the tree
#	[-p] - do not print site pairs
#	[-A] - Specifies a treating mode of ambiguous alleles in abbreviations of mutations in the XPARR file:
#		removes mutations with ambiguous alleles if '-B' and '-F' options are not specified
#		infers ambiguous alleles from site profiles if profile files for the background and for foreground are specified in the '-B' and '-F' options.
#		For the intragenic site pairs it is sufficient to specify one of two options: '-B' or '-F'.
#		For the intergenic site pairs the both options are required
#	[-B <FN>] - specifies a file in the IQTree *.sitefreq format with an allele profile for the background
#	[-F <FN>] - specifies a file in the IQTree *.sitefreq format with an allele profile for the foreground
#Params:
#<parameters> - The file with script params.
#		Contents:
#		[TAU="<ufloat>"] - timescale parameter (>0)
#		PairsType=("intra[gene]"|"inter[gene]")
#		[Alphabet="STR"] - the string of ordered one letter symbols for allele states in mutation codes: <state1><site><state2>.
#			This parameter is requred for calculating of allele statistics
#		[StatType=("all_pairs"|"all_allele_pairs"|"independent_pairs"|"indep_allele_pairs"|"branch_pairs"|"allele_branch_pairs")] - switches the way of convolution into a site pair epistatic statistics of epistatic statistics of corresponding mutation pairs
#			"all_pairs" - all consecutive  mutation pairs are treated independently of each other
#			"all_allele_pairs_A1" - all pairs of derived alleles are considered. The likelihood of a consecutive  mutationn pair is proportional to P(s1,A1|b2,s2,B2)
#				!!!In theory, _A1 is preffered on the _B2 for allele alphabets of small sizes. _B2 requires at least three alleles in a site.
#			"all_allele_pairs_B2" - all pairs of derived alleles are considered. The likelihood of a consecutive  mutationn pair is proportional to 1/2*(P(B2|s1,A1,b2,s2)+P(!B2|!(s1,A1),b2,s2))
#			"independent_pairs" - average epistatic statistics for mutation pairs having common background mutation
#			"indep_allele_pairs_A1" - average epistatic statistics for mutation pairs having common background mutation, all allele pairs are considered. The likelihood of a consecutive  mutationn pair is proportional to P(s1,A1|b2,s2,B2)
#				!!!In theory, _A1 is preffered on the _B2 for allele alphabets of small sizes. _B2 requires at least three alleles in a site.
#			"indep_allele_pairs_B2" - average epistatic statistics for mutation pairs having common background mutation, all allele pairs are considered. The likelihood of a consecutive  mutationn pair is proportional to 1/2*(P(B2|s1,A1,b2,s2)+P(!B2|!(s1,A1),b2,s2))
#			"branch_pairs" - only events occurred on same branches are counted. This option is adjusted for searching for genotype-phenotype associations (GWAS)
#			"allele_branch_pairs" - pairs of alleles derived on same branches are counted. This option is adjusted for searching for genotype-phenotype associations (GWAS)
#		[IgnoreBranchPairs="0|1"] - specifies whether to ignore site pairs on the same branches
#			Default=0 - count all site pairs. The value=1 is not sutable for StatType="branch_pairs" or StatType="allele_branch_pairs"
#		[IgnoreTerminals="0|1"] - specifies whether to ignore terminal branches
#			Default="1"
#		[OutSitePairs="Suffix"] - list of site pairs for which epistatic statistics is calculated
#			Default=".site_pairs"
#		[OutEpiStatistics="Suffix"] - file with epistatic statistics for pairs
#			Default=".stat.exp".$TAU
#		[BranchDistUnits="(DEFAULT|SYN|NSYN|TOTAL)"] - Units for branch length distances
#			Default="SYN"
#		[BranchConfidences="<FN>"] - A file containing confidences for tree branches
#		[BranchConfidenceThreshold="<ufloat>"] - A threshold to accept of a branch support
#		[Phenotypes="1|0"] - Use|Ignore phenotypes assigned to the tree branches

use strict;
#server inst
use lib "$ENV{EPISTAT_LIB}";
#######
use Bio::Phylo::IO;
use Getopt::Std;
use File::Basename;
use EpistatAnalysis::PhyloDrawSites qw(site_pairs2figtree_str sites2figtree_str subst_distrib2figtree_str);
use List::BinarySearch qw(binsearch);
use SitePairInfo;
use Pairs::ConsecutivePairs;
use Pairs::AlleleConsecutivePairs;
use Pairs::BranchPairs;
use Pairs::AlleleBranchPairs;
use ParallelSubstWeights;
use IO::EmpiricalMutationModel qw(print_branch_site_mutation_matrix);
use SiteModel::ProteinMixtureModelProfile;
use IO::XPARR qw(parse_xparr);

#constants
use constant {
	ALL_PAIRS => Pairs::ConsecutivePairs::ALL_PAIRS,
	INDEP_PAIRS => Pairs::ConsecutivePairs::INDEP_PAIRS,
	BRANCH_PAIRS => Pairs::BranchPairs::BRANCH_PAIRS,
	ALL_ALLELE_PAIRS_A1 => Pairs::AlleleConsecutivePairs::ALL_ALLELE_PAIRS_A1,
	ALL_ALLELE_PAIRS_B2 => Pairs::AlleleConsecutivePairs::ALL_ALLELE_PAIRS_B2,
	INDEP_ALLELE_PAIRS_A1 => Pairs::AlleleConsecutivePairs::INDEP_ALLELE_PAIRS_A1,
	INDEP_ALLELE_PAIRS_B2 => Pairs::AlleleConsecutivePairs::INDEP_ALLELE_PAIRS_B2,
	ALLELE_BRANCH_PAIRS => Pairs::AlleleBranchPairs::ALLELE_BRANCH_PAIRS,
};

my %args;
if(!getopts('x:m:a:ps:AB:F:',\%args)){
	die "\nError in option string!";
}

my $xpar_fn=$args{x};
die "\nThe input XPAR file is not specified! Use -x option!" unless defined $xpar_fn;
my $branch_conf_fn;
my %branch_confidences;
my $branch_conf_threshold;
my $stat_fn;
my $pairs_fn=".site_pairs";
my $tau;
my $f_intragene;
my $f_ignore_terminals=1;
my $f_branch_pairs;
my $f_mtx_mode=0;
my $f_print_pairs=1;
my $f_stat_type;
my %f_analysis;
my $dist_fn=".analysis.dist";
my $input_site_pair_fn;
my $npheno;
my @pheno_labels;
my $alphabet_str;
my $ra_alphabet;
my $rh_allele2idx;
my $f_treat_allleles=0;
my $bgr_treat_prof;
my $fgr_treat_prof;

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
			die "\nUnknown value ($value) for the IgnoreTerminals parameter in file $ARGV[0]!\n\tUINT in the range is [0..1] expected!" unless $value=~/^[01]\s*$/;
		}elsif($key eq "PairsType"){
			if($value=~m/^intra(gene)*/i){
				$f_intragene=1;
			}elsif($value=~m/^inter(gene)*/i){
				$f_intragene=0;
			}else{
				die "\nUnknown value ($value) for the PairsType parameter in file $ARGV[0]!"; 
			}
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
				die "\nUnknown value for 'BranchDistUnits'!";
			}
		}elsif($key eq "StatType"){
			if($value=~m/^all_pairs/i){
				$f_stat_type=ALL_PAIRS;
			}elsif($value=~m/^independent_pairs/i){
				$f_stat_type=INDEP_PAIRS;
			}elsif($value=~m/^branch_pairs/i){
				$f_stat_type=BRANCH_PAIRS;
			}elsif($value=~m/^all_allele_pairs_A1/i){
				$f_stat_type=ALL_ALLELE_PAIRS_A1;
			}elsif($value=~m/^all_allele_pairs_B2/i){
				$f_stat_type=ALL_ALLELE_PAIRS_B2;
			}elsif($value=~m/^indep_allele_pairs_A1/i){
				$f_stat_type=INDEP_ALLELE_PAIRS_A1;
			}elsif($value=~m/^indep_allele_pairs_B2/i){
				$f_stat_type=INDEP_ALLELE_PAIRS_B2;
			}elsif($value=~m/^allele_branch_pairs/i){
				$f_stat_type=ALLELE_BRANCH_PAIRS;
			}else{
				die "\nUndefined type of statistics: StatType=$value!";
			}
		}elsif($key eq "IgnoreBranchPairs"){
			if($value=~/^[01]$/){
				if($value==0){
					$f_branch_pairs=1;
				}else{
					$f_branch_pairs=0;
				}
			}
		}elsif($key eq "Alphabet"){
			$alphabet_str=uc $value;
			$alphabet_str=~s/\s//;
			$ra_alphabet=[];
			@{$ra_alphabet}=split "",$alphabet_str;
			$rh_allele2idx={};
			for(my $i=0;$i<@{$ra_alphabet};$i++){
				if(defined $rh_allele2idx->{$ra_alphabet->[$i]}){
					die "\nError: The alphabet string contains the character ".$ra_alphabet->[$i]." several times!";
				}else{
					$rh_allele2idx->{$ra_alphabet->[$i]}=$i;
				}
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

die "\nThe calculating statistics type is undefined! Please, specify the StatType parameter!" unless defined $f_stat_type;
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
die "\nError: The TAU parameter has no defaul value!" unless defined($tau)||($f_stat_type==BRANCH_PAIRS)||($f_stat_type==ALLELE_BRANCH_PAIRS);
$f_branch_pairs=1 unless defined $f_branch_pairs;
if($f_branch_pairs==0){
	die "\nError: Unable to ignore pairs of mutations occurred on the same branches if StatType=\"ALLELE_BRANCH_PAIRS\"|\"ALLELE_BRANCH_PAIRS\"!" if $f_stat_type==ALLELE_BRANCH_PAIRS||
		$f_stat_type==BRANCH_PAIRS;
}
unless(defined $stat_fn){
	$stat_fn=".stat";
	$stat_fn.=".exp.$tau" if defined $tau;
}

if(defined $args{A}){
	$f_treat_allleles=1;
	die "\nError: the alphabet is not specified! Use 'Alphabet' parameter in $ARGV[0]!" unless defined $ra_alphabet;
	my $bgr_pro_fn=$args{B};
	my $fgr_pro_fn=$args{F};
	unless(defined($bgr_pro_fn)||defined($fgr_pro_fn)){
		$f_treat_allleles=2;
	}elsif(defined($bgr_pro_fn)&&defined($fgr_pro_fn)){
		if($f_intragene){
			if($bgr_pro_fn eq $fgr_pro_fn){
				$fgr_pro_fn=undef;
			}else{
				die "\nError: different site profiles for the intragene analysis accounted!";
			}
		}
	}else{
		if($f_intragene){
			if(defined $fgr_pro_fn){
				$bgr_pro_fn=$fgr_pro_fn;
				$fgr_pro_fn=undef;
			}
		}else{
			die "\nError: both site profiles are required to treat ambiguous allleles in the intergenic mode!";
		}
	}
	$bgr_treat_prof=SiteModel::ProteinMixtureModelProfile->new($bgr_pro_fn,$ra_alphabet) if defined $bgr_pro_fn;
	$fgr_treat_prof=SiteModel::ProteinMixtureModelProfile->new($fgr_pro_fn,$ra_alphabet) if defined $fgr_pro_fn;
	$fgr_treat_prof=$bgr_treat_prof if $f_intragene;
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
my $rh_terminals;
if($f_ignore_terminals){
	$rh_terminals={};
	foreach my $node(@{$tree->get_terminals}){
		my $name=$node->get_name;
		$rh_terminals->{$name}=1;
	}
}

#my $str=Bio::Phylo::IO->unparse(-phylo => $tree,-format => 'newick');
#print $str;

my %bg_site2idx; #converts background site coordinate into array index
my @bg_idx2site;
my %tg_site2idx; #converts target site coordinate into array index
my @tg_idx2site;
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

sub calc_subst_numbers{
	my $sp_info=shift;
	my $compl_sp_info=shift;
	my %bg_branches;
	my %tg_branches;
	if(defined $sp_info){
		foreach my $sb_info(@{$sp_info->subst_info}){
			$bg_branches{$sb_info->bg_branch()->get_name()}=1;
			$tg_branches{$sb_info->tg_branch()->get_name()}=1;
		}
	}
	if(defined $compl_sp_info){
		foreach my $sb_info(@{$compl_sp_info->subst_info}){
			#This version calculates total number of different background branches and total number of different foreground branches
			#$bg_branches{$sb_info->bg_branch()->get_name()}=1;
			#$tg_branches{$sb_info->tg_branch()->get_name()}=1;
			#This version calculates total number of different branches in the first site and the total number of different branches in the secons site of the first site pair
			$bg_branches{$sb_info->tg_branch()->get_name()}=1;
			$tg_branches{$sb_info->bg_branch()->get_name()}=1;
		}
	}
	return (scalar(keys %bg_branches),scalar(keys %tg_branches));
}

#begin script
if($f_analysis{dist}||$f_analysis{figtree}||$f_analysis{subst_dens}){
	$rh_bg_subst_branches={};
	$rh_tg_subst_branches={};
}
{
	my $rh_syn_ncounts;
	my $rh_nsyn_ncounts;
	if($distance_mode==1||$distance_mode==3){
		$rh_syn_ncounts={};
	}
	if($distance_mode==2||$distance_mode==3){
		$rh_nsyn_ncounts={};
	}
	parse_xparr($xpar_fn,"-f_treat_alleles" => $f_treat_allleles,"-bgr_treat_profile" => $bgr_treat_prof,"-fgr_treat_profile" => $fgr_treat_prof,
		"-allele_indices" => $rh_allele2idx, "-f_ignore_terminals" => $f_ignore_terminals,"-term_branches" => $rh_terminals,
		"-orh_bgr_site_nsubst" => \%bg_nsubst,"-orh_fgr_site_nsubst" => \%tg_nsubst,
		"-orh_bgr_subst_map" => \%bg_subst_map, "-orh_fgr_subst_map" => \%tg_subst_map,
		"-orh_bgr_site2idx" => \%bg_site2idx, "-orh_fgr_site2idx" => \%tg_site2idx,
		"-ora_bgr_idx2site" => \@bg_idx2site, "-ora_fgr_idx2site" => \@tg_idx2site,
		"-orh_syn_counts" => $rh_syn_ncounts, "-orh_fgr_nsyn_counts" => $rh_nsyn_ncounts,
		"-orh_bgr_subst2branch" => $rh_bg_subst_branches, "-orh_fgr_subst2branch" => $rh_tg_subst_branches,
		"-orh_phenotypes" => \%phenotypes, "-ora_pheno_labels" => \@pheno_labels
	);
	if(defined $npheno){
		$npheno=@pheno_labels;
#DEBUG:
#print "\nBGR:";
#foreach my $site (keys %bg_site2idx){
#	print "\n$site\t$bg_site2idx{$site}\t$tg_site2idx{$site}" unless defined $tg_site2idx{$site};
#}
#print "\n\nFGR:";
#foreach my $site (keys %tg_site2idx){
#	print "\n$site\t$bg_site2idx{$site}\t$tg_site2idx{$site}" unless defined $bg_site2idx{$site};
#}
#exit;
#
	}
	#reweight tree
	if($distance_mode){
		foreach my $node($tree->get_nodes){
			my $name=$node->get_name();
			my $l=0;
			if($distance_mode==1||$distance_mode==3){
				$l+=$rh_syn_ncounts->{$name};
			}
			if($distance_mode==2||$distance_mode==3){
				$l+=$rh_nsyn_ncounts->{$name};
			}
			$node->set_branch_length($l);
		}
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
#calculating site pairs
sub is_allele_pairs_type{
	my $f_stat_type=shift;
	return $f_stat_type==ALL_ALLELE_PAIRS_A1||$f_stat_type==ALL_ALLELE_PAIRS_B2||$f_stat_type==INDEP_ALLELE_PAIRS_A1||$f_stat_type==INDEP_ALLELE_PAIRS_B2||$f_stat_type==ALLELE_BRANCH_PAIRS;
}

my @site_pair_info;
{
	my $sp_calc;
	if($f_stat_type==ALL_PAIRS||$f_stat_type==INDEP_PAIRS){
		$sp_calc=Pairs::ConsecutivePairs->new($f_intragene,$tree,$f_ignore_terminals,
			\%bg_subst_map,\%tg_subst_map,\%bg_site2idx,\%tg_site2idx,\@bg_idx2site,$f_branch_pairs);
		if(defined $branch_conf_fn){
			@site_pair_info=$sp_calc->calc_statistics($tau,$f_stat_type,\%bg_subst_map,\%tg_subst_map);
		}else{
			@site_pair_info=$sp_calc->calc_statistics($tau,$f_stat_type);
		}
	}elsif($f_stat_type==BRANCH_PAIRS){
		$sp_calc=Pairs::BranchPairs->new($f_intragene,$tree,$f_ignore_terminals,
			\%bg_subst_map,\%tg_subst_map,\%bg_site2idx,\%tg_site2idx);
		if(defined $branch_conf_fn){
			@site_pair_info=$sp_calc->calc_statistics($tau,\%bg_subst_map,\%tg_subst_map);
		}else{
			@site_pair_info=$sp_calc->calc_statistics($tau);
		}
	}elsif(is_allele_pairs_type($f_stat_type)){
		unless($f_stat_type==ALLELE_BRANCH_PAIRS){
			$sp_calc=Pairs::AlleleConsecutivePairs->new($f_intragene,$tree,$f_ignore_terminals,
				\%bg_subst_map,\%tg_subst_map,\%bg_site2idx,\%tg_site2idx,\@bg_idx2site,$f_branch_pairs,$ra_alphabet);
		}else{
			$sp_calc=Pairs::AlleleBranchPairs->new($f_intragene,$tree,$f_ignore_terminals,
				\%bg_subst_map,\%tg_subst_map,\%bg_site2idx,\%tg_site2idx,$ra_alphabet);
		}
		@site_pair_info=$sp_calc->calc_statistics($tau,$f_stat_type,\%bg_subst_map,\%tg_subst_map);
		#DEBUG:
		#$sp_calc->print_algn_profile(*STDOUT,2,\%tg_site2idx);
		#exit;
		####
	}
}
#all site pairs initialized
#Printing results
sub print_epistat{
	my ($sp_idx,$spi)=@_;
	print OPF_STAT "$sp_idx\t";
	my $val=$spi->epistat;
	$val=sprintf("%.6f",$val);
	print OPF_STAT "$val\t";
	$val=$spi->npair_subst;
	$val=sprintf("%.2f",$val);
	print OPF_STAT "$val\t";
	my ($bgn,$tgn)=calc_subst_numbers($spi);
	print OPF_STAT "$bgn\t$tgn\n";
}

sub print_allele_epistat{
	my ($sp_idx,$spi1,$spi2)=@_;
	#DEBUG
	die "\nAt least one of ordered site pairs required!" unless defined($spi1)||defined($spi2);
	if(defined($spi1)&&defined($spi2)){
		die "\nError print_allele_epistat(): not complementary pairs!" unless $spi1->bg_site==$spi2->tg_site&&$spi1->tg_site==$spi2->bg_site;
	}
	####
	print OPF_STAT "$sp_idx\t";
	my $val;
	if(defined $spi1){
		$val=$spi1->epistat;
	}else{
		$val=$spi2->epistat;
	}
	$val=sprintf("%.6f",$val);
	print OPF_STAT "$val\t";
	$val=0;
	$val+=$spi1->npair_subst if defined $spi1;
	$val+=$spi2->npair_subst if defined $spi2;
	$val=sprintf("%.2f",$val);
	print OPF_STAT "$val\t";
	my ($bgn,$tgn)=calc_subst_numbers($spi1,$spi2);
	print OPF_STAT "$bgn\t$tgn\n";
}

if($f_analysis{epistat}){
	if($f_print_pairs){
		open OPF_PAIRS, ">$pairs_fn" or die "\nUnable to open output file $pairs_fn!";
		print OPF_PAIRS "bg site\ttarget site\tbg site subst\ttarget site subst\n";
	};
	open OPF_STAT, ">$stat_fn" or die "\nUnable to open output file $stat_fn!";
	{
		print OPF_STAT "NMaxDataLines=\"";
		my $n=scalar(@bg_idx2site)*scalar(@tg_idx2site);
		if($f_intragene){
			if(is_allele_pairs_type($f_stat_type)){
				$n=scalar(@bg_idx2site);
				$n*=($n-1);
				$n/=2;
			}else{
				foreach my $i(@bg_idx2site){
					$n-- if defined $tg_site2idx{$i};
				}
			}
		}
		print OPF_STAT "$n\"";
	}
	print OPF_STAT "\nsp_index\tepistat\tnpair_subst\tbgrn\tfgrn\n";
	my $k=0;
	my $spi=$site_pair_info[$k];
	my $sp_idx=0;
	my %pair2idx;
	for(my $i=0;$i<@bg_idx2site;$i++){
		my $bgs=$bg_idx2site[$i];
		for(my $j=0;$j<@tg_idx2site;$j++){
			my $tgs=$tg_idx2site[$j];
			next if($f_intragene&&$bgs==$tgs);
			my $bgn=$bg_nsubst{$bgs};
			my $tgn=$tg_nsubst{$tgs};
			print OPF_PAIRS "$bgs\t$tgs\t$bgn\t$tgn\n" if($f_print_pairs);
			$sp_idx++;
			if(is_allele_pairs_type($f_stat_type)&&$f_intragene){
				if($bgs<$tgs){
					my $str="$bgs,$tgs";
					$pair2idx{$str}=[($sp_idx)] unless defined $pair2idx{$str};
				}
			}
			if($bgs==$spi->bg_site&&$tgs==$spi->tg_site){
				unless(is_allele_pairs_type($f_stat_type)&&$f_intragene){
					print_epistat($sp_idx,$spi);
				}else{
					my $str;
					my $i;
					if($bgs<$tgs){
						$str="$bgs,$tgs";
						$i=1;
					}else{
						$str="$tgs,$bgs";
						$i=2;
					}
					$pair2idx{$str}->[$i]=$k;
				}		
				$k++;
				$spi=$site_pair_info[$k] if $k<@site_pair_info;
			}
		}
	}
	if(is_allele_pairs_type($f_stat_type)&&$f_intragene){
		my @tmp;
		foreach my $str(keys %pair2idx){
			push @tmp,$pair2idx{$str} if defined($pair2idx{$str}->[1])||defined($pair2idx{$str}->[2]);
		}
		foreach my $si(sort {$a->[0]<=>$b->[0]} @tmp){
			my $spi1=$site_pair_info[$si->[1]] if defined $si->[1];
			my $spi2=$site_pair_info[$si->[2]] if defined $si->[2];
			print_allele_epistat($si->[0],$spi1,$spi2);
		}
		%pair2idx=();
	}
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
		print_branch_site_mutation_matrix($dir.$basename,$f_mtx_mode,$f_intragene,\@bg_idx2site,\@tg_idx2site,
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
				print_branch_site_mutation_matrix($ofn_pref."#$j",$f_mtx_mode,$f_intragene,
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

if($f_analysis{figtree}){
	my @fig_site_pairs;
	my @site_pos;
	for(my $i=0;$i<@in_site_pairs;$i++){
		if(defined($in_site_pairs[$i]->[0])&&defined($in_site_pairs[$i]->[1])){
			my $key=SitePairInfo->new();
			$key->bg_site($in_site_pairs[$i]->[0]);
			$key->tg_site($in_site_pairs[$i]->[1]);
			my $j=binsearch {$a->bg_site<=>$b->bg_site||$a->tg_site<=>$b->tg_site} $key, @site_pair_info;
			if(defined $j){
				push @fig_site_pairs, $site_pair_info[$j];
			}else{
				push @fig_site_pairs, $key;
			}
		}else{
			push @site_pos,$in_site_pairs[$i];
		}
	}
	my ($basename,$dir,$ext) = fileparse($input_site_pair_fn,'\.[^\.]*$');
	my $out_figtree_fn=$dir.$basename.".site_pairs.tre";
	open OPF, ">$out_figtree_fn" or die "\nUnable to open output file $out_figtree_fn!";
	print OPF "#NEXUS\nbegin trees;";
	my $f_prune_tree=$tree->get_ntax()<1000?0:1;
	if(@site_pos){
		print OPF "\n\ttree tree_sites = [&R] ";
		my $tree_str;
		if(defined $branch_conf_fn){
			$tree_str=sites2figtree_str(\@site_pos,$rh_bg_subst_branches,$rh_tg_subst_branches,\%bg_subst_map,\%tg_subst_map,$tree,
				'-ra_idx2allele'=>$ra_alphabet,'-rh_branch_confidences'=>\%branch_confidences,'-f_prune_tree'=>$f_prune_tree);
		}else{
			$tree_str=sites2figtree_str(\@site_pos,$rh_bg_subst_branches,$rh_tg_subst_branches,\%bg_subst_map,\%tg_subst_map,$tree,
				'-ra_idx2allele'=>$ra_alphabet,'-f_prune_tree'=>$f_prune_tree);
		}
		print OPF $tree_str;
	}
	if(@fig_site_pairs){
		print OPF "\n\ttree tree_pairs = [&R] ";
		my $tree_str;
		if(defined $branch_conf_fn){
			$tree_str=site_pairs2figtree_str(\@fig_site_pairs,$rh_bg_subst_branches,$rh_tg_subst_branches,\%bg_subst_map,\%tg_subst_map,$tree,
				'-ra_idx2allele'=>$ra_alphabet,'-rh_branch_confidences'=>\%branch_confidences,
				'-rh_phenotypes' => \%phenotypes, '-ra_pheno_labels' => \@pheno_labels,'-f_prune_tree'=>$f_prune_tree);
		}else{
			$tree_str=site_pairs2figtree_str(\@fig_site_pairs,$rh_bg_subst_branches,$rh_tg_subst_branches,\%bg_subst_map,\%tg_subst_map,$tree,
				'-ra_idx2allele'=>$ra_alphabet,
				'-rh_phenotypes' => \%phenotypes, '-ra_pheno_labels' => \@pheno_labels,'-f_prune_tree'=>$f_prune_tree);
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
