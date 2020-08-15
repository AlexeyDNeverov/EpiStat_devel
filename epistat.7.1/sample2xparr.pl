#!/usr/bin/env perl
#This script converts matrix containing distribution of mutations on a tree branches into XPARR format.
#[options:]
#	-b <file name> - Sets BgrMutMatrix parameter
#	-t <file name> - Sets FgrMutMatrix parameter
#Params:
#<parameters> - The file with script params.
#		Contents:
#		XPAR="FN" - input file for EpiStat application [Kryazhimsky11]. Contains tree in the adjacent file format
#		PairsType=("intra[gene]"|"inter[gene]")
#		MatrixMode="1|2|3" - mode of mutation matrix  (see below): 1 - background;2 - target; 3 - both
#		IgnoreTerminals="0|1" - Ignore terminal branches
#			Default=1 - Ignore
#		[NBranchSubsets="1|2"] - Number of subsets of branches for which distributions of mutations were defined separately.
#			!!!Note that for each subset binary matrices of mutations in the background or in the foreground are required
#			Default="1"
#		[SubsetFolder="FOLDER"] - The folder containing data of a subset. If specified used as a path prefix for the corresponding mutations matrix.
#			!!!The number of specification is required to be equal to the number of subsets or to the zero.
#			Default=""
#		[BgrMutMatrix="FN"] - Distribution mutations in background. Binary matrix: raws - branches; columns - sites
#			!!!If defined, the number of specification is required to be equal to the number of subsets or to the one.
#			!!!In the case of a single specification the matrix with the same file name is expected in each subset folder.
#			Default=undef
#			!!!Option -b has a priority
#		BgrIdx2Site="FN" - Distribution mutations in background. Table to convert columns' indexes of MutMatrix into site positions
#			!!!If defined, the number of specifications is required to be equal to the number of subsets
#			Default=undef
#		BgrIdx2Branch="FN" - Distribution mutations in background. Table to convert raws' indexes of MutMatrix into branch names
#			!!!If defined, the number of specifications is required to be equal to the number of subsets
#			Default=undef
#		[FgrMutMatrix="FN"] - Distribution mutations in foreground. Binary matrix: raws - branches; columns - sites
#			!!!If defined, the number of specifications is required to be equal to the number of subsets or to the one.
#			!!!In the case of a single specification the matrix with the same file name is expected in each subset folder.
#			Default=undef
#			!!!Option -t has a priority
#		FgrIdx2Site="FN" - Distribution mutations in foreground. Table to convert columns' indexes of MutMatrix into site positions
#			!!!If defined, the number of specifications is required to be equal to the number of subsets
#			Default=undef
#		FgrIdx2Branch="FN" - Distribution mutations in foreground. Table to convert raws' indexes of MutMatrix in
#			!!!If defined, the number of specifications is required to be equal to the number of subsets
#			Default=undef

use strict;
use Bio::Phylo::IO;
use Class::Struct;
use Getopt::Std;
use File::Basename;

my %args;
if(!getopts('b:t:',\%args)){
	die "\nError in option string!";
}

my $xpar_fn;
my @bg_mmtx_fn;
push @bg_mmtx_fn,split(',',$args{b}) if defined $args{b};
my @tg_mmtx_fn;
push @tg_mmtx_fn,split(',',$args{t}) if defined $args{t};
my @bg_idx2site_fn;
my @bg_idx2branch_fn;
my @tg_idx2site_fn;
my @tg_idx2branch_fn;
my $f_intragene;
my $f_mtx_mode;
my $f_ignore_terminals=1;
#####################
my $brsubsets_num=1;
my @subset_folders;

open INFILE, "$ARGV[0]" or die "\nUnable to open input file: $ARGV[0]!";
while(<INFILE>){
	if(/#/){
		my $str1=$`;
		my $str2=$';
		$_=$str1 unless($str1=~/\"[^\"]*$/&&$str2=~/^[^\"]*\"/);
	}
	if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
		my ($key,$value)=($1,$2);
		if($key eq "XPAR"){
			$xpar_fn=$value;
		}elsif($key eq "MatrixMode"){
			$f_mtx_mode=$value;
			die "\nUnknown value ($value) for the PrintMatrixMode parameter in file $ARGV[0]!\n\tUINT in range is [1..3] expected!" unless $value=~/^[123]\s*$/;
		}elsif($key eq "IgnoreTerminals"){
			$f_ignore_terminals=$value;
			die "\nUnknown value ($value) for the IgnoreTerminals parameter in file $ARGV[0]!\n\tUINT in range is [0..1] expected!" unless $value=~/^[01]\s*$/;
		}elsif($key eq "BgrMutMatrix"){
			$value=~s/^\s+//;
			push @bg_mmtx_fn,$value unless defined $args{b};
		}elsif($key eq "BgrIdx2Site"){
			push @bg_idx2site_fn,$value;
		}elsif($key eq "BgrIdx2Branch"){
			push @bg_idx2branch_fn,$value;
		}elsif($key eq "FgrMutMatrix"){
			$value=~s/^\s+//;
			push @tg_mmtx_fn,$value unless defined $args{t};
		}elsif($key eq "FgrIdx2Site"){
			push @tg_idx2site_fn,$value;
		}elsif($key eq "FgrIdx2Branch"){
			push @tg_idx2branch_fn,$value;
		}elsif($key eq "PairsType"){
			if($value=~m/^intra(gene)*/i){
				$f_intragene=1;
			}elsif($value=~m/^inter(gene)*/i){
				$f_intragene=0;
			}else{
				die "\nUnknown value ($value) for the PairsType parameter in file $ARGV[0]!"; 
			}
		}elsif($key eq "NBranchSubsets"){
			$value=~s/^\s+//;
			if($value=~/^\d+\s*/){
				$brsubsets_num=$value;
			}else{
				die "\nUnpropper value for the 'NBranchSubsets' parameter: $value! UINT is expected.";
			}
		}elsif($key eq "SubsetFolder"){
			$value=~s/\s+$//;
			$value=~s/\/$//;
			$value.="/";
			push @subset_folders,$value;
		}else{
			die "\nWrong parameter: $key!";
		}
	}
}
close INFILE;
die "\nError: Parameter PairsType is undefined in the file $ARGV[0]!" unless defined $f_intragene;
die "\nThe number of 'BgrIdx2Site' records must be equal to the number of branch subsets: $brsubsets_num!" unless @bg_idx2site_fn==0||$brsubsets_num==@bg_idx2site_fn;
die "\nThe number of 'BgrIdx2Branch' records must be equal to the number of branch subsets: $brsubsets_num!" unless @bg_idx2branch_fn==0||$brsubsets_num==@bg_idx2branch_fn;
die "\nThe number of 'FgrIdx2Site' records must be equal to the number of branch subsets: $brsubsets_num!" unless @tg_idx2site_fn==0||$brsubsets_num==@tg_idx2site_fn;
die "\nThe number of 'FgrIdx2Branch' records must be equal to the number of branch subsets: $brsubsets_num!" unless @tg_idx2branch_fn==0||$brsubsets_num==@tg_idx2branch_fn;
die "\nUndefined where has the mutations been perturbed in background or foreground?" unless @bg_mmtx_fn||@tg_mmtx_fn;
die "\nThe number of 'SubsetFolder' records must be equal to the number of branch subsets: $brsubsets_num!" unless @subset_folders==0||$brsubsets_num==@subset_folders;
if(@bg_mmtx_fn){
	if(@bg_mmtx_fn==1){
		for(my $i=1;$i<$brsubsets_num;$i++){
			push @bg_mmtx_fn,$bg_mmtx_fn[0];
		}
	}elsif(@bg_mmtx_fn!=$brsubsets_num){
		die "\nThe number of background mutation matrices must be equal to the number of branch subsets: $brsubsets_num!";
	}
}
if(@tg_mmtx_fn){
	if(@tg_mmtx_fn==1){
		for(my $i=1;$i<$brsubsets_num;$i++){
			push @tg_mmtx_fn,$tg_mmtx_fn[0];
		}
	}elsif(@tg_mmtx_fn!=$brsubsets_num){
		die "\nThe number of foreground mutation matrices must be equal to the number of branch subsets: $brsubsets_num!";
	}
}
if(@subset_folders){
	die "\nThe number of 'SubsetFolder' records must be equal to the number of branch subsets: $brsubsets_num!" unless $brsubsets_num==@subset_folders;
	for(my $i=0;$i<$brsubsets_num;$i++){
		$bg_mmtx_fn[$i]=$subset_folders[$i].$bg_mmtx_fn[$i] if @bg_mmtx_fn;
		$tg_mmtx_fn[$i]=$subset_folders[$i].$tg_mmtx_fn[$i] if @tg_mmtx_fn;
	}
}

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
	}
}

#my $str=Bio::Phylo::IO->unparse(-phylo => $tree,-format => 'newick');
#print $str;

#substitution counts
my %bg_nsubst;
my %tg_nsubst;
 
#Original distribution of substitutions on branches
my %tg_subst_map;
my %bg_subst_map;
my %scaling_subst_map; #synonimous substitutions
my %phenotype_map;
my $f_phenotypes=0;

#New distribution of substitutions on branches
my $hr_smpl_bg_subst;
my $hr_smpl_tg_subst;
my $hr_smpl_bg_subst_map;
my $hr_smpl_tg_subst_map;

sub parse_subst_abbr{
	my $str=shift;
	my $ra_bases=shift;
	my $spos;
	if($str=~/([A-Za-z])(\d+)([[A-Za-z])/){
		@{$ra_bases}=($1,$3) if(defined $ra_bases);
		$spos=$2;
	}else{
		$str=~s/\D//g;
		$spos=$str;
	};
	return $spos;
}
	
sub init_subst_map{
	my $xpar_fn=shift;
	open INPF, "<$xpar_fn" or die "\nUnable to open input file $xpar_fn!";
	$_=<INPF>; #skip header line
	die "\nHeader line was omitted in file: $xpar_fn!" unless /^child\tparent\tlength/;
	$_=<INPF>; #skip root node statement
	die "\nThe root node statement required in file: $xpar_fn!" unless /^[\w.]+\s*$/;
	while(<INPF>){
		chomp;
		my @line=split "\t";
		if($f_ignore_terminals&&$terminals{$line[0]}){
			$terminals{$line[0]}=$_;
		};
		if($line[4] ne ""){
			my @sites=split(/;/, $line[4]);
			$tg_subst_map{$line[0]}={};
			foreach my $str(@sites){
				my @bases;
				my $spos=parse_subst_abbr($str,\@bases);
				$tg_nsubst{$spos}++;
				$tg_subst_map{$line[0]}->{$spos}=1;
			};
		};
		if($line[5] ne ""){
			my @sites=split(/;/, $line[5]);
			$bg_subst_map{$line[0]}={};
			foreach my $str(@sites){
				my @bases;
				my $spos=parse_subst_abbr($str,\@bases);
				$bg_nsubst{$spos}++;
				$bg_subst_map{$line[0]}->{$spos}=1;
			};			
		};
		$scaling_subst_map{$line[0]}=$line[3];
		if(defined($line[6])&&($line[6] ne "")){
			$f_phenotypes=1 unless $f_phenotypes;
			my $str=join '\t',@line[6..$#line];
			$phenotype_map{$line[0]}=$str;
		}
	};
	close INPF;
}

sub init_subst_distribution{
	my ($nsubsets,$ra_mmtx_fn,$ra_branch_fn,$ra_site_fn,$rh_map,$rh_site)=@_;
	%{$rh_site}=();
	%{$rh_map}=();
	for(my $k=0;$k<$nsubsets;$k++){
		my $mmtx_fn=$ra_mmtx_fn->[$k];
		my $branch_fn=$ra_branch_fn->[$k];
		my $site_fn=$ra_site_fn->[$k];
		die "\nError init_mutation_distribution(): Undefined parameter(s)!" unless defined($mmtx_fn)&&defined($branch_fn)&&defined($site_fn);
		
		my @idx2site;
		my @idx2branch;
		my $i=0;
		open INPF, "<$branch_fn" or die "\nUnable to open input file: $branch_fn!";
		while(<INPF>){
			chomp;
			if($_ ne ""){
				my @line=split "\t";
				$idx2branch[$i++]=$line[0];
			};
		};
		close INPF;
		$i=0;
		open INPF, "<$site_fn" or die "\nUnable to open input file: $site_fn!";
		while(<INPF>){
			chomp;
			if($_ ne ""){
				my @line=split "\t";
				$rh_site->{$line[0]}=0 unless defined $rh_site->{$line[0]};
				$idx2site[$i++]=$line[0];
			};
		};
		close INPF;
		open INPF, "<$mmtx_fn" or die "\nUnable to open input file: $mmtx_fn!";
		$i=0;
		#DEBUG
		#while(<INPF>){
		#	chomp;
		#	my $bname=$idx2branch[$i++];
		#	next if $_ eq "";
		#	$rh_map->{$bname}={};
		#	while(/(1+)/g){
		#		my $idx=pos;
		#		my $n=length $1;
		#		for(my $j=$idx-$n;$j<$idx;$j++){
		#			my $spos=$idx2site[$j];
		#			$rh_site->{$spos}++;
		#			$rh_map->{$bname}->{$spos}=1;
		#		}
		#	}
		#}
		#END DEBUG			
		while(<INPF>){
			$_=~s/\s+$//;
			chomp;
			my $bname=$idx2branch[$i++];
			next if $_ eq "";
			$rh_map->{$bname}={};
			my @sites=split(/;/, $_);
			foreach my $j(@sites){
				my $spos=$idx2site[$j-1];
				$rh_site->{$spos}++;
				$rh_map->{$bname}->{$spos}=1;
			}
		}
		die "\nThe number of lines in the $mmtx_fn\n\tisn't equal to the number of branches in the $branch_fn!" unless $i==@idx2branch;
		close INPF;
	}
}
	
#begin script
init_subst_map($xpar_fn);
my $f_bgr_eq_fgr=0;
if($f_intragene&&$f_mtx_mode==3){
	if(@bg_mmtx_fn==0){
		@bg_mmtx_fn=@tg_mmtx_fn;
		@bg_idx2site_fn=@tg_idx2site_fn;
		@bg_idx2branch_fn=@tg_idx2branch_fn;
	}elsif(@tg_mmtx_fn==0){
		$f_bgr_eq_fgr=1;
		@tg_mmtx_fn=@bg_mmtx_fn;
		@tg_idx2site_fn=@bg_idx2site_fn;
		@tg_idx2branch_fn=@bg_idx2branch_fn;
	}
}
if($f_mtx_mode!=2){
	die "Error: one or more input files for background are undefined!" unless @bg_mmtx_fn&&@bg_idx2branch_fn&&@bg_idx2site_fn;
	$hr_smpl_bg_subst_map={};
	$hr_smpl_bg_subst={};
	init_subst_distribution($brsubsets_num,\@bg_mmtx_fn,\@bg_idx2branch_fn,\@bg_idx2site_fn,$hr_smpl_bg_subst_map,$hr_smpl_bg_subst);
}
if($f_mtx_mode!=1){
	die "Error: one or more input files for foreground are undefined!" unless @tg_mmtx_fn&&@tg_idx2branch_fn&&@tg_idx2site_fn;
	unless($f_bgr_eq_fgr){
		$hr_smpl_tg_subst_map={};
		$hr_smpl_tg_subst={};
		init_subst_distribution($brsubsets_num,\@tg_mmtx_fn,\@tg_idx2branch_fn,\@tg_idx2site_fn,$hr_smpl_tg_subst_map,$hr_smpl_tg_subst);
	}
}
if($f_intragene&&defined($hr_smpl_tg_subst)&&defined($hr_smpl_bg_subst)){
	#check that permutated subsets of sites don't overlap
	foreach my $site(keys %{$hr_smpl_tg_subst}){
		if(exists $hr_smpl_bg_subst->{$site}){
			die "\nError: Independent permutations of both background and foreground overlapped subsets of sites are impossible!";
		}
	}
}
#print XPARR
print "child\tparent\tlength";
foreach my $node($tree->get_nodes){
	my $name=$node->get_name;
	if($f_ignore_terminals&&$node->is_terminal){
		print "\n$terminals{$name}";
		next;
	}
	print "\n$name";
	if(!$node->is_root){
		my $pname=$node->get_parent->get_name;
		my $length=$node->get_branch_length;
		print "\t$pname\t$length\t$scaling_subst_map{$name}\t";
		my @sites;
		my @tmp;
		if(defined $tg_subst_map{$name}){
			if(defined $hr_smpl_tg_subst_map){
				@tmp=keys %{$hr_smpl_tg_subst_map->{$name}};
			}else{
				@tmp=keys %{$tg_subst_map{$name}};
			}
			if(!$f_intragene || defined($hr_smpl_tg_subst_map)){
				@sites=@tmp;
			}else{
				@sites=();
				for(my $i=0;$i<@tmp;$i++){
					push @sites,$tmp[$i] unless defined $bg_nsubst{$tmp[$i]};
				}
				@tmp=keys %{$hr_smpl_bg_subst_map->{$name}};
				for(my $i=0;$i<@tmp;$i++){
					push @sites,$tmp[$i] if(defined $tg_nsubst{$tmp[$i]});
				}
			}
			die "\nError: no mutations in foreground on the branch: $name!" unless scalar @sites;
			@sites=sort {$a<=>$b} @sites;
			for(my $i=0;$i<@sites;$i++){
				print $sites[$i];
				print ";" unless $i==@sites-1;
			}
		};
		print "\t";
		if(defined $bg_subst_map{$name}){
			if(defined $hr_smpl_bg_subst_map){
				@tmp=keys %{$hr_smpl_bg_subst_map->{$name}};
			}else{
				@tmp=keys %{$bg_subst_map{$name}};
			}
			if(!$f_intragene || defined($hr_smpl_bg_subst_map)){
				@sites=@tmp;
			}else{
				@sites=();
				for(my $i=0;$i<@tmp;$i++){
					push @sites,$tmp[$i] if(!defined $tg_nsubst{$tmp[$i]});
				}
				@tmp=keys %{$hr_smpl_tg_subst_map->{$name}};
				for(my $i=0;$i<@tmp;$i++){
					push @sites,$tmp[$i] if(defined $bg_nsubst{$tmp[$i]});
				}
			}
			die "\nError: no mutations in foreground on the branch: $name!" unless scalar @sites;
			@sites=sort {$a<=>$b} @sites;
			for(my $i=0;$i<@sites;$i++){
				print $sites[$i];
				print ";" unless $i==@sites-1;
			}
		}
		if($f_phenotypes){
			print "\t";
			print $phenotype_map{$name} if defined $phenotype_map{$name};
		}
	}
}
			