#!/usr/bin/env perl
#This script calculates epistatic statistics and performs permutation test [Kryazhimsky 11]
#options:
#	-x <FN> - input file for EpiStat application in XPAR format [Kryazhimsky11]. Contains a tree in the adjacent file format
#	-m <"1|2|3">  Permutation Mode -  which set of mutations need to be randomized: 1 - background; 2 - foreground; 3 - both
#	[-p <uint>] Number of permutations
#		Default=0
#	[-R <uint>] Retry. Number of attempts to finish parallel calculations
#		Default=0 Do not retry
#	[-r] Rerun. Proceed previous unsucsessfully finished job
#	[-P] <FN> - input file with options strings for GNU Parallel
#		file format: epistat_script\tgnu_parallel_options
#			epistat_script=("shuffle_incidence_matrix.R"|"epistat.pl"|"sample2xparr.pl"|"allele_rate.sample2xparr.pl"|"fake_sample.pl")
#Params:
#<parameters> - epistat.pl configuration file
#<parameters> - configuration file
#		Contents:
#		[FDRSampleNum="<uint>"] - Number of samples used to calculate FDR
#			Default="0" - Do not calculate FDR
#		[TMPClean="0|1"] - Remove after finishing all temporary files and folders
#			Default="0" - Do not remove temp files
#		OutUpperPValues="Suffix" - Suffix used for pvalues' output file
#		OutLowerPValues="Suffix" - Suffix used for pvalues' output file
#		OutAvgEpiStats="Suffix" - Suffix used for average epistatic statistics output file
#		OutVarEpiStat="Suffix" - Suffix used for variances of epistatic statistics output file
#		OutSiteStats="Suffix" - Suffix used for an output file with average and variance of marginal epistatic statistics'  values for background and foreground sites
#		[OutCovOrdPairsEpiStat="Suffix"] - Suffix used for output file with covariances of epistatic statistics of two ordered site pairs. Covariances are calculated for internal epistasis only.
#			Default=".ord_site_pars.cov"
#		[OutFDR="Suffix"] - Suffix used for files with statistics for FDR samples
#		[AggrPhenotypesBy]="(best(_association)?|mixt(_associations)?)" - Method to aggregate samples from null distributions conditioned by phenotypes. 
#			!!!Required, if the input XPAR file declares two or more phenotype traits.
#			"best_association" - For each site pair the phenotype-conditioned null-model having the best unconditional associations with both sites of the site pair is preferred
#			"mixt_associations" - For each site pair the phenotype-conditioned null-model is selected proportional to unconditional associations of both sites of the site pair
#		[SplitBgrAndFgr="0|1"] - This flag allows independent permutations of mutations on the tree branches in background and foreground sites searching for intragenic associations.
#			This flag is always "1" for intergenic cases. The default value for intragenic cases is "0".
##Phenotypes:
#		[BgrSite2PhenoStat]="FN" - A file with association statistics of sites in background with phenotypes, e.g. two side pvalues (min(lower.pvalue,upper.pvalue))
#			!!!Required, if the input XPAR file declares two or more phenotype traits.
#		[FgrSite2PhenoStat]="FN" - A file with association statistics of sites in foreground with phenotypes, e.g. two side pvalues (min(lower.pvalue,upper.pvalue))
#			!!!Required, if the input XPAR file declares two or more phenotype traits.
#		[Site2PhenSortingOrder]="<(asc)ending|(desc)ending>" - A sorditing order of statistics from the best to the worst]
#			!!!Required, if the input XPAR file declares two or more phenotype traits.
#############
##Allele statistics:
#		[BgrTreatProfile="FN"] - The profile that is used to treat alleles for abbreviations of mutations in the background of the input XPARR file. 
#			The IQTree *.sitefreq format is expected!
#			!If 'TreatProfile' is not specified some scripts may fault trying to generate allele distributions
#		[BgrRootSequence="PROTEIN_SEQ"] - the string with protein sequence in the root of the tree specified in XPARR file.
#			This sequence is used to initialize alleles. If not specified the root sequence is generated from the profile.
#		[FgrTreatProfile="FN"] - The profile that is used to treat alleles for abbreviations of mutations in the foreground of the input XPARR file. 
#			The IQTree *.sitefreq format is expected!
#			!If 'TreatProfile' is not specified some scripts may fault trying to generate allele distributions
#		[FgrRootSequence="PROTEIN_SEQ"] - the string with protein sequence in the root of the tree specified in XPARR file.
#			This sequence is used to initialize alleles. If not specified the root sequence is generated from the profile.
#		[MkRootSequenceFrom=("xparr"|"align(ment)?")] - specifies how to generate root sequences if either of 'BgrRootSequence' or 'FgrRootSequence' parameters has not been specified.
#			"xparr" - root sequences are generated from the mutation distribution on the tree branches.
#			"align" - root sequences are sampled from the alignment profile.
#		[NullMutationModel=("neutral"|"allele_rates")] - specifies the type of mutational model used for generating the random distribution of mutations on the tree branches (null model).
#			"neutral" - this model saves the numbers of mutations in sites and on the tree branches. The mutation rates in sites are independent on alleles.
#			"allele_rates" - this model uses the estimated site and allele specific mutation rates for generating a random mutation distribution.
#			!This parameter is required
#############

use strict;
die "Set up the system variable \$EPISTAT_HOME" unless defined $ENV{EPISTAT_HOME};
use lib "$ENV{EPISTAT_LIB}";
use Getopt::Std;
use File::Basename;
use File::Path qw(make_path rmtree);
use EpistatNullModel::General;
use EpistatNullModel::Symmetric;
use IndexedFakeSamplesHeap;
use Phenotype::AggregateStatistics;

#server inst
my %commands;
$commands{"shuffle_incidence_matrix.R"}=[("$ENV{EPISTAT_HOME}/shuffle_incidence_matrix.R","")];
$commands{"epistat.pl"}=[("$ENV{EPISTAT_HOME}/epistat.pl","")];
$commands{"sample2xparr.pl"}=[("$ENV{EPISTAT_HOME}/sample2xparr.pl","")];
$commands{"allele_rate.sample2xparr.pl"}=[("$ENV{EPISTAT_HOME}/allele_rate.sample2xparr.pl","")];
$commands{"fake_sample.pl"}=[("$ENV{EPISTAT_HOME}/fake_sample.pl","")];
my $epistat_cmd=$commands{"epistat.pl"}->[0];
my $sample2xparr_cmd;
my $fake_sample_cmd=$commands{"fake_sample.pl"}->[0];

my %args;
if(!getopts('x:m:p:R:rP:',\%args)){
	die "\nError in option string!";
}
my $max_items_indir=1000;
my $files_heap;
my $xpar_fn=$args{x};
die "\nThe input XPAR file is not specified! Use -x option!" unless defined $xpar_fn;
my $f_mtx_mode;
if($args{m}=~/^[123]\s*$/){
	$f_mtx_mode=$args{m};
};
die "\nThe permutation mode is unknown or is not specified at all! Use -m option!" unless defined $f_mtx_mode;
my $nperm=0;
if(defined $args{p}){
	$nperm=$args{p};
	die "\nWrong value $nperm is specified for the option -p! Expects UINT!" unless $nperm=~/^\d+\s*$/;
};
my $stats_fn;
my $pairs_fn=".site_pairs";
my $idx2site_fn=".idx2site";
my $idx2branch_fn=".idx2branch";
my $mmtx_fn=".mmtx";
my $perm_fn=".mlist";
my $tau;
my $f_stat_type;
my $f_intragene;
my $f_splited_bgtg;
my $f_ignore_terminals=1;
my $f_rerun;
my $ntries=1;
my $f_mutation_model;
#for phenotypes:
my @pheno_labels;
my $npheno;
my $f_phenotype_pairs=0;
#specifyes the method of phenotypes aggregation: 0 - best_association, 1 - mean_association
my $f_pheno_aggr_method;
my $bg_sites2pheno_fn;
my $fg_sites2pheno_fn;
#0- ascending; 1- descending
my $site2phen_order;
my ($bgr_treat_prof_fn,$bgr_root_seq,$fgr_treat_prof_fn,$fgr_root_seq,$f_mk_root_seq);
my $alphabet_str;
$f_rerun=1 if defined $args{r};
if($args{R}=~/^\d+/){
	$ntries=$args{R}+1;
}else{
	die "\nWrong value for an option -R, UINT is expected!" if defined $args{R};
}
my $distance_mode;

if(defined $args{P}){
	open INPF, $args{P} or die "\nUnable to open input file: $args{P}!";
	while(<INPF>){
		chomp;
		my @line=split '\t';
		my $cmd=$line[0];
		my $run_parallel_opts=$line[1];
		die "\nError in the GNU Parallele option file format! Two columns are expected!" unless @line==2;
		die "\nThe valid program name expected in the 1st column of line:\n\t$_\nin -P ".$args{P} unless defined $commands{$cmd};
		if($run_parallel_opts=~/^\s*-{1,2}\S/){
			$run_parallel_opts=~s/^\s+//;
			$run_parallel_opts=~s/\s+$//;
			$run_parallel_opts.=" ";
			$commands{$cmd}->[1]=$run_parallel_opts;
		}else{
			die "\nUnable to parse an option string: $_!";
		}
	}
	close INPF;
}

open INFILE, "$ARGV[0]" or die "\nUnable to open input file: $ARGV[0]!";
while(<INFILE>){
	$_=$` if(/#/);
	if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
		my ($key,$value)=($1,$2);
		if($key eq "OutSitePairs"){
			$pairs_fn=$value;
		}elsif($key eq "StatType"){
			if($value=~m/^all_pairs/i){
				$f_stat_type=0;
			}elsif($value=~m/^independent_pairs/i){
				$f_stat_type=1;
			}elsif($value=~m/^branch_pairs/i){
				$f_stat_type=2;
			}elsif($value=~m/^all_allele_pairs/i){
				$f_stat_type=3;
			}elsif($value=~m/^indep_allele_pairs/i){
				$f_stat_type=3;
			}elsif($value=~m/^allele_branch_pairs/i){
				$f_stat_type=3;
			}else{
				die "\nUndefined type of statistics: StatType=$value!";
			}
		}elsif($key eq "OutEpiStatistics"){
			$stats_fn=$value;
		}elsif($key eq "IgnoreTerminals"){
			$f_ignore_terminals=$value;
			die "\nUnknown value ($value) for the IgnoreTerminals parameter in file $ARGV[0]!\n\tUINT in range is [0..1] expected!" unless $value=~/^[01]\s*$/;		
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
			$distance_mode=$value;
		}elsif($key eq "Phenotypes"){
			if($value>0&&$value<=2){
				$npheno=0;
				$f_phenotype_pairs=1 if $value==2;
			}else{
				die "\nUnpropper value for the 'Phenotype' parameter! Valid values: 0 or 1" unless $value==0;
			}
		}elsif($key eq "Alphabet"){
			$alphabet_str=$value;
			$alphabet_str=~s/\s//;
		}
	}
}
close INFILE;
die "\nError: Parameter PairsType is undefined in the file $ARGV[0]!" unless defined $f_intragene;
unless(defined $stats_fn){
	$stats_fn=".stat.exp";
	$stats_fn.=".$tau" if defined $tau;
}
if(defined $npheno){
	open INFILE, "<$xpar_fn" or die "\nUnable to open input file:$xpar_fn!";
	$_=<INFILE>; #skip header line
	chomp;
	s/\s+$//;
	my @line=split "\t";
	if($line[6]=~/^phenotypes:/){
		@pheno_labels=split ",", $';
		$npheno=@pheno_labels;
	}else{
		"\nNo phenotypes defined on branches of a tree in the file: $xpar_fn!";
	}
	close INFILE;
}

my $stats_ext=$stats_fn;
my $pairs_ext=$pairs_fn;

my $nfdr_samples=0;
my $f_tmp_clean=0;
my ($upper_pval_fn,$lower_pval_fn,$avg_fn,$var_fn,$fdr_fn);
my $site_stat_fn;
my $cov_fn=".ord_site_pars.cov";
open INFILE, "$ARGV[1]" or die "\nUnable to open input file: $ARGV[1]!";
while(<INFILE>){
	$_=$` if(/#/);
	if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
		my ($key,$value)=($1,$2);
		if($key eq "FDRSampleNum"){
			$nfdr_samples=$value;
			die "\nWrong value for FDRSampleNum accounted! Expects UINT!" unless $nfdr_samples=~/^\d+\s*$/;
		}elsif($key eq "TMPClean"){
			$f_tmp_clean=$value;
			die "\nWrong value for TMPClean accounted! Expects 0|1!" unless $f_tmp_clean=~/^[01]\s*$/;
		}elsif($key eq "OutUpperPValues"){
			$upper_pval_fn=$value;
		}elsif($key eq "OutLowerPValues"){
			$lower_pval_fn=$value;
		}elsif($key eq "OutAvgEpiStats"){
			$avg_fn=$value;
		}elsif($key eq "OutVarEpiStat"){
			$var_fn=$value;
		}elsif($key eq "OutCovOrdPairsEpiStat"){
			$cov_fn=$value;
		}elsif($key eq "OutFDR"){
			$fdr_fn=$value;
		}elsif($key eq "OutSiteStats"){
			$site_stat_fn=$value;
		}elsif($key eq "SplitBgrAndFgr"){
			$value=~s/\s+$//;
			$value=~s/^\s+//;
			$f_splited_bgtg=$value if $value=~/^0|1$/;
		}elsif($key eq "AggrPhenotypesBy"){
			if($value=~/best/){
				$f_pheno_aggr_method=0;
			}elsif($value=~/mixt/){
				$f_pheno_aggr_method=1;
			}else{
				die "\nUnknown phenotype aggregation method: $value!";
			}
		}elsif($key eq "BgrSite2PhenoStat"){
			$bg_sites2pheno_fn=$value;
		}elsif($key eq "FgrSite2PhenoStat"){
			$fg_sites2pheno_fn=$value;
		}elsif($key eq "Site2PhenSortingOrder"){
			if($value=~/^desc(ending)*/i){
				$site2phen_order=1;
			}elsif($value=~/^asc(ending)*/i){
				$site2phen_order=0;
			}
		}elsif($key eq "BgrTreatProfile"){
			$bgr_treat_prof_fn=$value;
		}elsif($key eq "BgrRootSequence"){
			$bgr_root_seq=uc $value;
			$bgr_root_seq=~s/\s//;
		}elsif($key eq "FgrTreatProfile"){
			$fgr_treat_prof_fn=$value;
		}elsif($key eq "FgrRootSequence"){
			$fgr_root_seq=uc $value;
			$fgr_root_seq=~s/\s//;
		}elsif($key eq "MkRootSequenceFrom"){
			if($value=~/^xparr/i){
				$f_mk_root_seq="xparr";
			}elsif($value=~/^align/i){
				$f_mk_root_seq="align";
			}else{
				die "\nUnknown method for generation of root sequences: $value!";
			}
		}elsif($key eq "NullMutationModel"){
			if($value=~/^neutral/i){
				$f_mutation_model=0;
			}elsif($value=~/^allele_rates/i){
				$f_mutation_model=1;
			}else{
				die "\nUnknown type of the null mutation model: $value!";
			}
		}else{
			die "\nUnknown parameter: $key in the input file $ARGV[1]!";
		}
	}
}
close INFILE;
if(defined $f_splited_bgtg){
	unless($f_intragene){
		die "\nThe parameter 'SplitBgrAndFgr' could not be set up to the 0 value!" if $f_splited_bgtg==0;
	}
}else{
	if($f_intragene){
		$f_splited_bgtg=0;
	}else{
		$f_splited_bgtg=1;
	}
}
if(defined($alphabet_str)){
	die "\nThe type of mutational model is not specified! Use the 'NullMutationModel' parameter!" unless defined $f_mutation_model;
	if($f_mutation_model==0){
		unless(defined($bgr_root_seq)&&defined($fgr_root_seq)){
			die "\nThe method to generate root sequence is not specified, use the 'MkRootSequenceFrom' parameter!" unless defined $f_mk_root_seq;
		}elsif(defined $f_mk_root_seq){
			warn "\nThe generating method for a root seqequence is specified but will be ignored!";
			$f_mk_root_seq=undef;
		}
	}
}else{
	unless(defined $f_mutation_model){
		$f_mutation_model=0;
	}elsif($f_mutation_model==1){
		die "\nThe specified type of mutational model is applicable only for the input data with alleles! May be 'Alphabet' is missed in the first parameter file!";
	}
	if(defined $f_mk_root_seq){
		warn "\nThe generating method for a root seqequence is specified but will be ignored!";
		$f_mk_root_seq=undef;
	}
}
die "\nError: The OutUpperPValues parameter has no defaul value!" unless defined $upper_pval_fn;
die "\nError: The OutLowerPValues parameter has no defaul value!" unless defined $lower_pval_fn;
die "\nError: The OutAvgEpiStats parameter has no defaul value!" unless defined $avg_fn;
die "\nError: The OutVarEpiStat parameter has no defaul value!" unless defined $var_fn;
die "\nError: The OutFDR parameter has no defaul value!" if $nfdr_samples>0&&!defined($fdr_fn);
die "\nError: The OutSiteStats parameter has no defaul value!" unless defined $site_stat_fn;
if(defined $npheno){
	if($npheno>1){
		die "\nThe NullMutationModel=\"allele_rates\" is not supported yet for input data with phenotypes!" if $f_mutation_model==1;
		die "\nUndefined aggregation method for phenotype-constrained fake samples!" unless defined $f_pheno_aggr_method;
		die "\nNo files with site to phenotype associations has been specified!" unless defined($bg_sites2pheno_fn)||defined($fg_sites2pheno_fn);
		die "\nThe best to worst sorting order for site to phenotype associations is required!" unless defined $site2phen_order;
		unless($f_intragene){
			die "\nThe file with association statistics of sites in background to phenotype hasn't been specified!" unless defined($bg_sites2pheno_fn);
			die "\nThe file with association statistics of sites in foreground to phenotype hasn't been specified!" unless defined($fg_sites2pheno_fn);
		}
	}else{
		die "\nError: at lest two phenotypes required to use the paired-phenotype model!" if $f_phenotype_pairs;
	}
}
if($f_mutation_model==0){
	$sample2xparr_cmd=$commands{"sample2xparr.pl"}->[0];
}else{
	$sample2xparr_cmd=$commands{"allele_rate.sample2xparr.pl"}->[0];
}

my ($basename,$dir,$ext) = fileparse($xpar_fn,'\.[^\.]*$');
$pairs_fn=$dir.$basename.$pairs_ext;
$stats_fn=$dir.$basename.$stats_ext;
if(!$f_rerun){
	#Deleting previous versions of files
	unlink $pairs_fn if(-e $pairs_fn);
	unlink $stats_fn if(-e $stats_fn);
}

#indices in @bg_files and @tg_files: 0-matrix; 1-idx2branch; 2-idx2site
my @bg_files;
my @tg_files;
my @mmtx_prefs;
if($f_mutation_model==0){
	if(defined $npheno){
		if($f_phenotype_pairs){
			for(my $i=0;$i<$npheno;$i++){
				for(my $j=$i+1;$j<$npheno;$j++){
					push @mmtx_prefs,[($pheno_labels[$i]."_".$pheno_labels[$j],"#0")];
					push @mmtx_prefs,[($pheno_labels[$i]."_".$pheno_labels[$j],"#1")];
					push @mmtx_prefs,[($pheno_labels[$i]."_".$pheno_labels[$j],"#2")];
					push @mmtx_prefs,[($pheno_labels[$i]."_".$pheno_labels[$j],"#3")];
				}
			}
		}else{
			for(my $i=0;$i<$npheno;$i++){
				push @mmtx_prefs,[($pheno_labels[$i],"#0")];
				push @mmtx_prefs,[($pheno_labels[$i],"#1")];
			}
		}
	}else{
		push @mmtx_prefs,[];
	}
	foreach my $rpref(@mmtx_prefs){
		my $fn=$dir.$basename;
		$fn.=".".join "",@{$rpref} if defined $npheno;
		if($f_intragene){
			$fn.=".intra";
		}else{
			$fn.=".inter";
		};
		my $str="";
		#$f_mtx_mode is a permutation mode: 1- bgr; 2- fgr; 3- both
		if($f_mtx_mode!=2){
			if($f_intragene&&$f_mtx_mode==3){
				$str=".both";
			}else{
				$str=".bgr";
			}
			push @bg_files,[($fn.$str.$mmtx_fn,$fn.$str.$idx2branch_fn,$fn.$str.$idx2site_fn)];
		}
		if($f_mtx_mode!=1){
			if($f_intragene&&$f_mtx_mode==3){
				$str=".both";
			}else{
				$str=".fgr";
			}
			if($f_splited_bgtg){
				#$f_splited_bgtg==1 for the intergene mode
				#$f_splited_bgtg==0 for the infragene mode by default
				push @tg_files,[($fn.$str.$mmtx_fn,$fn.$str.$idx2branch_fn,$fn.$str.$idx2site_fn)];
			}
		}
	}
	if(!$f_rerun){
		foreach my $fset(@bg_files){
			foreach my $fn(@{$fset}){
				unlink $fn if -e $fn;
			}
		}
		foreach my $fset(@tg_files){
			foreach my $fn(@{$fset}){
				unlink $fn if -e $fn;
			}
		}
	}
}else{
	push @mmtx_prefs,[];
}

sub print_child_termination_status{
	if ($? == -1) {
		print "failed to execute: $!\n";
		exit 1;
	}elsif ($? & 127) {
		printf "\tchild died with signal %d, %s coredump\n",
					($? & 127),  ($? & 128) ? 'with' : 'without';
	}else {
		if($? >> 8!=0){
			printf "\tchild exited with value %d\n", $? >> 8;
		};
	}
}

my $str=$epistat_cmd." -x $xpar_fn ";
$str.="-m $f_mtx_mode " if $f_mutation_model==0;
if(defined $fgr_treat_prof_fn){
	$str.="-A -F $fgr_treat_prof_fn ";
}
if(defined $bgr_treat_prof_fn){
	$str.="-A -B $bgr_treat_prof_fn ";
}
$str.=$ARGV[0];
unless($f_rerun&&(-e $pairs_fn)&&(-e $stats_fn)){
	print STDERR "Calculating epistatic statistics for the input tree: $xpar_fn\n$str";
	system($str);
	print_child_termination_status();
}
die "\n\tFatal Error: The expected file does not exist: $pairs_fn" unless -e $pairs_fn;
die "\n\tFatal Error: The expected file does not exist: $stats_fn" unless -e $stats_fn;
my @out_file_types;
my @out_fake_file_types;
{
	my $pvalue_mode_str="";
	if($f_intragene){
		$pvalue_mode_str=".intragene.unord_pairs";
		push @out_file_types,$pvalue_mode_str.$lower_pval_fn;
		push @out_file_types,$pvalue_mode_str.$upper_pval_fn;
		@out_fake_file_types=@out_file_types;
		if($f_stat_type==3){
			$pvalue_mode_str="";
			push @out_file_types,$site_stat_fn;
		}else{
			$pvalue_mode_str=".intragene.ord_pairs";
			push @out_file_types,$cov_fn;
		}
	}
	push @out_file_types,$avg_fn;
	push @out_file_types,$var_fn;
	unless(($f_stat_type==3)&&$f_intragene){
		push @out_file_types,$pvalue_mode_str.$lower_pval_fn;
		push @out_file_types,$pvalue_mode_str.$upper_pval_fn;
		push @out_fake_file_types,$pvalue_mode_str.$lower_pval_fn;
		push @out_fake_file_types,$pvalue_mode_str.$upper_pval_fn;
		push @out_file_types,".bgr".$site_stat_fn;
		push @out_file_types,".fgr".$site_stat_fn;
	}
}

sub check_outfiles{
	my $fn_prefix=shift;
	my @fn_suffixes=@_;
	foreach my $suff(@fn_suffixes){
		my $fn=$fn_prefix.$suff;
		return 0 unless (-e $fn);
	}
	return 1;
}

sub print_statistics{
	my ($out_file_name,$header)=@_;
	die "\nError print_statistics(): Nothing to print. No arrays with statistics!" unless scalar(@_)-2;
	my $nlines=@{$_[2]};
	for(my $j=3;$j<@_;$j++){
		my $ra_stat=$_[$j];
		my $n=@{$ra_stat};
		die "\nError print_statistics(): different number of lines $n vs. $nlines in an input array: argument no.= $j!" unless $n==$nlines;
	}
	open OPF, ">$out_file_name" or die "\nUnable to open output file:$out_file_name!";
	if(defined $header){
		print OPF "$header\n";
	}
	my $prec=int(log($nperm)/log(10))+2;
	for(my $i=0;$i<$nlines;$i++){
		my $str="";
		for(my $j=2;$j<@_;$j++){
			my $ra_stat=$_[$j];
			$str.="\t" if $j>2;
			if(!defined $ra_stat->[$i]){
				$str.="NA";
			}else{
				if(abs($ra_stat->[$i]-int($ra_stat->[$i]))>0){
					$str.= sprintf("%.".$prec."f", $ra_stat->[$i]);
				}else{
					$str.=$ra_stat->[$i];
				}
			}
		}
		print OPF "$str\n";
	}
	close OPF;
}

sub read_statistics{
	my ($stats_fn,$ra_stat)=@_;
	open INPF, "<$stats_fn" or die "\nUnable to open input file: $stats_fn!";
	@{$ra_stat}=();
	while(<INPF>){
		chomp;
		$_=~s/\s*$//;
		my @line=split "\t";
		if($line[0] ne ""){
			push @{$ra_stat},$line[0];
		}
	}
	close INPF;
}

sub gen_tempname{
	my $nchar=shift;
	my @chars = ( "A" .. "Z", "a" .. "z", 0 .. 9 );
	return join("", @chars[ map { rand @chars } ( 1 .. $nchar ) ]);
}

sub count_samples{
	my ($path,$ext,$ra_samples_storage)=@_;
	$path=~s/\s*$//;
	$path=~s/\/$//;
	my $n=0;
	foreach my $si(@{$ra_samples_storage}){
		my $folder=$path.$si->[2]."/";
		my $from=$si->[0];
		my $to=$si->[1];
		for(my $i=$from;$i<=$to;$i++){
			my $str=$folder.$i.$ext;
			if(-e $str){
				my $filesize=-s $str;
				$n++ if($filesize);
			}
		}
	}
	return $n;
}

sub find_missed_samples{
	my ($path,$ext,$ra_samples_storage_node)=@_;
	my @gaps;
	$path=~s/\s*$//;
	$path=~s/\/$//;
	$path.=$ra_samples_storage_node->[2]."/";
	my $from=$ra_samples_storage_node->[0];
	my $to=$ra_samples_storage_node->[1];
	for(my $i=$from;$i<=$to;$i++){
		my $str=$path.$i.$ext;
		if(-e $str){
			my $filesize=-s $str;
			push @gaps,$i unless($filesize);
		}else{
			push @gaps,$i;
		}
	}
	return @gaps;
}

sub read_sites2pheno{
	my ($infile,$ra_pheno_labels,$nsamples,$rh_sites2pheno)=@_;
	my $alpha=0.5;
	my $lp0=log($alpha)-log($nsamples);
	open INPF, "<$infile" or die "\nUnable to open input file:$infile!";
	#reading the header line
	$_=<INPF>;
	chomp;
	s/\s+$//;
	my @line=split "\t";
	my $npheno=@{$ra_pheno_labels};
	die "\nError in the header line of the $infile:\n$_!" unless ($line[0] eq "site")&&(@line==$npheno+1);
	for(my $i=0;$i<$npheno;$i++){
		die "\nInconsistence in the phenotypes labels:\ni=$i\t$line[$i+1] expected ".$ra_pheno_labels->[$i] unless $line[$i+1] eq $ra_pheno_labels->[$i];
	}
	%{$rh_sites2pheno}=();
	my $I=2;
	while(<INPF>){
		chomp;
		s/^\s+//;
		s/\s+$//;
		if(/\S+/){
			my @line=split "\t";
			die "\nWrong number of columns in the data row $I in file $infile!" unless @line==$npheno+1;
			$rh_sites2pheno->{$line[0]}=[];
			foreach my $val(@line[1..$npheno]){
				if($val>0){
					push @{$rh_sites2pheno->{$line[0]}},log($val);
				}elsif($val==0){
					push @{$rh_sites2pheno->{$line[0]}},$lp0;
				}else{
					die "\nError in read_sites2pheno(): in $infile line:\n\t $_\n\tthe negative value accounted!";
				}
			}
			for(my $j=1;$j<@line;$j++){
				die "\nUnable to use negative association statistics: site=$line[0],phenotype=$j" if $line[$j]<0;
			}
			$I++;
		}
	}
	close INPF;
}

sub make_sites2pheno_pairs{
	my ($npheno,$rh_sites2pheno,$rh_sites2pheno_pairs)=@_;
	%{$rh_sites2pheno_pairs}=();
	foreach my $site (keys %{$rh_sites2pheno}){
		$rh_sites2pheno_pairs->{$site}=[];
		for(my $i=0;$i<$npheno;$i++){
			my $lph1=$rh_sites2pheno->{$site}->[$i];
			for(my $j=$i+1;$j<$npheno;$j++){
				my $lph2=$rh_sites2pheno->{$site}->[$j];
				push @{$rh_sites2pheno_pairs->{$site}},$lph1+$lph2;
			}
		}
	}
}

sub dump_into{
	my ($ra_data,$tmp_fname)=@_;
	open OPF, ">$tmp_fname" or die "\nUnable to open output file:$tmp_fname!";
	my $i=shift @{$ra_data};
	print OPF "$i";
	foreach $i(@{$ra_data}){
		print OPF "\n$i";
	}
	close OPF;
}

sub run_shuffle_incidence_matrix{
	#Contents of @{$ini_files} has been already generated by the call of epistat.pl!
	my $ini_files=$_[0];
	my $outdir=$_[1];
	my @samples_storage=@{$_[2]};
	my $perm_fn=$_[3];
	my $brw_R=$commands{"shuffle_incidence_matrix.R"}->[0];
	my $parallel_opts=$commands{"shuffle_incidence_matrix.R"}->[1];
	if(@{$ini_files}){
		my $n=0;
		if(!$f_rerun&&(-d $outdir)){
			rmtree($outdir);
		}
		my $t=1;
		foreach my $fn(@{$ini_files}){
			if(! -e $fn){
				warn "\n\tFatal Error: The required file is not exists: $fn";
				$t=0;
			}
		}
		die unless $t;
		for(my $i=0;$i<@samples_storage;$i++){
			my $from=$samples_storage[$i]->[0];
			my $to=$samples_storage[$i]->[1];
			my $mmtx_fn=$ini_files->[0];
			my $out_path=$outdir.$samples_storage[$i]->[2];
			my $str;
			my $tmp_fname;
			if($f_rerun){
				my @gap=find_missed_samples($outdir,$perm_fn,$samples_storage[$i]);
				if(@gap){
					$tmp_fname=gen_tempname(10).$perm_fn.".missed";
					dump_into(\@gap,$tmp_fname);
					$str="cat $tmp_fname";
				}
			}else{$str="seq $from $to";};
			if(defined $str){
				$str.="| parallel ".$parallel_opts."$brw_R $mmtx_fn $out_path {}";
				print STDERR "\n\nGenerating permutations\n$str";
				make_path($out_path) unless -d $out_path;
				system($str);
				print_child_termination_status();
				unlink $tmp_fname if defined $tmp_fname;
			}
		}
	}
}

if($nperm){
	my @samples_storage;
	if($nperm>$max_items_indir){
		$files_heap=IndexedFakeSamplesHeap->new($nperm,$max_items_indir);
		my @tmp=$files_heap->get_storage_path_list;
		for(my $i=0;$i<@tmp-1;$i++){
			push @samples_storage,[($max_items_indir*$i+1,$max_items_indir*($i+1),"/".$tmp[$i])];
		}
		push @samples_storage,[($max_items_indir*(@tmp-1)+1,$nperm,"/".$tmp[-1])];
	}else{
		push @samples_storage,[(1,$nperm,"")];
	}
	my $outdir=$dir.$basename;
	my $samples_dir=$outdir."/samples";
	my $pvalue_mode_str="";
	unless($f_rerun && check_outfiles($outdir,@out_file_types)){
		if($f_mutation_model==0){
			for(my $ids=0;$ids<@mmtx_prefs;$ids++){
				my $outdir=$dir.$basename;
				$outdir.="/pheno/".join "/",@{$mmtx_prefs[$ids]} if(defined $npheno);
				#$perm_fn - an extension for a fake mutation distribution file ('.mlist' by default)
				#scalar(@bg_files) equals to the total number of background subsets of branches
				if(@bg_files){
					run_shuffle_incidence_matrix($bg_files[$ids],$outdir."/bgr",\@samples_storage,$perm_fn);
					my $n=count_samples($outdir."/bgr",$perm_fn,\@samples_storage);
					die "\nError: There are less than expected *$perm_fn files in: $outdir/bgr!" unless($n>=$nperm);
				}
				#scalar(@bg_files) equals to the total number of foreground subsets of branches
				if(@tg_files){
					run_shuffle_incidence_matrix($tg_files[$ids],$outdir."/fgr",\@samples_storage,$perm_fn);
					my $n=count_samples($outdir."/fgr",$perm_fn,\@samples_storage);
					die "\nError: There are less than expected *$perm_fn files in: $outdir/fgr!" unless($n>=$nperm);
				}
			}
		}
		#die("\nWrong number of file sets!") unless scalar(@bg_files)==scalar(@tg_files)||scalar(@bg_files)==0||scalar(@tg_files)==0;
		my $step;
		if(defined $npheno){
			if($f_phenotype_pairs){
				$step=4;
			}else{
				$step=2;
			}
		}else{
			$step=1;
		}
		#iterate phenotypes
		for(my $fsi=0;$fsi<@mmtx_prefs;$fsi+=$step){
			my $sample2xparr_prm=$outdir;
			#$mmtx_prefs[$fsi]->[0] - a phenotype label
			$sample2xparr_prm.=".".$mmtx_prefs[$fsi]->[0] if defined $npheno;
			$sample2xparr_prm.=".sample2xparr.prm";
			open OPF, ">$sample2xparr_prm" or die "\nUnable to open output file:$sample2xparr_prm!";
			print OPF "XPAR=\"$xpar_fn\"";
			print OPF "\nPairsType=";
			if($f_intragene){
				if($f_splited_bgtg){
					print OPF "\"intergene\"";
				}else{
					print OPF "\"intragene\"";
				}
			}else{
				print OPF "\"intergene\"";
			}
			print OPF "\nMatrixMode=\"$f_mtx_mode\""; #$f_mtx_mode: 1- bgr; 2- fgr; 3- both
			if(defined $alphabet_str){
				print OPF "\nAlphabet=\"$alphabet_str\"";
				print OPF "\nBgrTreatProfile=\"$bgr_treat_prof_fn\"" if defined $bgr_treat_prof_fn;
				print OPF "\nBgrRootSequence=\"$bgr_root_seq\"" if(defined $bgr_root_seq);
				print OPF "\nFgrTreatProfile=\"$fgr_treat_prof_fn\"" if defined $fgr_treat_prof_fn;
				print OPF "\nFgrRootSequence=\"$fgr_root_seq\"" if(defined $fgr_root_seq);
			}
			if($f_mutation_model==0){
				if(defined $alphabet_str){
					print OPF "\nMkRootSequenceFrom=\"$f_mk_root_seq\"" if defined $f_mk_root_seq;
				}
				print OPF "\nIgnoreTerminals=\"$f_ignore_terminals\"";
				print OPF "\nNBranchSubsets=\"$step\"";
				#print parameters for each branch subset within a phenotype
				for(my $i=0;$i<$step;$i++){
					print OPF "\nSubsetFolder=\"".$outdir;
					print OPF "/pheno/".join("/",@{$mmtx_prefs[$fsi+$i]}) if defined $npheno;
					print OPF "\"";
					if(@bg_files>0){
						print OPF "\nBgrIdx2Site=\"".$bg_files[$fsi+$i]->[2]."\"";
						print OPF "\nBgrIdx2Branch=\"".$bg_files[$fsi+$i]->[1]."\"";
					}
					if(@tg_files>0){
						print OPF "\nFgrIdx2Site=\"".$tg_files[$fsi+$i]->[2]."\"";
						print OPF "\nFgrIdx2Branch=\"".$tg_files[$fsi+$i]->[1]."\"";
					}
				}
			}else{
				print OPF "\nBranchDistUnits=\"$distance_mode\"";
			}
			close OPF;
			my $samples_dir=$outdir;
			$samples_dir.="/pheno/".$mmtx_prefs[$fsi]->[0] if defined $npheno;
			$samples_dir.="/samples";
			for(my $i=0;$i<@samples_storage;$i++){
				my $from=$samples_storage[$i]->[0];
				my $to=$samples_storage[$i]->[1];
				my $out_path=$samples_storage[$i]->[2];
				make_path($samples_dir.$out_path);
				my @gap;
				my $tmp_fname;
				my $cmd_str;
				$str="";
				if($f_mutation_model==0){
					if(@bg_files>0){
						$str.="-b bgr".$out_path."/{}$perm_fn";
					}
					if(@tg_files>0){
						$str.=" " if(@bg_files>0);
						$str.="-t fgr".$out_path."/{}$perm_fn";
					}
				}
				if($f_rerun){
					@gap=find_missed_samples($samples_dir,".xpar",$samples_storage[$i]);
					if(@gap){
						$tmp_fname=gen_tempname(10).".xpar.missed";
						dump_into(\@gap,$tmp_fname);
						$cmd_str="cat $tmp_fname";
					}
				}else{$cmd_str="seq $from $to";};
				if($cmd_str){
					my $run_parallel_opts;
					if($f_mutation_model==0){
						$run_parallel_opts=$commands{"sample2xparr.pl"}->[1];
					}elsif($f_mutation_model==1){
						$run_parallel_opts=$commands{"allele_rate.sample2xparr.pl"}->[1];
					}
					$str=$cmd_str." | parallel ".$run_parallel_opts."$sample2xparr_cmd $str $sample2xparr_prm"."\">\"".$samples_dir.$out_path."/{}.xpar";
					print STDERR "\n\nGenerating XPAR files for randomized trees!\n$str";
					system($str);
					print_child_termination_status();
					unlink $tmp_fname if defined $tmp_fname;
				}
			}
			my $n=count_samples($samples_dir,".xpar",\@samples_storage);
			die "\nError: There are less than expected *.xpar files in: $samples_dir!" unless $n>=$nperm;
			for(my $i=0;$i<@samples_storage;$i++){
				my $from=$samples_storage[$i]->[0];
				my $to=$samples_storage[$i]->[1];
				my @gap;
				my $tmp_fname;
				my $cmd_str;
				$str="-x ".$samples_dir.$samples_storage[$i]->[2]."/{}.xpar -p";
				if($f_rerun){
					@gap=find_missed_samples($samples_dir,$stats_ext,$samples_storage[$i]);
					if(@gap){
						$tmp_fname=gen_tempname(10).".stat.missed";
						dump_into(\@gap,$tmp_fname);
						$cmd_str="cat $tmp_fname";
					}
				}else{$cmd_str="seq $from $to";}
				if($cmd_str){
					my $run_parallel_opts=$commands{"epistat.pl"}->[1];
					$str=$cmd_str." | parallel ".$run_parallel_opts."$epistat_cmd $str ".$ARGV[0];
					print STDERR "\n\nCalculating epistatic statistics for randomized trees!\n$str";
					system($str);
					print_child_termination_status();
					unlink $tmp_fname if defined $tmp_fname;
				}
			}
			$n=count_samples($samples_dir,$stats_ext,\@samples_storage);
			die "\nError: There are less than expected *$stats_ext files in: $samples_dir!" unless $n>=$nperm;
		}
		#Mixture phenotype distributions
		if(defined $npheno){
			if((!$f_rerun)&&(-d $samples_dir)){
				rmtree($samples_dir);
			}
			my @pheno_samples_dirs;
			if($f_phenotype_pairs){
				for(my $i=0;$i<$npheno;$i++){
					for(my $j=$i+1;$j<$npheno;$j++){
						push @pheno_samples_dirs, $pheno_labels[$i]."_".$pheno_labels[$j];
					}
				}
			}else{
				@pheno_samples_dirs=@pheno_labels;
			}
			for(my $i=0;$i<@pheno_samples_dirs;$i++){
				$pheno_samples_dirs[$i].="/samples";
			}
			if(@pheno_samples_dirs>1){
				my $rh_bg_sites2pheno;
				my $rh_fg_sites2pheno;
				if(defined $bg_sites2pheno_fn){
					$rh_bg_sites2pheno={};
					if($f_phenotype_pairs){
						my %tmp;
						read_sites2pheno($bg_sites2pheno_fn,\@pheno_labels,$nperm,\%tmp);
						make_sites2pheno_pairs($npheno,\%tmp,$rh_bg_sites2pheno);
					}else{
						read_sites2pheno($bg_sites2pheno_fn,\@pheno_labels,$nperm,$rh_bg_sites2pheno);
					}
				}
				if(defined $fg_sites2pheno_fn){
					unless($f_intragene&&($bg_sites2pheno_fn eq $fg_sites2pheno_fn)){
						$rh_fg_sites2pheno={};
						if($f_phenotype_pairs){
							my %tmp;
							read_sites2pheno($fg_sites2pheno_fn,\@pheno_labels,$nperm,\%tmp);
							make_sites2pheno_pairs($npheno,\%tmp,$rh_fg_sites2pheno);
						}else{
							read_sites2pheno($fg_sites2pheno_fn,\@pheno_labels,$nperm,$rh_fg_sites2pheno);
						}
					}
				}
				if($f_intragene){
					$rh_fg_sites2pheno=$rh_bg_sites2pheno unless defined $rh_fg_sites2pheno;
					$rh_bg_sites2pheno=$rh_fg_sites2pheno unless defined $rh_bg_sites2pheno;
				}
				#$site2phen_order - the best to worst sorting order of site to phenotype association statistics
				if($f_pheno_aggr_method==0){
					Phenotype::AggregateStatistics::aggregateByBest($pairs_fn,
						$f_intragene,$outdir."/pheno",\@pheno_samples_dirs,\@samples_storage,$stats_ext,
						$rh_bg_sites2pheno,$rh_fg_sites2pheno,$site2phen_order,$samples_dir,$ntries);
				}else{
					Phenotype::AggregateStatistics::aggregateByMixture($pairs_fn,
						$f_intragene,$outdir."/pheno",\@pheno_samples_dirs,\@samples_storage,$stats_ext,
						$rh_bg_sites2pheno,$rh_fg_sites2pheno,$site2phen_order,$samples_dir,$ntries);
				}
				%{$rh_bg_sites2pheno}=();
				%{$rh_fg_sites2pheno}=();
			}else{
				my $indir=$outdir."/pheno/".$pheno_samples_dirs[0];
				make_path($samples_dir) unless -d $samples_dir;
				my $str="mv $indir/* $samples_dir";
				system($str);
				print_child_termination_status();
			}
			my $n=count_samples($samples_dir,$stats_ext,\@samples_storage);
			die "\nError: There are less than expected *$stats_ext files in: $samples_dir!" unless $n>=$nperm;
		}
		my @stat1;
		my @stat2;
		my @bgr_sites;
		my @fgr_sites;
		print STDERR "\n\nBuilding zero model distribution and calculating statistics!";
		my $h0;
		#'$f_stat_type==3' - allele based epistatic statistics
		unless(($f_stat_type==3)&&$f_intragene){
			if($nperm>$max_items_indir){
				$h0=EpistatNullModel::General->new("-size" => $nperm,"-samples_dir" => $samples_dir,"-stats_ext" => $stats_ext,
					"-obs_stats_fn" => $stats_fn,"-obs_pairs_fn" => $pairs_fn, 
					"-is_intragene" => $f_intragene, "-max_samples_indir" => $max_items_indir);
			}else{
				$h0=EpistatNullModel::General->new("-size" => $nperm,"-samples_dir" => $samples_dir,"-stats_ext" => $stats_ext,
					"-obs_stats_fn" => $stats_fn,"-obs_pairs_fn" => $pairs_fn, "-is_intragene" => $f_intragene);
			}
			#Print marginal epistatic statistics for sites
			$h0->get_sites(\@bgr_sites,\@fgr_sites);
			$h0->calc_bgr_moments12(\@stat1,\@stat2);
			$str=$outdir.".bgr".$site_stat_fn;
			print_statistics($str,"site\tavg\tvar",\@bgr_sites,\@stat1,\@stat2);
			$h0->calc_fgr_moments12(\@stat1,\@stat2);
			$str=$outdir.".fgr".$site_stat_fn;
			print_statistics($str,"site\tavg\tvar",\@fgr_sites,\@stat1,\@stat2);
			#########
			if($f_intragene){
				$h0->calc_ordered_site_pairs_covariances(\@stat1);
				$str=$outdir.$cov_fn;
				print_statistics($str,undef,\@stat1);
				$pvalue_mode_str=".intragene.ord_pairs";
			}
			$h0->calc_pvalues(\@stat1,\@stat2);
			$str=$outdir.$pvalue_mode_str.$lower_pval_fn;
			print_statistics($str,undef,\@stat1);
			$str=$outdir.$pvalue_mode_str.$upper_pval_fn;
			print_statistics($str,undef,\@stat2);
		}else{
			if($nperm>$max_items_indir){
				$h0=EpistatNullModel::Symmetric->new("-size" => $nperm,"-samples_dir" => $samples_dir,"-stats_ext" => $stats_ext,
					"-obs_stats_fn" => $stats_fn,"-obs_pairs_fn" => $pairs_fn, 
					"-max_samples_indir" => $max_items_indir);
			}else{
				$h0=EpistatNullModel::Symmetric->new("-size" => $nperm,"-samples_dir" => $samples_dir,"-stats_ext" => $stats_ext,
					"-obs_stats_fn" => $stats_fn,"-obs_pairs_fn" => $pairs_fn);
			}
			#Print marginal epistatic statistics for sites
			$h0->get_sites(\@bgr_sites,\@fgr_sites);
			$h0->calc_site_moments12(\@stat1,\@stat2);
			$str=$outdir.$site_stat_fn;
			print_statistics($str,"site\tavg\tvar",\@bgr_sites,\@stat1,\@stat2);
			########
		}
		$h0->calc_moments12(\@stat1,\@stat2);
		$str=$outdir.$avg_fn;
		print_statistics($str,undef,\@stat1);
		$str=$outdir.$var_fn;
		print_statistics($str,undef,\@stat2);
		if($f_intragene){
			$h0->calc_unordered_pairs_pvalues(\@stat1,\@stat2);
			$pvalue_mode_str=".intragene.unord_pairs";
			$str=$outdir.$pvalue_mode_str.$lower_pval_fn;
			print_statistics($str,undef,\@stat1);
			$str=$outdir.$pvalue_mode_str.$upper_pval_fn;
			print_statistics($str,undef,\@stat2);
		}
		undef($h0);
		@stat1=();
		@stat2=();
	}
	if($nfdr_samples){
		my $fdrdir=$outdir."/"."fdr";
		$_=make_path($fdrdir);
		my %fdr_samples;
		my $n;
		if($f_rerun&&$_==0){
			opendir(DIR, $fdrdir) or die "can't opendir $fdrdir: $!";
			my $file;
			while (defined($file = readdir(DIR))) {
				my ($basename,$dir,$ext) = fileparse($file,'\..*$');
				if($ext eq $stats_ext){
					$fdr_samples{$basename}=1;
				}
			}
			closedir(DIR);
			$n=$nfdr_samples - scalar(keys %fdr_samples);
		}else{
			$n=$nfdr_samples;
		};
		print STDERR "\n\nStarting FDR estimation for $nfdr_samples samples:";
		my $k=0;
		while($k<$n){
			$str=int(rand($nperm)) +1;
			if(!defined $fdr_samples{$str}){
				$fdr_samples{$str}=1;
				$k++;
			}
		}
		my @fdr_samples_keys=keys %fdr_samples;
		%fdr_samples=();
		my $tmp_fname=gen_tempname(10);
		$tmp_fname=$fdrdir."/".$tmp_fname;
		open OPF, ">$tmp_fname" or die "\nUnable to open output file: $tmp_fname!";
		for(my $j=0;$j<$nfdr_samples;$j++){
			print OPF $fdr_samples_keys[$j];
			print OPF "\n" if $j<$nfdr_samples-1;
		}
		close OPF;
		$str="-d $outdir -p $nperm -s $pairs_fn";
		$str.=" -m $max_items_indir" if $nperm>$max_items_indir;
		$str.=" -i {}";
		{
			my $run_parallel_opts=$commands{"fake_sample.pl"}->[1];
			$str="cat $tmp_fname | parallel ".$run_parallel_opts."$fake_sample_cmd $str ".$ARGV[0]." ".$ARGV[1];
		}
		for(my $I=0;$I<$ntries;$I++){
			print STDERR "\n\nTry$I of FDR estimation for $nfdr_samples samples in parallel:!\n$str";
			system($str);
			print_child_termination_status();
			$n=0;
			for(my $j=0;$j<$nfdr_samples;$j++){
				$k=$fdr_samples_keys[$j];
				$n++ if -e $fdrdir."/".$k.$stats_ext && check_outfiles($fdrdir."/".$k."_fake",@out_fake_file_types);
			}
			last if $n==$nfdr_samples;
		}
		die "\nSome of FDR samples were not processed: $n\t$nfdr_samples!" unless $n==$nfdr_samples;
		unlink $tmp_fname;
	}
	if($f_tmp_clean){
		#unlink $sample2xparr_prm;
		if(defined $npheno){
			rmtree($outdir."/pheno");
		}else{
			rmtree($outdir."/"."bgr") if(@bg_files>0);
			rmtree($outdir."/"."fgr") if(@tg_files>0);
		}
		rmtree($samples_dir);
		rmtree($outdir) unless $nfdr_samples;
	}
}
