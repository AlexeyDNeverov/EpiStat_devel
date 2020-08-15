#!/usr/bin/env perl
#This script calculates Peacks Over the Threshold approximation for empirical pvalues
#Usage: options <parameters>
#options
#	-N <SamplesNum> - A number of samples from a null model
#	-d <SamplesRoot> - A root of a tree of storage directories with fake samples
#	-n <MaxSamplesNum> - Maximal number of samples in a storage directory
#	-p <PvalueTheshold> in [0,0.01] - Upper bound for pvalues to select pairs for POT approximation
#<parameters> - configuration file
#		Contents:
#		[PairsOrdering]="0"|"1" - The type of site pairs: 0 - unordered or 1 - ordered
#		PairsType=("intra[gene]"|"inter[gene]")
#		Pairs="FN" - A file with site pairs
#		Epistat_Ext="STR" - An extation of files with epistatic statistics
#		ObsEpistat="FN" - A file name with observed epistatic statistics
#		ObsPvalues="FN" - A file with pvalues
#Note: A parameter file for estimate_fdr.pl could be used here!

#server inst
use lib "$ENV{EPISTAT_HOME}";
#desctop inst
#use lib "/cygdrive/d/EvoLab/devel/epistat";

use strict;
use Getopt::Std;
use PValuesPOTapprox;

my $f_intragene;
my $f_ordered_pairs=1;
my ($pvalues_fn,$pairs_fn,$epistat_ext,$epistat_fn);
my ($nsamples,$samples_dir,$max_samples_indir,$max_pvalue);

my %args;
if(!getopts('d:n:p:N:',\%args)){
	die "\nError in option string!";
}
$nsamples=$args{N};
die "\nNumber of samples from a null model is required: use -N option!" unless defined $nsamples;
$samples_dir=$args{d};
die "\nA directory with fake samples was not specified: use -d option!" unless defined $samples_dir;
$max_samples_indir=$args{n};
die "\nA maximal number of samples in subdirectories of $samples_dir was not specified: use -n option!" unless defined $max_samples_indir;
$max_pvalue=$args{p};
die "\nA POT approximation could not be applied for all pvalues please specify a threshold <=1%: use -p option!" unless defined($max_pvalue)&&($max_pvalue>=0)&&($max_pvalue<=0.01);

open INFILE, "$ARGV[0]" or die "\nUnable to open input file: $ARGV[0]!";
while(<INFILE>){
	$_=$` if(/#/);
	if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
		my ($key,$value)=($1,$2);
		if($key eq "ObsPvalues"){
			$pvalues_fn=$value;
		}elsif($key eq "Pairs"){
			$pairs_fn=$value;
		}elsif($key eq "Epistat_Ext"){
			$epistat_ext=$value;
		}elsif($key eq "ObsEpistat"){
			$epistat_fn=$value;
		}elsif($key eq "PairsType"){
			if($value=~m/^intra(gene)*/i){
				$f_intragene=1;
			}elsif($value=~m/^inter(gene)*/i){
				$f_intragene=0;
			}else{
				die "\nUnknown value ($value) for the PairsType parameter in file $ARGV[0]!"; 
			}
		}elsif($key eq "PairsOrdering"){
			if($value==0||$value==1){
				$f_ordered_pairs=$value;
			}else{
				die "\nUnknown value ($value) for the 'PairsOrdering' parameter in file $ARGV[0]!"; 
			}
		}
	}
}
close INFILE;
my $smpl=PValuesPOTapprox->new("-size" => $nsamples, "-samples_dir" => $samples_dir, "-max_samples_indir" => $max_samples_indir,
	"-obs_stats_fn" => $epistat_fn, "-obs_pairs_fn" => $pairs_fn, "-pvalues_fn" => $pvalues_fn,
	"-pvalue_POThreshold" => $max_pvalue, "-stats_ext" => $epistat_ext, "-is_intragene" => $f_intragene, 
	"-pairs_ordering" => $f_ordered_pairs);
$smpl->print();