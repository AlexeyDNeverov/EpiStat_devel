#!/usr/bin/env perl
#This script calculates pvalue statistics for a fake sample with specified ID using samples from NULL model
#options:
#	-d <PATH> - input folder with samples
#	-p <uint> Total number of samples
#	-i <uint> Sample ID, must be in range 1..<total number of samples>
#	[-s] <FN> Site pairs file name is required to calculate pvalues for unordered site pairs
#	[-m] <uint> - Maximal number of samples in storage folders. 
#		If it is less than total number of samples the input folder is a root of tree of folders storing fake samples' files. 
#Params:
#<parameters> - epistat.pl configuration file
#<parameters> - configuration file
#		Contents:
#		[FDRSampleNum="<uint>"] - Number of samples used to calculate FDR
#			Default="0" - Do not calculate FDR
#		[TMPClean="0|1"] - Remove after finishing all temporary files and folders
#			Default="0" - Do not remove temp files
#		OutUpperPValues="Suffix" - Suffix used for a pvalues' output file
#		OutLowerPValues="Suffix" - Suffix used for a pvalues' output file
#		OutAvgEpiStats="Suffix" - Suffix used for an average epistatic statistics output file
#		OutVarEpiStat="Suffix" - Suffix used for a variance of epistatic statistics output file
#		OutBgrSiteStats="Suffix" - Suffix used for an output file with average and variance of marginal epistatic statistics'  values for background sites
#		OutFgrSiteStats="Suffix" - Suffix used for an output file with average and variance of marginal epistatic statistics'  values for foreground sites
#		[OutCovOrdPairsEpiStat="Suffix"] - Suffix used for output file with covariances of epistatic statistics of two ordered site pairs. Covariances are calculated for internal epistasis only.
#			Default=".ord_site_pars.cov"
#		[OutFDR="Suffix"] - Suffix used for files with statistics for FDR samples

use strict;
use lib "$ENV{EPISTAT_LIB}";

use Getopt::Std;
use EpistatNullModel::General;
use EpistatNullModel::Symmetric;
use IndexedFakeSamplesHeap;

my %args;
if(!getopts('d:p:i:s:m:',\%args)){
	die "\nError in option string!";
}
my $max_items_indir=10000000;
my $files_heap;
my $outdir=$args{d};
die "\nThe input folder not specified! Use -d option!" unless defined $outdir;
$outdir=~s/\/$//;
my $nperm=$args{p};
die "\nThe the number of samples not specified! Use -p option!" unless defined $nperm;
if(defined $args{m}){
	if($args{m}=~/^[1-9]\d*/){
		$max_items_indir=$args{m};
	}else{
		die "\nIllegal value for the '-m' option!";
	}
}
my $sample_id=$args{i};
die "\nThe input sample not specified! Use -i option!" unless defined $sample_id;
die "\nSample ID must be not greater than total number of samples $nperm!" unless $sample_id<=$nperm;
my $pairs_fn=$args{s};
my $stats_ext;
my $xparr_ext=".xpar";
my $tau;
my $f_intragene;
my $f_stat_type;
open INFILE, "$ARGV[0]" or die "\nUnable to open input file: $ARGV[0]!";
while(<INFILE>){
	$_=$` if(/#/);
	if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
		my ($key,$value)=($1,$2);
		if($key eq "OutEpiStatistics"){
			$stats_ext=$value;
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
			};
		}elsif($key eq "PrintAllPairs"){
			die "\nUnknown value ($value) for the PrintAllPairs parameter in file $ARGV[0]!\n\tUINT in range is [0..1] expected!" unless $value=~/^[01]\s*$/;
			die "\nThe PrintAllPairs=\"0\" in $ARGV[0] is not supported!\n\tUse PrintAllPairs=\"1\"" unless $value==1;
		}
	}
}
close INFILE;
die "\nError: Parameter PairsType is undefined in the file $ARGV[0]!" unless defined $f_intragene;

unless(defined $stats_ext){
	$stats_ext=".stat.exp";
	$stats_ext.=".$tau" if defined $tau;
}
my ($upper_pval_fn,$lower_pval_fn,$avg_fn,$var_fn,$fdr_fn);
my $cov_fn=".ord_site_pars.cov";
open INFILE, "$ARGV[1]" or die "\nUnable to open input file: $ARGV[1]!";
while(<INFILE>){
	$_=$` if(/#/);
	if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
		my ($key,$value)=($1,$2);
		if($key eq "OutUpperPValues"){
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
		}
	}
}
close INFILE;
die "\nError: The OutUpperPValues parameter has no defaul value!" unless defined $upper_pval_fn;
die "\nError: The OutLowerPValues parameter has no defaul value!" unless defined $lower_pval_fn;
die "\nError: The OutAvgEpiStats parameter has no defaul value!" unless defined $avg_fn;
die "\nError: The OutVarEpiStat parameter has no defaul value!" unless defined $var_fn;
die "\nError: The OutFDR parameter has no defaul value!" unless defined($fdr_fn);

my @out_file_types;
{
	my $pvalue_mode_str="";
	if($f_intragene){
		$pvalue_mode_str=".intragene.unord_pairs";
		push @out_file_types,$pvalue_mode_str.$lower_pval_fn;
		push @out_file_types,$pvalue_mode_str.$upper_pval_fn;
		$pvalue_mode_str=".intragene.ord_pairs";
	}
	unless(($f_stat_type==3)&&$f_intragene){
		push @out_file_types,$pvalue_mode_str.$lower_pval_fn;
		push @out_file_types,$pvalue_mode_str.$upper_pval_fn;
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

if($nperm){
	my $fdrdir=$outdir."/"."fdr";
	my $samples_dir=$outdir."/samples";
	my $name=$fdrdir."/".$sample_id.$fdr_fn;
	my $pvalue_mode_str="";
	unless(check_outfiles($name,@out_file_types)){
		my @stat1=();
		my @stat2=();
		my $h0;
		my $sample_id_path=$samples_dir;
		unless(($f_stat_type==3)&&$f_intragene){
			if($nperm>$max_items_indir){
				my $samples_heap=IndexedFakeSamplesHeap->new($nperm,$max_items_indir);
				$sample_id_path.="/".$samples_heap->sample_idx2path($sample_id);
				$h0=EpistatNullModel::General->new("-size" => $nperm,"-samples_dir" => $samples_dir,
						"-stats_ext" => $stats_ext,"-obs_pairs_fn" => $pairs_fn,
						"-skip_obs_id" => $sample_id, "-is_intragene" => $f_intragene, 
						"-max_samples_indir" => $max_items_indir);
			}else{
				$h0=EpistatNullModel::General->new("-size" => $nperm,"-samples_dir" => $samples_dir,
						"-stats_ext" => $stats_ext,"-obs_pairs_fn" => $pairs_fn,
						"-skip_obs_id" => $sample_id, "-is_intragene" => $f_intragene);
			}
			if($f_intragene){
				$pvalue_mode_str=".intragene.ord_pairs";
			}
			$h0->calc_pvalues(\@stat1,\@stat2);
			print_statistics($name.$pvalue_mode_str.$lower_pval_fn,undef,\@stat1);
			print_statistics($name.$pvalue_mode_str.$upper_pval_fn,undef,\@stat2);
			@stat1=();
			@stat2=();
		}else{
			if($nperm>$max_items_indir){
				my $samples_heap=IndexedFakeSamplesHeap->new($nperm,$max_items_indir);
				$sample_id_path.="/".$samples_heap->sample_idx2path($sample_id);
				$h0=EpistatNullModel::Symmetric->new("-size" => $nperm,"-samples_dir" => $samples_dir,
						"-stats_ext" => $stats_ext,"-obs_pairs_fn" => $pairs_fn,
						"-skip_obs_id" => $sample_id,"-max_samples_indir" => $max_items_indir);
			}else{
				$h0=EpistatNullModel::Symmetric->new("-size" => $nperm,"-samples_dir" => $samples_dir,
						"-stats_ext" => $stats_ext,"-obs_pairs_fn" => $pairs_fn,
						"-skip_obs_id" => $sample_id);
			}
		}
		if($f_intragene){
			$h0->calc_unordered_pairs_pvalues(\@stat1,\@stat2);
			$pvalue_mode_str=".intragene.unord_pairs";
			print_statistics($name.$pvalue_mode_str.$lower_pval_fn,undef,\@stat1);
			print_statistics($name.$pvalue_mode_str.$upper_pval_fn,undef,\@stat2);
			@stat1=();
			@stat2=();
		}
		undef($h0);
		
		if(!(-e $fdrdir."/".$sample_id.$stats_ext)){
			my $str="cp $sample_id_path/$sample_id$stats_ext $fdrdir";
			system($str);
			print_child_termination_status();
		}
		if((-e "$sample_id_path/$sample_id$xparr_ext")&&(!(-e $fdrdir."/".$sample_id.$xparr_ext))){
			my $str="cp $sample_id_path/$sample_id$xparr_ext $fdrdir";
			system($str);
			print_child_termination_status();
		}
	}
}
