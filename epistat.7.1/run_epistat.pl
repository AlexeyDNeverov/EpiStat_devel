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
#	[-P] <FN> - input file with options string for GNU Parallel
#	[-J] <uint> - Number of parallel jobs for running 'shuffle_incidence_matrix.R' script
#		Default=0 auto defined
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
#		OutBgrSiteStats="Suffix" - Suffix used for an output file with average and variance of marginal epistatic statistics'  values for background sites
#		OutFgrSiteStats="Suffix" - Suffix used for an output file with average and variance of marginal epistatic statistics'  values for foreground sites
#		[OutCovOrdPairsEpiStat="Suffix"] - Suffix used for output file with covariances of epistatic statistics of two ordered site pairs. Covariances are calculated for internal epistasis only.
#			Default=".ord_site_pars.cov"
#		[OutFDR="Suffix"] - Suffix used for files with statistics for FDR samples

use strict;
#server inst
use lib "$ENV{EPISTAT_HOME}";
use Getopt::Std;
use File::Basename;
use File::Path qw(make_path rmtree);
use EpistatNullModel;
use IndexedFakeSamplesHeap;
use Phenotype::AggregateStatistics;

#server inst
my $rcmd="Rscript";
my $brw_R="$ENV{EPISTAT_HOME}/shuffle_incidence_matrix.R";
my $epistat_cmd="$ENV{EPISTAT_HOME}/epistat.pl";
my $sample2xparr_cmd="$ENV{EPISTAT_HOME}/sample2xparr.pl";
my $fake_sample_cmd="$ENV{EPISTAT_HOME}/fake_sample.pl";

my %args;
if(!getopts('x:m:p:R:rJ:P:',\%args)){
	die "\nError in option string!";
}
my $max_items_indir=1000;
my $no_brw_jobs=0;
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
my $f_intragene;
my $f_ignore_terminals=1;
my $f_rerun;
my $ntries=1;
#for phenotypes:
my @pheno_labels;
my $npheno;
#specifyes the method of phenotypes aggregation: 0 - best_association, 1 - mean_association
my $f_pheno_aggr_method;
my $bg_sites2pheno_fn;
my $fg_sites2pheno_fn;
#0- ascending; 1- descending
my $site2phen_order;

$f_rerun=1 if defined $args{r};
if($args{R}=~/^\d+/){
	$ntries=$args{R}+1;
}else{
	die "\nWrong value for an option -R, UINT is expected!" if defined $args{R};
}
if(defined $args{J}){
	die "\nError wrong option's value: -J ".$args{J}."!" unless $args{J}=~/^\d+/;
	$no_brw_jobs=$args{J};
}
my $run_parallel_opts="";
if(defined $args{P}){
	open INPF, $args{P} or die "\nUnable to open input file: $args{P}!";
	while(<INPF>){
		chomp;
		if(/^\s*-{1,2}\S/){
			$run_parallel_opts=$_;
			$run_parallel_opts=~s/^\s+//;
			$run_parallel_opts=~s/\s+$//;
			$run_parallel_opts.=" ";
			last;
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
		}elsif($key eq "OutEpiStatistics"){
			$stats_fn=$value;
		}elsif($key eq "IgnoreTerminals"){
			$f_ignore_terminals=$value;
			die "\nUnknown value ($value) for the IgnoreTerminals parameter in file $ARGV[0]!\n\tUINT in range is [0..1] expected!" unless $value=~/^[01]\s*$/;		}elsif($key eq "OutMutMatrix"){
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
			}
		}elsif($key eq "TAU"){
			if($value=~/\d+(\.\d+)*/){
				$tau=$value;
			}else{
				die "\nWrong value for the TAU parameter in file $ARGV[0]!\n\tPositive number expected!";
			}
		}elsif($key eq "Phenotypes"){
			if($value==1){
				$npheno=0;
			}else{
				die "\nUnpropper value for the 'Phenotype' parameter! Valid values: 0 or 1" unless $value==0;
			}
		}elsif($key eq "PrintAllPairs"){
			die "\nUnknown value ($value) for the PrintAllPairs parameter in file $ARGV[0]!\n\tUINT in range is [0..1] expected!" unless $value=~/^[01]\s*$/;
			die "\nThe PrintAllPairs=\"0\" in $ARGV[0] is not supported!\n\tUse PrintAllPairs=\"1\"" unless $value==1;
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
my ($bgr_stat_fn,$fgr_stat_fn);
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
		}elsif($key eq "OutBgrSiteStats"){
			$bgr_stat_fn=$value;
		}elsif($key eq "OutFgrSiteStats"){
			$fgr_stat_fn=$value;
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
		}else{
			die "\nUnknown parameter: $key in the input file $ARGV[1]!";
		}
	}
}
close INFILE;
die "\nError: The OutUpperPValues parameter has no defaul value!" unless defined $upper_pval_fn;
die "\nError: The OutLowerPValues parameter has no defaul value!" unless defined $lower_pval_fn;
die "\nError: The OutAvgEpiStats parameter has no defaul value!" unless defined $avg_fn;
die "\nError: The OutVarEpiStat parameter has no defaul value!" unless defined $var_fn;
die "\nError: The OutFDR parameter has no defaul value!" if $nfdr_samples>0&&!defined($fdr_fn);
die "\nError: The OutBgrSiteStats parameter has no defaul value!" unless defined $bgr_stat_fn;
die "\nError: The OutFgrSiteStats parameter has no defaul value!" unless defined $fgr_stat_fn;
if(defined($npheno)&&($npheno>1)){
	die "\nUndefined aggregation method for phenotype-constrained fake samples!" unless defined $f_pheno_aggr_method;
	die "\nNo files with site to phenotype associations has been specified!" unless defined($bg_sites2pheno_fn)||defined($fg_sites2pheno_fn);
	die "\nThe best to worst sorting order for site to phenotype associations is required!" unless defined $site2phen_order;
	unless($f_intragene){
		die "\nThe file with association statistics of sites in background to phenotype hasn't been specified!" unless defined($bg_sites2pheno_fn);
		die "\nThe file with association statistics of sites in foreground to phenotype hasn't been specified!" unless defined($fg_sites2pheno_fn);
	}
}

my ($basename,$dir,$ext) = fileparse($xpar_fn,'\.[^\.]*$');
$pairs_fn=$dir.$basename.$pairs_ext;
$stats_fn=$dir.$basename.$stats_ext;
if(!$f_rerun){
	#Deleting previous versions of files
	unlink $pairs_fn if(-e $pairs_fn);
	unlink $stats_fn if(-e $stats_fn);
}

#0-matrix; 1-idx2branch; 2-idx2site
my @bg_files;
my @tg_files;
my @mmtx_prefs;
{
	if(defined $npheno){
		for(my $i=0;$i<$npheno;$i++){
			push @mmtx_prefs,[($pheno_labels[$i],"#0")];
			push @mmtx_prefs,[($pheno_labels[$i],"#1")];
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
		if($f_mtx_mode!=2){
			my $str;
			if($f_intragene&&$f_mtx_mode==3){
				$str.=".both";
			}else{
				$str.=".bgr";
			}
			push @bg_files,[($fn.$str.$mmtx_fn,$fn.$str.$idx2branch_fn,$fn.$str.$idx2site_fn)];
		}
		if($f_mtx_mode!=1){
			if(!($f_intragene&&$f_mtx_mode==3)){
				push @tg_files,[($fn.".fgr".$mmtx_fn,$fn.".fgr".$idx2branch_fn,$fn.".fgr".$idx2site_fn)];
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
};

my $str=$epistat_cmd." -x $xpar_fn -m $f_mtx_mode $ARGV[0]";
unless($f_rerun&&(-e $pairs_fn)&&(-e $stats_fn)){
	print STDERR "Calculating epistatic statistics for the input tree: $xpar_fn\n$str";
	system($str);
	print_child_termination_status();
}
die "\n\tFatal Error: The expected file is not exist: $pairs_fn" unless -e $pairs_fn;
die "\n\tFatal Error: The expected file is not exist: $stats_fn" unless -e $stats_fn;
my @out_file_types;
my @out_fake_file_types;
{
	my $pvalue_mode_str="";
	if($f_intragene){
		$pvalue_mode_str=".intragene.unord_pairs";
		push @out_file_types,$pvalue_mode_str.$lower_pval_fn;
		push @out_file_types,$pvalue_mode_str.$upper_pval_fn;
		$pvalue_mode_str=".intragene.ord_pairs";
		@out_fake_file_types=@out_file_types;
		push @out_file_types,$cov_fn;
	}
	push @out_file_types,$pvalue_mode_str.$lower_pval_fn;
	push @out_file_types,$pvalue_mode_str.$upper_pval_fn;
	push @out_fake_file_types,$pvalue_mode_str.$lower_pval_fn;
	push @out_fake_file_types,$pvalue_mode_str.$upper_pval_fn;
	push @out_file_types,$avg_fn;
	push @out_file_types,$var_fn;
	push @out_file_types,$bgr_stat_fn;
	push @out_file_types,$fgr_stat_fn;
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
	my ($infile,$ra_pheno_labels,$rh_sites2pheno)=@_;
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
			$rh_sites2pheno->{$line[0]}=[@line[1..$npheno]];
			for(my $j=1;$j<@line;$j++){
				die "\nUnable to use negative association statistics: site=$line[0],phenotype=$j" if $line[$j]<0;
			}
			$I++;
		}
	}
	close INPF;
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
		for(my $ids=0;$ids<@mmtx_prefs;$ids++){
			my $outdir=$dir.$basename;
			$outdir.="/pheno/".join "/",@{$mmtx_prefs[$ids]} if(defined $npheno);
			if(@bg_files){
				my $n=0;
				if($f_rerun&&(-d $outdir."/"."bgr")){
					$n=count_samples($outdir."/bgr",$perm_fn,\@samples_storage);
					rmtree($outdir."/"."bgr") unless($n==$nperm);
				}
				unless(-d $outdir."/"."bgr"){
					my $t=1;
					foreach my $fn(@{$bg_files[$ids]}){
						if(! -e $fn){
							warn "\n\tFatal Error: The expected file is not exists: $fn";
							$t=0;
						}
					}
					die unless $t;
					for(my $i=0;$i<@samples_storage;$i++){
						my $from=$samples_storage[$i]->[0];
						my $to=$samples_storage[$i]->[1];
						my $mmtx_fn=$bg_files[$ids]->[0];
						my $out_path=$outdir."/bgr".$samples_storage[$i]->[2];
						if($no_brw_jobs==1){
							$str="$rcmd $brw_R $mmtx_fn $out_path $from $to";
						}else{
							$str=$run_parallel_opts;
							$str.="-j $no_brw_jobs " if $no_brw_jobs;
							$str="seq $from $to | parallel ".$str."$rcmd $brw_R $mmtx_fn $out_path {}";
						}
						print STDERR "\n\nGenerating permutations for background\n$str";
						make_path($out_path);
						system($str);
						print_child_termination_status();
					}
					$n=count_samples($outdir."/bgr",$perm_fn,\@samples_storage);
					die "\nError: There are less than expected *.$perm_fn files in: $outdir/bgr!" unless($n==$nperm);
				}
			}
			if(@tg_files){
				my $n=0;
				if($f_rerun&&(-d $outdir."/"."fgr")){
					$n=count_samples($outdir."/fgr",$perm_fn,\@samples_storage);
					rmtree($outdir."/"."fgr") unless($n==$nperm);
				}
				unless(-d $outdir."/"."fgr"){
					my $t=1;
					foreach my $fn(@{$tg_files[$ids]}){
						if(! -e $fn){
							warn "\n\tFatal Error: The expected file is not exists: $fn";
							$t=0;
						}
					};
					die unless $t;
					for(my $i=0;$i<@samples_storage;$i++){
						my $from=$samples_storage[$i]->[0];
						my $to=$samples_storage[$i]->[1];
						my $mmtx_fn=$tg_files[$ids]->[0];
						my $out_path=$outdir."/fgr".$samples_storage[$i]->[2];
						if($no_brw_jobs==1){
							$str="$rcmd $brw_R $mmtx_fn $out_path $from $to";
						}else{
							$str=$run_parallel_opts;
							$str.="-j $no_brw_jobs " if $no_brw_jobs;
							$str="seq $from $to | parallel ".$str."$rcmd $brw_R $mmtx_fn $out_path {}";
						}
						print STDERR "\n\nGenerating permutations for foreground\n$str";
						make_path($out_path);
						system($str);
						print_child_termination_status();
					}
					$n=count_samples($outdir."/fgr",$perm_fn,\@samples_storage);
					die "\nError: There are less than expected *$perm_fn files in: $outdir/fgr!" unless($n==$nperm);
				}
			}
		}
		#die("\nWrong number of file sets!") unless scalar(@bg_files)==scalar(@tg_files)||scalar(@bg_files)==0||scalar(@tg_files)==0;
		my $step=1;
		$step=2 if defined $npheno;
		for(my $fsi=0;$fsi<@mmtx_prefs;$fsi+=$step){
			my $sample2xparr_prm=$outdir;
			$sample2xparr_prm.=".".$mmtx_prefs[$fsi]->[0] if defined $npheno;
			$sample2xparr_prm.=".sample2xparr.prm";
			open OPF, ">$sample2xparr_prm" or die "\nUnable to open output file:$sample2xparr_prm!";
			print OPF "XPAR=\"$xpar_fn\"";
			print OPF "\nPairsType=";
			if($f_intragene){
				print OPF "\"intragene\"";
			}else{
				print OPF "\"intergene\"";
			};
			print OPF "\nMatrixMode=\"$f_mtx_mode\"";
			print OPF "\nIgnoreTerminals=\"$f_ignore_terminals\"";
			print OPF "\nNBranchSubsets=\"$step\"";
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
				if(@bg_files>0){
					$str.="-b bgr".$out_path."/{}$perm_fn";
				}
				if(@tg_files>0){
					$str.=" " if(@bg_files>0);
					$str.="-t fgr".$out_path."/{}$perm_fn";
				};
				if($f_rerun){
					@gap=find_missed_samples($samples_dir,".xpar",$samples_storage[$i]);
					if(@gap){
						$tmp_fname=gen_tempname(10).".xpar.missed";
						dump_into(\@gap,$tmp_fname);
						$cmd_str="cat $tmp_fname";
					}
				}else{$cmd_str="seq $from $to";};
				if($cmd_str){
					$str=$cmd_str." | parallel ".$run_parallel_opts."$sample2xparr_cmd $str $sample2xparr_prm"."\">\"".$samples_dir.$out_path."/{}.xpar";
					print STDERR "\n\nGenerating XPAR files for randomized trees!\n$str";
					system($str);
					print_child_termination_status();
					unlink $tmp_fname if defined $tmp_fname;
				}
			}
			my $n=count_samples($samples_dir,".xpar",\@samples_storage);
			die "\nError: There are less than expected *.xpar files in: $samples_dir!" unless $n==$nperm;
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
					$str=$cmd_str." | parallel ".$run_parallel_opts."$epistat_cmd $str ".$ARGV[0];
					print STDERR "\n\nCalculating epistatic statistics for randomized trees!\n$str";
					system($str);
					print_child_termination_status();
					unlink $tmp_fname if defined $tmp_fname;
				}
			}
			$n=count_samples($samples_dir,$stats_ext,\@samples_storage);
			die "\nError: There are less than expected *$stats_ext files in: $samples_dir!" unless $n==$nperm;
		}
		if(defined $npheno){
			if((!$f_rerun)&&(-d $samples_dir)){
				rmtree($samples_dir);
			}
			my @pheno_samples_dirs=@pheno_labels;
			for(my $i=0;$i<$npheno;$i++){
				$pheno_samples_dirs[$i].="/samples";
			}
			if($npheno>1){
				my $rh_bg_sites2pheno;
				my $rh_fg_sites2pheno;
				if(defined $bg_sites2pheno_fn){
					$rh_bg_sites2pheno={};
					read_sites2pheno($bg_sites2pheno_fn,\@pheno_labels,$rh_bg_sites2pheno);
				}
				if(defined $fg_sites2pheno_fn){
					unless($f_intragene&&($bg_sites2pheno_fn eq $fg_sites2pheno_fn)){
						$rh_fg_sites2pheno={};
						read_sites2pheno($fg_sites2pheno_fn,\@pheno_labels,$rh_fg_sites2pheno);
					}
				}
				if($f_intragene){
					$rh_fg_sites2pheno=$rh_bg_sites2pheno unless defined $rh_fg_sites2pheno;
					$rh_bg_sites2pheno=$rh_fg_sites2pheno unless defined $rh_bg_sites2pheno;
				}
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
			die "\nError: There are less than expected *$stats_ext files in: $samples_dir!" unless $n==$nperm;
		}
		my $nlines;
		my @stat1;
		my @stat2;
		print STDERR "\n\nBuilding zero model distribution and calculating statistics!";
		my $h0;
		if($nperm>$max_items_indir){
			$h0=EpistatNullModel->new("-size" => $nperm,"-samples_dir" => $samples_dir,"-stats_ext" => $stats_ext,
				"-obs_stats_fn" => $stats_fn,"-obs_pairs_fn" => $pairs_fn, 
				"-is_intragene" => $f_intragene, "-max_samples_indir" => $max_items_indir);
		}else{
			$h0=EpistatNullModel->new("-size" => $nperm,"-samples_dir" => $samples_dir,"-stats_ext" => $stats_ext,
				"-obs_stats_fn" => $stats_fn,"-obs_pairs_fn" => $pairs_fn, "-is_intragene" => $f_intragene);
		}
		{#Print marginal epistatic statistics for sites
			my @bgr_sites;
			my @fgr_sites;
			$h0->get_sites(\@bgr_sites,\@fgr_sites);
			$h0->calc_bgr_moments12(\@stat1,\@stat2);
			$str=$outdir.$bgr_stat_fn;
			print_statistics($str,"site\tavg\tvar",\@bgr_sites,\@stat1,\@stat2);
			$h0->calc_fgr_moments12(\@stat1,\@stat2);
			$str=$outdir.$fgr_stat_fn;
			print_statistics($str,"site\tavg\tvar",\@fgr_sites,\@stat1,\@stat2);
		}
		$h0->calc_moments12(\@stat1,\@stat2);
		$str=$outdir.$avg_fn;
		print_statistics($str,undef,\@stat1);
		$str=$outdir.$var_fn;
		print_statistics($str,undef,\@stat2);
		if($f_intragene){
			$h0->calc_ordered_site_pairs_covariances(\@stat1);
			$str=$outdir.$cov_fn;
			print_statistics($str,undef,\@stat1);
			$h0->calc_unordered_pairs_pvalues(\@stat1,\@stat2);
			$pvalue_mode_str=".intragene.unord_pairs";
			$str=$outdir.$pvalue_mode_str.$lower_pval_fn;
			print_statistics($str,undef,\@stat1);
			$str=$outdir.$pvalue_mode_str.$upper_pval_fn;
			print_statistics($str,undef,\@stat2);
			$pvalue_mode_str=".intragene.ord_pairs";
		}
		$h0->calc_pvalues(\@stat1,\@stat2);
		$str=$outdir.$pvalue_mode_str.$lower_pval_fn;
		print_statistics($str,undef,\@stat1);
		$str=$outdir.$pvalue_mode_str.$upper_pval_fn;
		print_statistics($str,undef,\@stat2);
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
		};
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
		$str="cat $tmp_fname | parallel ".$run_parallel_opts."$fake_sample_cmd $str ".$ARGV[0]." ".$ARGV[1];
		for(my $I=0;$I<$ntries;$I++){
			print STDERR "\n\nTry$I of FDR estimation for $nfdr_samples samples in parallel:!\n$str";
			system($str);
			print_child_termination_status();
			$n=0;
			for(my $j=0;$j<$nfdr_samples;$j++){
				$k=$fdr_samples_keys[$j];
				$n++ if -e $fdrdir."/".$k.$stats_ext && check_outfiles($fdrdir."/".$k,@out_fake_file_types);
			}
			last if $n==$nfdr_samples;
		}
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
