#!/usr/bin/env perl
#Make Benjamini-Hochberg correction procedure for a specified pvalue column in the dataframe
#Usage: [options] <table_filename>
#options:
#	-s uint - index of a column with p-values
#	-a ufloat - max accepted FDR value
#		if accounted the input table is splitted into accepted and rejected subsets

use strict;
use Getopt::Std;
use File::Basename;
use Statistics::Multtest qw(BH);

my %args;
if(!getopts('s:a:',\%args)){
	die "\nError in option string!";
}

my $idx=0;
if($args{s}=~m/^(\d+)/){
	$idx=$1-1;
	die "\nError: Columns need to be numbered starting from one!" if $idx<0;
}else{
	die "\nWrong value in the -s parameter: $args{s}!";
}
my $max_fdr;
if($args{a}=~m/^([0-9]*\.?[0-9]+)/){
	$max_fdr=$1;
	die "\nError: maximal FDR value is upper than one!" unless $max_fdr<1.0;
}
die "\nThe input file name is not accounted!" unless defined $ARGV[0];
my ($basename,$dir,$ext) = fileparse($ARGV[0],'\.[^\.]*$');
my @table;
my @pvals;
open INPF, "<$ARGV[0]" or die "\nUnable to open input file: $ARGV[0]!";
my $header;
my $sel_col_name;
my $i=0;
while(<INPF>){
	chomp;
	if(/\S+/){
		my @tmp=split "\t",$_,-1;
		unless(defined $header){
			$header=$_;
			$sel_col_name=$tmp[$idx];
		}else{
			if($tmp[$idx]<=1.0){
				push @pvals,$tmp[$idx];
				push @table,$_;
			}else{
				die "\nError in row $i in column ".($idx+1)." pvalue=$tmp[$idx]!";
			}
		}
	}
}
close INPF;
my $out_passed_fn;
my $out_failed_fn;
if(defined $max_fdr){
	$out_passed_fn=$dir.$basename;
	$out_failed_fn=$out_passed_fn;
	my $FDR=sprintf("%.0f",$max_fdr*100);
	$out_failed_fn.=".BH_less_FDR$FDR".$ext;
	$out_passed_fn.=".BH_greater_FDR$FDR".$ext;
	open OPF_PASS,">$out_passed_fn" or die "\nUnable to open output file $out_passed_fn!";
	open OPF_FAIL,">$out_failed_fn" or die "\nUnable to open output file $out_failed_fn!";
	print OPF_PASS "$header\tBH_adj.$sel_col_name";
	print OPF_FAIL "$header\tBH_adj.$sel_col_name";
}
my $adj_pvals=BH(\@pvals);
print "$header\tBH_adj.$sel_col_name";
for(my $i=0;$i<@table;$i++){
	print "\n$table[$i]\t".$adj_pvals->[$i];
	if(defined $max_fdr){
		if($adj_pvals->[$i]<$max_fdr){
			print OPF_FAIL "\n$table[$i]\t".$adj_pvals->[$i];
		}else{
			print OPF_PASS "\n$table[$i]\t".$adj_pvals->[$i];
		}
	}
}
if(defined $max_fdr){
	close OPF_FAIL;
	close OPF_PASS;
}