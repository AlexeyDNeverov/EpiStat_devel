#!/usr/bin/env perl
#add FDR estimates for pvalues into the table
#Usage: options table_file_name
#options:
#	-f <FN> - FDR estimates file name
#	-p uint - pvalue column number
use strict;
use Getopt::Std;
use List::BinarySearch qw(binsearch_pos);

my %args;
if(!getopts('f:p:',\%args)){
	die "\nError in option string!";
}
my $pval_column=$args{p};
my $fdr_fn=$args{f};
die "\nThe number of the column containing pvalues is required!" unless defined $pval_column;
$pval_column--;
die "\nThe file with FDR estimates for pvalues is required!" unless defined $fdr_fn;
my @pvalues;
my @fdr;
open INPF, "<$fdr_fn" or die "\nUnable to open input file $fdr_fn!";
while(<INPF>){
	chomp;
	if(/^\d+([.Ee][+-]?\d+)?\t/){
		my @line=split "\t";
		if(@line==4){
			push @pvalues,$line[0];
			push @fdr,[@line[1 .. 3]];
		}
	}
}
close INPF;
open INPF, "<$ARGV[0]"  or die "\nUnable to open input file $ARGV[0]!";
my $is_head;
while(<INPF>){
	chomp;
	if(/\S+/){
		$is_head=1 unless defined $is_head;
		my @line=split "\t",$_,-1;
		print join("\t",@line[0 .. $pval_column]);
		if($is_head){
			print "\tFDR\t";
			$is_head=0;
		}else{
			my $i = binsearch_pos { $a <=> $b } $line[$pval_column], @pvalues;
			die "\nUnable to find a pvalue=$line[$pval_column] in the $fdr_fn!" if $i==@pvalues;
			print "\t";
			printf("%.3f", $fdr[$i]->[1]/$fdr[$i]->[0]);
			print "\t";
		}
		print join("\t",@line[$pval_column+1 .. $#line]) if @line-$pval_column-1;
		print "\n";
	}
}
close INPF;