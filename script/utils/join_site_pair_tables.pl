#!/usr/bin/env perl
#Join two site pair tables
#Usage: [options] table1 table2
#options:
#	-s uint(,uint)* - Select specified columns from the table 2
#	-p STR - The prefix added to the colunm names in the table 2 
#	-n STR(,STR)* - Rename columns in the table 2
#		number of renamed columns need to be equal to the number of selected columns by -s option or to the total number of columns in the table 2

use strict;
use Getopt::Std;

my %args;
if(!getopts('s:p:n:',\%args)){
	die "\nError in option string!";
}
my $colname_prefix;
$colname_prefix=$args{p} if($args{p});

my @col_indices;
if($args{s}){
	@col_indices=split ",", $args{s};
	for(my $i=0;$i<@col_indices;$i++){
		$col_indices[$i]--;
	}
}
my @col_names;
if($args{n}){
	@col_names=split ",", $args{n};
	die "\nError: The number of colunm names is not equal to the number of selected columns in the table 2!" unless @col_names=@col_indices;
}

my %table2;
my $header2;
my $is_head;
open INPF, "<$ARGV[1]" or die "\nUnable to open input file $ARGV[1]!";
while(<INPF>){
	chomp;
	if(/\S+/){
		$is_head=1 unless defined $is_head;
		my @line=split "\t",$_,-1;
		my $str=$_;
		if(@col_indices){
			#select columns
			my @tmp;
			foreach my $i (@col_indices){
				push @tmp,$line[$i];
			}
			$str=join "\t",@tmp;
		}
		if($is_head){
			$header2=$str;
			$is_head=0;
		}else{
			if(/^(\d+\t\d+)\t/){
				$table2{$1}=$str;
			}
		}
	}
}
close INPF;
if(defined $colname_prefix || @col_names){
	my @tmp= split "\t", $header2;
	if(@col_names){
		die "\nError: The number of colunm names is not equal to the total number of columns in the table 2!" unless @col_names=@tmp;
		@tmp=@col_names;
	}
	if(defined $colname_prefix){
		for(my $i=0;$i<@tmp;$i++){
			$tmp[$i]=$colname_prefix.$tmp[$i];
		}
	}
	$header2=join "\t",@tmp;
}

my $is_head=undef;
open INPF, "<$ARGV[0]" or die "\nUnable to open input file $ARGV[0]!";
while(<INPF>){
	chomp;
	if(/\S+/){
		$is_head=1 unless defined $is_head;
		if($is_head){
			print $_."\t$header2";
			$is_head=0;
		}else{
			if(/^(\d+\t\d+)\t/){
				print "\n$_";
				print "\t".$table2{$1} if(defined $table2{$1});
			}
		}
	}
}
close INPF;