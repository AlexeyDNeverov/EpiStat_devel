#!/usr/bin/env perl
#Join two site pair tables
#Usage: [options] table1 table2
#options:
#	-s uint(,uint)* - Select specified columns from the table 2
use strict;
use Getopt::Std;

my %args;
if(!getopts('s:',\%args)){
	die "\nError in option string!";
}
my @col_indices;
if($args{s}){
	@col_indices=split ",", $args{s};
	for(my $i=0;$i<@col_indices;$i++){
		$col_indices[$i]--;
	}
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