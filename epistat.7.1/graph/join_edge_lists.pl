#!/usr/bin/env perl
#This script joines two branch tables excluding duplicates
#Usage: <table1> <table2>
use strict;
my @tables_fn=@ARGV;
my %branches;
my $header;
foreach my $fn(@tables_fn){
	open INPF,"<$fn" or die "\nUnable to open input file $fn!";
	while(<INPF>){
		chomp;
		s/^\s+//;
		s/\s+$//;
		if(/\S+/){
			my @line=split '\t';
			my $str=$line[0]."\t".$line[1];
			if($line[0]=~/^\d+/ && $line[1]=~/^\d+/){
				my $w=$line[-1];
				if(defined $branches{$str}){
					warn "\nDifferent branch weights for the same site pairs: $str!" unless $branches{$str}==$w;
				}else{
					$branches{$str}=$w;
				}
			}else{
				$header=$str."\t".$line[-1] unless defined $header;
			}
		}
	}
	close INPF;
}
#print results
print $header;
foreach my $sp (keys %branches){
	print "\n$sp\t".$branches{$sp};
}
