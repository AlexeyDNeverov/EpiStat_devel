#!/usr/bin/env perl
use strict;
open INPF, "<$ARGV[0]" or die "\nUnable to open input file $ARGV[0]!";
my %subst_map;
while(<INPF>){
	chomp;
	my @line=split "\t",$_,-1;
	if(defined($line[0])&&defined($line[1])){
		if(defined($line[4])){
			my $pname=$line[1];
			$subst_map{$pname}={} unless defined $subst_map{$pname};
			my @subst=split ";",$line[4];
			foreach my $ss(@subst){
				if($ss=~/[A-Z](\d+)[A-Z]/){
					my $site=$1;
					$subst_map{$pname}->{$ss}++;
				}
			}
		}
	}
}
close INPF;
foreach my $name(keys %subst_map){
	foreach my $ss(keys %{$subst_map{$name}}){
		print "\n$name\t$ss\t".$subst_map{$name}->{$ss} if $subst_map{$name}->{$ss}>1;
	}
}