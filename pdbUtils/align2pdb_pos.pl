#!/usr/bin/env perl
#The script converts alignment coordinates into positions on a primary sequence of PDB
#usage: options
#options:
#	-B <FN> - A name of file with a table for mapping of sites of marker sequence on positions of a primary sequence retrieved from ATOM records in the PDB file
#	-A <FN> - A name of file with the marker sequence taken from the initial alignment file as is (with gaps)
use strict;
use Getopt::Std;

my %args;
if(!getopts('B:A:',\%args)){
	die "\nError in option string!";
}
my $align_mseq_file=$args{A};
my $site2pdb_file=$args{B};

my @site2align_pos;
my $str;
open INPF, "<$align_mseq_file" or die "open <$align_mseq_file $!";
while (<INPF>){
	s/^\s+//;
	if(/^>/){
		$str="";
		last;
	}
}
die "\nThe file $align_mseq_file is not in the FASTA format!" unless defined $str;
while (<INPF>){
	last if /^>/;
	s/^\s+//;
	s/\s+$//;
	$str.=$_;
}
close INPF;
{
	my @seq=split "",$str;
	my $j=0;
	for(my $i=0;$i<@seq;$i++){
		push @site2align_pos,($i+1) unless $seq[$i] eq "-";
	}
}
my %site2pdb;
open INPF, "<$site2pdb_file" or die "open <$site2pdb_file $!";
while (<INPF>){
	chomp;
	s/^\s+//;
	s/\s+$//;
	my @splitter = split /\s+/;
	if(/^(\d+)\s+(\d+)/){
		$site2pdb{$site2align_pos[$1-1]}=$2 if ($1>0)&&($2>0);
		unless($1<=@site2align_pos){
			my $str=@site2align_pos;
			die "\nThe marker sequence Nsites=$str is inconsistent to the sequence used for mapping on the PDB (pos=$1)!";
		}
	}
}
close INPF;
print "align_pos\tpdb_pos\n";
foreach my $apos(@site2align_pos){
	print "$apos\t$site2pdb{$apos}\n";
}