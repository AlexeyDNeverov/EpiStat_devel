#!/usr/bin/env perl
#This script draws mutations on tree branches
#options:
#	-x <FN> - input file for EpiStat application in XPAR format [Kryazhimsky11]. Contains tree in the adjacent file format
#	-p <epistat_prm> - input file with parameters for running epistat.pl script
#	[-u] unordered site pairs. Print both ordered site pairs on one tree
#Params:
#<site_pairs_list> - a file with list of site pairs in TAB delimited format. Mutations for each pair are drawn on separate tree

use strict;
use Getopt::Std;
use File::Basename;
die "Set up the system variable \$EPISTAT_HOME" unless defined $ENV{EPISTAT_HOME};
my $epistat_cmd="$ENV{EPISTAT_HOME}/epistat.pl";

my %args;
if(!getopts('x:p:u',\%args)){
	die "\nError in option string!";
}

my $xparr_fn=$args{x};
die "\nNo default parameter for -x option!" unless defined $xparr_fn;
my $epistat_prm=$args{p};
die "\nNo default parameter for -p option!" unless defined $epistat_prm;
my $site_pairs_fn=$ARGV[0];
die "\nThe file with a list of site pairs to draw is not accounted!" unless defined $site_pairs_fn;

sub gen_tempname{
	my $nchar=shift;
	my @chars = ( "A" .. "Z", "a" .. "z", 0 .. 9 );
	return join("", @chars[ map { rand @chars } ( 1 .. $nchar ) ]);
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

my @tmp_fn;
open INPF, "<$site_pairs_fn" or die "\nUnable to open input file: $site_pairs_fn!";
my %site_pairs;
while(<INPF>){
	my @line=split '\t';
	$line[0]=~s/\s+//;
	$line[1]=~s/\s+//;
	if($line[0]=~/\d+/ || $line[1]=~/\d+/){
		my ($bgs,$tgs)=($line[0],$line[1]);
		my $str;
		if(defined $args{u}){
			if($bgs>$tgs){
				$str="$tgs,$bgs";
			}else{
				$str="$bgs,$tgs";
			}
		}else{
			$str="$bgs,$tgs";
		}
		$site_pairs{$str}++;
	}
}
close INPF;
my $tmpfiles_fn="sp_temp_list.txt";
open OUTF, ">$tmpfiles_fn" or die "\Unable to open output file: $tmpfiles_fn!";
my @sp_keys=keys %site_pairs;
foreach my $sp(@sp_keys){
	my $tmp_fn=gen_tempname(10);
	push @tmp_fn,$tmp_fn;
	$tmp_fn.=".csv";
	print OUTF "$tmp_fn\n";
	open OPF, ">$tmp_fn" or die "\nUnable to open output file: $tmp_fn!";
	print OPF $sp;
	if(defined $args{u}){
		my ($bgs,$tgs)=split ",",$sp;
		print OPF "\n$tgs,$bgs";
	}
	close OPF;
}
close OUTF;

my $str="cat $tmpfiles_fn|parallel $epistat_cmd -x $xparr_fn -a figtree_site_pairs -s {} $epistat_prm";
print STDERR "\n\nPrinting trees!\n$str";
system($str);
print_child_termination_status();

my ($basename,$dir,$ext) = fileparse($site_pairs_fn,'\.[^.]*$');
my $outfn=$dir.$basename.".mutations.tre";
open OUTF, ">$outfn" or die "\Unable to open output file: $outfn!";
print OUTF "#NEXUS\nbegin trees;";
for(my $i=0;$i<@tmp_fn;$i++){
	my $tmp_fn=$tmp_fn[$i];
	open INPF, "<$tmp_fn.site_pairs.tre" or warn "\nUnable to open input file $tmp_fn.tre!";
	while(<INPF>){
		chomp;
		if(/^\ttree\s+(\w+)\s*=\s*/){
			my $tree_str=$';
			my $tree_name=$sp_keys[$i];
			$tree_name=~s/,/_/;
			my $j=$i+1;
			$tree_name.="_n".$j;
			print OUTF "\n\ttree $tree_name = $tree_str";
		}
	}
	close INPF;
	unlink glob $tmp_fn.".*";
}
print OUTF "\nend;";
close OUTF;
unlink $tmpfiles_fn;


