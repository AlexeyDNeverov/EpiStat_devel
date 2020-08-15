#!/usr/bin/env perl
#This script draws mutations on tree branches
#options:
#	-x <FN> - input file for EpiStat application in XPAR format [Kryazhimsky11]. Contains tree in the adjacent file format
#	-p <epistat_prm> - input file with parameters for running epistat.pl script
#Params:
#<site_pairs_list> - CSV file with list of site pairs. Mutations for each pair are drawn on separate tree

use strict;
use Getopt::Std;
use File::Basename;

my $epistat_cmd="$ENV{EPISTAT_HOME}/epistat.pl";

my %args;
if(!getopts('x:p:',\%args)){
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
my $tmpfiles_fn="sp_temp_list.txt";
open OUTF, ">$tmpfiles_fn" or die "\Unable to open output file: $tmpfiles_fn!";
while(<INPF>){
	my @line=split ',';
	if($line[0]=~/^\d+/ || $line[1]=~/\d+/){
		my $tmp_fn=gen_tempname(10);
		push @tmp_fn,$tmp_fn;
		$tmp_fn.=".csv";
		print OUTF "$tmp_fn\n";
		open OPF, ">$tmp_fn" or die "\nUnable to open output file: $tmp_fn!";
		print OPF $_;
		close OPF;
	}
}
close INPF;
close OUTF;

my $str="cat $tmpfiles_fn|parallel $epistat_cmd -x $xparr_fn -a figtree -s {} $epistat_prm";
print STDERR "\n\nPrinting trees!\n$str";
system($str);
print_child_termination_status();

my ($basename,$dir,$ext) = fileparse($site_pairs_fn,'\.[^.]*$');
my $outfn=$dir.$basename.".mutations.tre";
open OUTF, ">$outfn" or die "\Unable to open output file: $outfn!";
print OUTF "#NEXUS\nbegin trees;";
for(my $i=0;$i<@tmp_fn;$i++){
	my $tmp_fn=$tmp_fn[$i];
	open INPF, "<$tmp_fn.tre" or warn "\nUnable to open input file $tmp_fn.tre!";
	while(<INPF>){
		chomp;
		if(/^\ttree\s+(\w+)\s*=\s*/){
			my $tree_name=$1;
			my $tree_str=$';
			my $j=$i+1;
			$tree_name.="_".$j;
			print OUTF "\n\ttree $tree_name = $tree_str";
		}
	}
	close INPF;
	unlink glob $tmp_fn.".*";
}
print OUTF "\nend;";
close OUTF;
unlink $tmpfiles_fn;


