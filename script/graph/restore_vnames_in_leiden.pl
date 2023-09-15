#!/usr/bin/env perl
#converts vertex indices into vertex names in the leiden output
#Usage: <index2vertex_name_fn> <leiden_communities_fn>
open INPF, $ARGV[0] or die "\nUnable to open input file $ARGV[0]!";
my %idx2vname;
my $i=0;
while(<INPF>){
	$idx2vname{$i}=$1 if(/^(\d+)$/);
	$i++;
}
close INPF;
open INPF, $ARGV[1] or die "\nUnable to open input file $ARGV[1]!";
my %vname2group;
print "vertex\tgroup\n";
while(<INPF>){
	my @line=split '\t',$_,-1;
	my $vname=$idx2vname{$line[0]};
	if(defined $vname){
		print "$vname\t$line[1]";
	}else{
		die "\nThe file with vertex names is inconsistent with the leiden output file!";
	}
}
close INPF;