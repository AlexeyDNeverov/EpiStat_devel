#!/usr/bin/env perl
#The script calculates distances on a protein structure between sites in intragenic site pairs
#usage: options <site_pairs_file>
#options
#	-p <FN> - A file with protein structures in PDB format
#	-c <CHAINS> - A string of chain labels in PDB delimited by ','
#	[-o] <FN> - A name of output file
#	-B <FN> - A name of file with a table for mapping of sites of marker sequence on positions of a primary sequence retrieved from ATOM records in the PDB file
#	-A <FN> - A name of file with the marker sequence taken from the initial alignment file as is (with gaps)

use strict;
use Getopt::Std;
use IO::File;

my %args;
if(!getopts('p:c:o:B:A:',\%args)){
	die "\nError in option string!";
}

sub print_site_pair_distances {
	my ($pdb_file,$chains_str,$site2pdb_file,$align_mseq_file,$pairs_file,$out_file)=@_;
	my @chain_ids;
	if(defined $chains_str){
		@chain_ids=split ',',$chains_str;
	}else{
		die "\nThe chains corresponded to site pairs are not specified!";
	}
	my %atoms_CA;
	open PDB, "<$pdb_file" or die "\nopen <$pdb_file $!";
	while(<PDB>){
 		if (/^ATOM\s+[0-9]+\s+CA\s/){
 			my @splitter = split /\s+/;
			$atoms_CA{$splitter[4]}=[] unless defined $atoms_CA{$splitter[4]};
			push @{$atoms_CA{$splitter[4]}},[$splitter[6], $splitter[7], $splitter[8]];
 		}
	}
	close PDB;
	foreach my $chid(@chain_ids){
		die "\nThe chain $chid not detected in the input PDB: $pdb_file!" unless defined($atoms_CA{$chid});
	}
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
			$site2pdb{$site2align_pos[$1-1]}=$2-1 if ($1>0)&&($2>0);
			unless($1<=@site2align_pos){
				my $str=@site2align_pos;
				die "\nThe marker sequence Nsites=$str is inconsistent to the sequence used for mapping on the PDB (pos=$1)!";
			}
		}
	}
	close INPF;
	open PAIRS,  "<$pairs_file" or die "open <$pairs_file $!";
	my $out_fh;
	if(defined $out_file){
		$out_fh=IO::File->new(">$out_file") or die "open >$out_file $!";
	}else{
		$out_fh=*STDOUT;
	}
 	while(<PAIRS>){
 		chomp;
		s/^\s+//;
		if(/^(\d+)\s+(\d+)/){
			my $bg = $site2pdb{$1};
			my $fg = $site2pdb{$2};
			my @same_dist;
			my @inter_dist;
			my $samemin;
			my $intermin;
			if(defined($bg)&&defined($fg)){
				foreach my $chid(@chain_ids){
					my $dist=distance($atoms_CA{$chid}->[$bg], $atoms_CA{$chid}->[$fg]) if $bg<@{$atoms_CA{$chid}}&&$fg<@{$atoms_CA{$chid}};
					push @same_dist,$dist;
				}
				for(my $i=0;$i<@chain_ids;$i++){
					for(my $j=0;$j<@chain_ids;$j++){
						next if $j==$i;
						my $dist=distance($atoms_CA{$chain_ids[$i]}->[$bg], $atoms_CA{$chain_ids[$j]}->[$fg]) if $bg<@{$atoms_CA{$chain_ids[$i]}}&&$fg<@{$atoms_CA{$chain_ids[$j]}};
						push @inter_dist,$dist;
					}
				}
				
				if(@same_dist){
					$samemin=$same_dist[0];
					for(my $i=1;$i<@same_dist;$i++){
						$samemin=$same_dist[$i] if $same_dist[$i]<$samemin;
					}
				}
				
				if(@inter_dist){
					$intermin=$inter_dist[0];
					for(my $i=1;$i<@inter_dist;$i++){
						$intermin=$inter_dist[$i] if $inter_dist[$i]<$intermin;
					}
				}
			}
			my $min;
			if(defined($samemin)){
				if(!defined($intermin)||($samemin < $intermin)){
					print $out_fh "same\t";
					$min=$samemin;
				}
			}
			if(defined($intermin)){ 
				if(!defined($samemin)||($intermin < $samemin)){
					print $out_fh "inter\t";
					$min=$intermin;
				}
			}
			if(defined $min){
				print $out_fh $min;
			}else{
				print $out_fh "\t";
			}
			print $out_fh "\n";
		}
 	}
 	close PAIRS;
 	close $out_fh;
}

 sub distance{
 	my @coord1 = @{$_[0]};
 	my @coord2 = @{$_[1]};
 	my $dist = sqrt(($coord1[0] - $coord2[0])**2 + ($coord1[1] - $coord2[1])**2 + ($coord1[2] - $coord2[2])**2);
 	return $dist;
 }
 
print_site_pair_distances($args{p},$args{c},$args{B},$args{A},$ARGV[0],$args{o});