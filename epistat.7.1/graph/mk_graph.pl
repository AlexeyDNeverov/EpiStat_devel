#!/usr/bin/env perl
#This script creates two graph slices with positive and negative weights
#Options:
#	-p <pos_edge_list> - A positive graph slice edge list
#	-n <neg_edge_list> - A negative graph slice edge list
#	[-u] - Print indirected graph
#	-o <ouput_fn_prefix> - A prefix string for two output files (negative and positive slices)
#	-f <"Pajek"|"GraphML"|"Simap"> - An output graph format
#edge list file format:
#	bg_site\tfg_site\tweight

use strict;
use Getopt::Std;
use Class::Struct;

my %args;
if(!getopts('p:n:o:f:u',\%args)){
	die "\nError in option string!";
}
my $pos_graph_fname;
my $neg_graph_fname;
my $out_graph_fname;
my $f_directed=1; #1-directed 0-indirected
my $f_out_graph_format;

if($args{p}){
	$pos_graph_fname=$args{p};
}

if($args{n}){
	$neg_graph_fname=$args{n};
}

if($args{o}){
	$out_graph_fname=$args{o};
}else{
	die "\nThe prefix of output files wasn't specified!";
}

if($args{u}){
	$f_directed=0;
}

if($args{f}){
	if($args{f}=~m/pajek/i){
		$f_out_graph_format="pajek";
	}elsif($args{f}=~m/graphml/i){
		$f_out_graph_format="graphml";
	}elsif($args{f}=~m/simap/i){
		$f_out_graph_format="simap";
	}else{
		die "\nUnknown output graph format: ".$args{f};
	}
}else{
	die "\nThe output graph format wasn't specified!";
}

struct EdgeInfo =>{
	node1 => '$',
	node2 => '$',
	weight => '$'
};

my %pos_edge_weights;
my %neg_edge_weights;
my %vertices;
if(defined $pos_graph_fname){
	open INFILE, "<$pos_graph_fname" or die "\nUnable to open input file: $pos_graph_fname!";
	my $wcol_idx=-1;
	while(<INFILE>){
		chomp;
		my @line=split '\t';
		next unless ($line[0]=~m/^\d+$/)&&($line[1]=~m/^\d+$/);
		next unless ($line[$wcol_idx]>0)&&($line[0] != $line[1]);
		$vertices{$line[0]}=1;
		$vertices{$line[1]}=1;
		my $key=$line[0].",".$line[1];
		if($f_directed){
			$pos_edge_weights{$key}=$line[$wcol_idx];
		}else{
			my $key1=$line[1].",".$line[0];
			if(!(defined $pos_edge_weights{$key}||defined $pos_edge_weights{$key1})){
				$pos_edge_weights{$key}=$line[$wcol_idx];
			}
		}
	}
	close INFILE;
}else{
	die "\nFile with positive associatins is required! Use -p option"
}

if(defined $neg_graph_fname){
	open INFILE, "<$neg_graph_fname" or die "\nUnable to open input file:$neg_graph_fname!";
	my $wcol_idx=-1;
	while(<INFILE>){
		chomp;
		my @line=split '\t';
		next unless ($line[0]=~m/^\d+$/)&&($line[1]=~m/^\d+$/);
		next unless ($line[$wcol_idx]<0)&&($line[0] != $line[1])&&defined($vertices{$line[0]})&&defined($vertices{$line[1]});
		my $key=$line[0].",".$line[1];
		if($f_directed){
			$neg_edge_weights{$key}=$line[$wcol_idx];
		}else{
			my $key1=$line[1].",".$line[0];
			if(!(defined $neg_edge_weights{$key}||defined $neg_edge_weights{$key1})){
				$neg_edge_weights{$key}=$line[$wcol_idx];
			}
		}
	}
	close INFILE;
}

sub print_graphML{
	my ($fh,$gr_name,$rh_vertices,$rh_edges)=@_;
	print {$fh} "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"; 
	print {$fh} "\n<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\">";
	print {$fh} "\n<key id=\"weight\" for=\"edge\" attr.name=\"weight\" attr.type=\"double\"/>";
	print {$fh} "\n\t<graph id=\"$gr_name\" edgedefault=\"";
	if($f_directed){
		print {$fh} "directed";
	}else{
		print {$fh} "undirected";
	}
	print {$fh} "\">";
	foreach my $vname(keys %{$rh_vertices}){
		print {$fh} "\n\t\t<node id=\"$vname\"/>";
	}
	foreach my $str(keys %{$rh_edges}){
		my @line=split ',', $str;
		my ($site1,$site2)=($line[0],$line[1]);
		if(!$f_directed){
			if($site1>$site2){
				$site1=$site2;
				$site2=$line[0];
			}
		}
		my $w=$rh_edges->{$str};
		print {$fh} "\n\t\t<edge source=\"$site1\" target=\"$site2\">";
		print {$fh} "\n\t\t\t<data key=\"weight\">$w</data>";
		print {$fh} "\n\t\t</edge>";
	}
	print {$fh} "\n\t</graph>";
	print {$fh} "\n</graphml>";
}

sub print_pajek{
	my $fh=shift;
	my $rh_vertices=shift;
	my @edge_slices=@_;
	my %site2id;
	my $n=keys %{$rh_vertices};
	print {$fh} "*Vertices $n";
	$n=1;
	foreach my $site(sort {$a<=>$b} keys %{$rh_vertices}){
		print {$fh} "\n$n \"$site\"";
		$site2id{$site}=$n++;
	}
	if($f_directed){
		print {$fh} "\n*Arcs";
	}else{
		print {$fh} "\n*Edges";
	}
	my @edges;
	foreach my $rh_edges(@edge_slices){
		foreach my $str(keys %{$rh_edges}){
			my @line=split ',', $str;
			my $w=$rh_edges->{$str};
			my ($site1,$site2)=($site2id{$line[0]},$site2id{$line[1]});
			my $ei=EdgeInfo->new();
			if((!$f_directed)&&$site1>$site2){
				$ei->node1($site2);
				$ei->node2($site1);
			}else{
				$ei->node1($site1);
				$ei->node2($site2);
			}
			$ei->weight($w);
			push @edges,$ei;
		}
	}
	@edges=sort {$a->node1<=>$b->node1||$a->node2<=>$b->node2} @edges;
	foreach my $ei(@edges){
		my ($site1,$site2)=($ei->node1,$ei->node2);
		print {$fh} "\n".$site1." ".$site2." ".$ei->weight;
	}
}

sub print_simap{
	my $fh=shift;
	my $rh_vertices=shift;
	my @edge_slices=@_;
	my @edges;
	foreach my $rh_edges(@edge_slices){
		foreach my $str(keys %{$rh_edges}){
			my @line=split ',', $str;
			my $w=$rh_edges->{$str};
			my ($site1,$site2)=($line[0],$line[1]);
			my $ei=EdgeInfo->new();
			if((!$f_directed)&&$site1>$site2){
				$ei->node1($site2);
				$ei->node2($site1);
			}else{
				$ei->node1($site1);
				$ei->node2($site2);
			}
			$ei->weight($w);
			push @edges,$ei;
		}
	}
	@edges=sort {$a->node1<=>$b->node1||$a->node2<=>$b->node2} @edges;
	foreach my $ei(@edges){
		my ($site1,$site2)=($ei->node1,$ei->node2);
		print {$fh} $site1." ".$site2." ".$ei->weight."\n";
	}
}

if($f_out_graph_format eq "graphml"){
	my $ofname=$out_graph_fname.".positive_edges.graphml";
	open OPF, ">$ofname" or die "\nUnable to open output file: $ofname!";
	print_graphML(*OPF,"Gpos",\%vertices,\%pos_edge_weights);
	close OPF;
	$ofname=$out_graph_fname.".negative_edges.graphml";
	open OPF, ">$ofname" or die "\nUnable to open output file: $ofname!";
	print_graphML(*OPF,"Gneg",\%vertices,\%neg_edge_weights);
	close OPF;
}elsif($f_out_graph_format eq "pajek"){
	my $ofname=$out_graph_fname.".net";
	open OPF, ">$ofname" or die "\nUnable to open output file: $ofname!";
	print_pajek(*OPF,\%vertices,\%pos_edge_weights,\%neg_edge_weights);
	close OPF;
}elsif($f_out_graph_format eq "simap"){
	my $ofname=$out_graph_fname.".positive_edges.simap";
	open OPF, ">$ofname" or die "\nUnable to open output file: $ofname!";
	print_simap(*OPF,\%vertices,\%pos_edge_weights);
	close OPF;
	$ofname=$out_graph_fname.".all_edges.simap";
	open OPF, ">$ofname" or die "\nUnable to open output file: $ofname!";
	print_simap(*OPF,\%vertices,\%pos_edge_weights,\%neg_edge_weights);
	close OPF;
}

