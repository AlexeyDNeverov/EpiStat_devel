#!/usr/bin/env perl
#This script calculates epistatic statistics [Kryazhimsky11]
#options:
#	-t  - Ignore terminal nodes
#	-b <0<Q<=1> - Split tree on branches with relative lengths upper than Q percentile
#		Default=1
#	[-u (DEFAULT|SYN|NSYN|TOTAL)] - Branch length units
#		default=DEFAULT
#Params:
# <FN> - input file for EpiStat application in XPAR format [Kryazhimsky11]. Contains tree in the adjacent file format

use strict;
use Bio::Phylo::IO;
use Getopt::Std;

my %args;
if(!getopts('tu:b:',\%args)){
	die "\nError in option string!";
}
		
open INFILE, "$ARGV[0]" or die "\nUnable to open input file: $ARGV[0]!";
my $xpar_fn=$ARGV[0];
my $f_ignore_terminals=0;
$f_ignore_terminals=1 if defined $args{t};
#0 - default
#1 - number of synonimous substitutions
#2 - number of non synonimous substitutions
#3 - total number of substitutions
my $distance_mode=0;
if(defined $args{u}){
	if($args{u}=~m/^default/i){
		$distance_mode=0;
	}elsif($args{u}=~m/^syn/i){
		$distance_mode=1;
	}elsif($args{u}=~m/^nsyn/i){
		$distance_mode=2;
	}elsif($args{u}=~m/^total/i){
		$distance_mode=3;
	}else{
		die "\nUndefined distance mode -m $args{u}";
	}
}
my $Q=1;
$Q=$args{b} if defined $args{b};
die "\nWrong value for -b option (0,1] is expected!" unless $Q>0&&$Q<=1;

my $tree=Bio::Phylo::IO->parse(
    -file => $xpar_fn,
    -format => 'adjacency');
#!!!The return value of Adjacency parser not corresponded to interface of the Bio::Phylo::IO->parse
$tree=$tree->[0]->first;

#my $str=Bio::Phylo::IO->unparse(-phylo => $tree,-format => 'newick');
#print $str;

open INPF, "<$xpar_fn" or die "\nUnable to open input file $xpar_fn!";
$_=<INPF>; #skip header line
die "\nHeader line was omitted in file: $xpar_fn!" unless /^child\tparent\tlength/;
$_=<INPF>; #skip root node statement
die "\nThe root node statement required in file: $xpar_fn!" unless /^[\w.]+\s*$/;
my %syn_counts;
my %nsyn_counts;
while(<INPF>){
	chomp;
	my @line=split "\t";
	my $n=0;
	if($line[4] ne ""){
		while($line[4]=~m/;/g){$n++};
		$n++;
	};
	$nsyn_counts{$line[0]}=$n;
	if($distance_mode){
		$n=0;
		if($line[3] ne ""){
			while($line[3]=~m/;/g){$n++};
			$n++;
		};
		$syn_counts{$line[0]}=$n;
	}
};
close INPF;
if($distance_mode){
	foreach my $node($tree->get_nodes){
		if($distance_mode==1){
			$node->set_branch_length($syn_counts{$node->get_name});
		}elsif($distance_mode==2){
			$node->set_branch_length($nsyn_counts{$node->get_name});
		}elsif($distance_mode==3){
			$node->set_branch_length($syn_counts{$node->get_name}+$nsyn_counts{$node->get_name});
		}
	}
}
my $max_branch_length;
my $median_branch_length;
{
	my @brlengths;
	foreach my $node($tree->get_nodes){
		push @brlengths,$node->get_branch_length;
	}
	@brlengths=sort {$a<=>$b} @brlengths;
	my $i=int($Q*@brlengths)-1;
	$max_branch_length=$brlengths[$i];
	$i=int($Q*@brlengths/2);
	$median_branch_length=$brlengths[$i];
	$median_branch_length=($median_branch_length+$brlengths[$i-1])/2 unless int($Q*@brlengths)%2;
}
my ($n1,$n2)=(0,0);
$tree->visit_breadth_first(
	-in => sub{
		my $node=shift;
		if(!$node->is_root){
			my $brlength=$node->get_branch_length;
			if($brlength>$max_branch_length){
				$node->set_generic('bref'=>undef)
			}else{
				my $pbref=$node->get_parent->get_generic('bref');
				if($pbref){
					$node->set_generic('bref'=>$pbref);
				}else{
					$node->set_generic('bref'=>$node);
				}
			}
			if($n1<$nsyn_counts{$node->get_name}){
				$n1=$nsyn_counts{$node->get_name};
			}elsif($n2<$nsyn_counts{$node->get_name}){
				$n2=$nsyn_counts{$node->get_name};
			}
			my $len=$node->get_parent()->get_generic('time')+$node->get_branch_length();
			$node->set_generic('time'=>$len);
		}else{
			$node->set_generic('time'=>0);
			$node->set_generic('bref'=>undef);
		}
	}
);


#begin script
my $scale=$n1*$n2;
my $norm=0;
my $tau=0;

sub get_distance{
	my ($node1,$node2)=@_;
	return -1 if $node1->is_root();
	return -1 unless $node1->get_generic('bref')==$node2->get_generic('bref');
	my $dist=$node2->get_generic('time')-$node1->get_generic('time');
	return -1 if($dist<0);
	if($node1==$node2){
		$dist=$node2->get_branch_length()/3;
	}else{
		$dist-=$node2->get_branch_length()/2;
		$dist+=$node1->get_branch_length()/2;
	}
	return $dist;
}

#calculate site pairs
my $nbranches=$tree->calc_number_of_nodes();
if($f_ignore_terminals){
	$nbranches-=$tree->calc_number_of_terminals();
}
print STDERR "\nStart calculations (NBranches=$nbranches):";
$tree->visit_breadth_first(
	-in   => sub {
		my $node=shift;
		if(!($node->is_root||($f_ignore_terminals&&$node->is_terminal)||(!defined($node->get_generic('bref'))))){
			my $ptr=$node;
			while(!$ptr->is_root){
				my $l=get_distance($ptr,$node);
				last if $l<0;
				my $n=$nsyn_counts{$node->get_name};
				if($ptr!=$node){
					$n*=$nsyn_counts{$ptr->get_name};
				}else{
					#$n*=($n-1)/2;
					$n*=$n;
				}
				$tau+=$l*$n/$scale;
				$norm+=$n/$scale;
				$ptr=$ptr->get_parent;
			}
		}	
	}
);
$tau/=$norm;
print "XPARR=\"$xpar_fn\"";
print "\nDistances=";
if($distance_mode==0){
	print "\"default\"";
}elsif($distance_mode==1){
	print "\"synonymous\"";
}elsif($distance_mode==2){
	print "\"nonsynonymous\"";
}elsif($distance_mode==3){
	print "\"total\"";
}
if($Q<1){
	print "\nBranch filter Q=\"$Q\"";
}
print "\nMax branch length=\"$max_branch_length\"";
print "\nMedian branch length=\"$median_branch_length\"";
print "\nTAU=\"$tau\"";