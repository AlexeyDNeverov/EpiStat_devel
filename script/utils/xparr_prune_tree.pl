#!/usr/bin/env perl
#This script removes specified taxa or clades from the xparr file
#Usage: options <xparr>
#options:
#	-l <FN> - List of nodes to prune
use strict;
use Bio::Phylo::IO;
use Getopt::Std;
my %args;
if(!getopts('l:',\%args)){
	die "\nError in option string!";
}
my $node_list=$args{l};
my %prune_nodes;
die "\nThe file with list of tree nodes to prune is required!" unless defined $node_list;
open INPF, "<$node_list" or die "\nUnable to open node list file: $node_list!";
while(<INPF>){
	chomp;
	if(/\S+/){
		$_=~s/^\s+//;
		$_=~s/\s+$//;
		$prune_nodes{$_}=undef;
	}
}
close INPF;

open INPF, "<$ARGV[0]" or die "\nUnable to open the tree file in XPARR format: $ARGV[0]!";
my $header;
my %xparr_lines;
while(<INPF>){
	if(/^(\S+)\s*[\t\n]/){
		chomp;
		unless(defined $header){
			$header=$_;
		}else{
			my $str=$1;
			$xparr_lines{$str}=[split("\t", $_,-1)];
		}
	}
}
close INPF;
my $tree=Bio::Phylo::IO->parse(
    -file => $ARGV[0],
    -format => 'adjacency');
#!!!The return value of Adjacency parser not corresponded to interface of the Bio::Phylo::IO->parse
$tree=$tree->[0]->first;
foreach my $node($tree->get_nodes){
	my $name=$node->get_name();
	if(exists $prune_nodes{$name}){
		$prune_nodes{$name}=$node;
	}
}
my @tips;
my $str;
foreach my $name(keys %prune_nodes){
	if(defined $prune_nodes{$name}){
		my $node=$prune_nodes{$name};
		if($node->is_terminal){
			push @tips,$name;
		}else{
			push @tips,@{$node->get_terminals()};
		}
	}else{
		$str.=";" if defined $str;
		$str.=$name;
	}
}
{
	my @tmp;
	@tips=sort @tips;
	push @tmp,$tips[0];
	for(my $i=1;$i<@tips;$i++){
		push @tmp,$tips[$i] if($tmp[-1] ne $tips[$i]);
	}
	@tips=@tmp;
}
unless(defined $str){
	$tree->prune_tips(\@tips);
	print $header;
	$tree->visit_breadth_first(
		-in => sub {
			my $node=shift;
			my $name=$node->get_name;
			if(defined $xparr_lines{$name}){
				my @tmp=@{$xparr_lines{$name}};
				my @sites=({},{},{});
				for(my $i=3;$i<@tmp;$i++){
					my @muts=split ";",$tmp[$i];
					foreach my $mut(@muts){
						if($mut=~/([A-Z]?)(\d+)([A-Z]?)/){
							$sites[$i-3]->{$2}=[($1,$2,$3)];
						}
					}
				}
				unless($node->is_root){
					my $pnode=$node->get_parent;
					my $pname=$pnode->get_name;
					my $pname1=$tmp[1];
					my $n=0;
					while($pname ne $pname1){
						#join mutation lists on consecutive branches
						if(defined $xparr_lines{$pname1}){
							my @tmp1=@{$xparr_lines{$pname1}};
							for(my $i=3;$i<@tmp1;$i++){
								my @muts=split ";",$tmp1[$i];
								foreach my $mut(@muts){
									if($mut=~/([A-Z]?)(\d+)([A-Z]?)/){
										my ($aa,$da)=($1,$3);
										my $pos=$2;
										if(defined $sites[$i-3]->{$pos}){
print STDERR "\n$name\t$pname1\t$pname\t".join "", @{$sites[$i-3]->{$pos}};
print STDERR "\n$pname1\t$tmp1[1]\t$mut";
											if(($aa ne "")&&($da ne "")){
												if($da eq $sites[$i-3]->{$pos}->[0]){
													if($aa ne $sites[$i-3]->{$pos}->[2]){
														$sites[$i-3]->{$pos}->[0]=$aa;
													}else{
														delete $sites[$i-3]->{$pos};
													}
												}else{
													die "\nInconsistence in the site $pos in consecutive branches: $name, $pname1!";
												}
											}
										}else{
											$sites[$i-3]->{$pos}=[($aa,$pos,$da)];
										}
									}
								}
							}
							$pname1=$tmp1[1];
							$n++;
						}else{
							die "\nUnable to find a node name: $pname1 in the input XPARR file!";
						}
					}
					if($n){
						for(my $i=0;$i<3;$i++){
							my $str;
							foreach my $pos (sort {$a<=>$b} keys %{$sites[$i]}){
								$str.=";" if defined $str;
								$str.=join "",@{$sites[$i]->{$pos}};
							}
							$tmp[3+$i]=$str;
						}
						$tmp[1]=$pname;
						$tmp[2]=$node->get_branch_length;
					}
				}
				print "\n".join "\t",@tmp;
			}else{
				die "\nUnable to find a node name: $name in the input XPARR file!";
			}
		}
	)
}else{
	die "\nSome node names requested to prune are absent in the input XPARR file:\n\t$str!"
}