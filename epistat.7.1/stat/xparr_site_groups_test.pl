#!/usr/bin/env perl
#This script identifies branches on the phylogenetic tree those subtrees have unusual distribution of mutation counts in initially defined subgroups of sites.
#options:
#	-x <FN> - input file for EpiStat application in XPAR format [Kryazhimsky11]. Contains tree in the adjacent file format
#	-a <mult[inomial]|crosstab|no> - test applied for branches. 
#		If 'no' test applied the script displays mutation counts on the tree nodes and branches
#	[-b] <FN> - input file in TAB format with sites in background to groups assignment
#	[-t] <FN> - input file in TAB format with sites in foreground to groups assignment
#		Note: '-b' and '-t' options are mutual exclusive. At least one is required!
#	[-c] <FN> - input file with groups' colour codes 
use strict;
use lib "$ENV{EPISTAT_LIB}";
my $multinom_r="$ENV{EPISTAT_HOME}/stat/multinomial_test.R";
my $crosstab_r="$ENV{EPISTAT_HOME}/stat/fishers_exact_test.R";

use Bio::Phylo::IO;
use Class::Struct;
use Getopt::Std;
use Time::Progress;
use File::Basename;
use TreeUtils::Phylo::FigTree qw(tree2str);
my %args;
if(!getopts('x:b:t:c:a:',\%args)){
	die "\nError in option string!";
}

my $alpha=0.05;
my %site2group;
my $ngroups=0;
my $xparr_fn=$args{x};
my $colors_fn=$args{c};
my $groups_fn;
my $branch_test;
my @branch_tests=("no","multinomial","crosstab");
if(defined $args{a}){
	if($args{a}=~/^mult(inomial)*/){
		$branch_test=1;
	}elsif($args{a}=~/^crosstab/){
		$branch_test=2;
	}elsif($args{a}=~/^no/){
		$branch_test=0;
	}else{
		die "\nUnproper value $args{a} for parameter '-a'!"; 
	}
}else{
	die "\nThe type of test for branches isn't accounted!";
}
if(defined $args{b}){
	unless(defined($args{t})){
		$groups_fn=$args{b};
	}else{
		die "\nOnly one option: '-b' or '-t' is accepted!";
	}
}
if(defined $args{t}){
	unless(defined($args{b})){
		$groups_fn=$args{t};
	}else{
		die "\nOnly one option: '-b' or '-t' is accepted!";
	}
}

die "\nThe input XPAR file is not specified! Use -x option!" unless defined $xparr_fn;
die "\nThe input TAB file with sites' ID to groups' ID table is not specified! Use -b or -t options!" unless defined $groups_fn;

my $tree=Bio::Phylo::IO->parse(
    -file => $xparr_fn,
    -format => 'adjacency');
#!!!The return value of Adjacency parser not corresponded to interface of the Bio::Phylo::IO->parse
$tree=$tree->[0]->first;
my $ntests=@{$tree->get_internals};
#substitutions
my %subst_map;
my %groups_counts;
	
struct SubstInfo => {
	site => '$',
	weight => '$',
	bases => '@'
};

sub parse_subst_abbr{
	my $str=shift;
	my $ra_bases=shift;
	my $spos;
	if($str=~/([A-Za-z])(\d+)([[A-Za-z])/){
		@{$ra_bases}=($1,$3) if(defined $ra_bases);
		$spos=$2;
	}else{
		$str=~s/\D//g;
		$spos=$str;
	};
	return $spos;
}

sub init_subst_map{
	my ($xpar_fn,$tree,$rh_bg_subst_map,$rh_tg_subst_map,$rh_syn_counts,$rh_nsyn_counts)=@_;
	die "\nAt least one hashref must be defined!" unless defined($rh_bg_subst_map)||defined($rh_tg_subst_map);
	open INPF, "<$xpar_fn" or die "\nUnable to open input file $xpar_fn!";
	$_=<INPF>; #skip header line
	die "\nHeader line was omitted in file: $xpar_fn!" unless /^child\tparent\tlength/;
	$_=<INPF>; #skip root node statement
	die "\nThe root node statement required in file: $xpar_fn!" unless /^\w+\s*$/;
	my %nsyn_counts;
	my $I=0;
	my $J=0;
	while(<INPF>){
		chomp;
		my @line=split "\t";
		#next if $f_ignore_terminals&&$terminals{$line[0]};
		if(defined($rh_tg_subst_map)){
			my $n=0;
			if($line[4] ne ""){
				my @sites=split(/;/, $line[4]);
				$n=@sites;
				$rh_tg_subst_map->{$line[0]}={};
				foreach my $str(@sites){
					my @bases;
					my $spos=parse_subst_abbr($str,\@bases);
					my $si=SubstInfo->new();
					$si->site($spos);
					$si->bases([@bases]) if @bases;
					$rh_tg_subst_map->{$line[0]}->{$spos}=$si;
					$n++;
				}
			}
			$rh_nsyn_counts->{$line[0]}=$n if defined $rh_nsyn_counts;
		}
		if(defined($rh_bg_subst_map)){
			my $n=0;
			if($line[5] ne ""){
				my @sites=split(/;/, $line[5]);
				$rh_bg_subst_map->{$line[0]}={};
				foreach my $str(@sites){
					my @bases;
					my $spos=parse_subst_abbr($str,\@bases);
					my $si=SubstInfo->new();
					$si->site($spos);
					$si->bases([@bases]) if @bases;
					$rh_bg_subst_map->{$line[0]}->{$spos}=$si;
					$n++;
				}
			}
			$rh_nsyn_counts->{$line[0]}=$n if defined $rh_nsyn_counts;
		}
		if(defined $rh_syn_counts){
			my $n=0;
			chomp $line[3];
			if($line[3] ne ""){
				while($line[3]=~m/;/g){$n++};
				$n++;
			}
			$rh_syn_counts->{$line[0]}=$n;
		}
	}
	close INPF;
}
#possible error here:
#   counts on branches include mutations on the upward edges!!!
#fixed
sub calc_groups_counts{
	my ($tree,$ngroups,$rh_site2group,$rh_subst_map,$rh_groups_counts)=@_;
	%{$rh_groups_counts}=();
	$tree->visit_depth_first(
		-pre => sub {
				my $node=shift;
				my $name=$node->get_name;
				die "\nUnnamed nodes in the tree are not allowed!" unless defined $name;
				$rh_groups_counts->{$name}={"branch" => [(0) x $ngroups], "node" => [(0) x $ngroups]};
			},
		-in => sub{
				my $node=shift;
				my $name=$node->get_name;
				foreach my $chnode(@{$node->get_children}){
					my $chname=$chnode->get_name;
					for(my $i=0;$i<$ngroups;$i++){
						$rh_groups_counts->{$name}->{node}->[$i]+=$rh_groups_counts->{$chname}->{node}->[$i];
					}
					if(defined $rh_subst_map->{$chname}){
						foreach my $site(keys %{$rh_subst_map->{$chname}}){
							my $grid=$rh_site2group->{$site};
							if(defined $grid){
								$rh_groups_counts->{$name}->{node}->[$grid]++;
								$rh_groups_counts->{$chname}->{branch}->[$grid]++;
							}
						}
					}
				}
			}
	);
}

sub gen_tempname{
	my $nchar=shift;
	my @chars = ( "A" .. "Z", "a" .. "z", 0 .. 9 );
	return join("", @chars[ map { rand @chars } ( 1 .. $nchar ) ]);
}

sub print_child_termination_status{
	my ($node_name)=@_;
	if ($? == -1) {
		print "$node_name: failed to execute: $!\n";
		exit 1;
	}elsif ($? & 127) {
		printf "\t$node_name: child died with signal %d, %s coredump\n",
					($? & 127),  ($? & 128) ? 'with' : 'without';
	}else {
		if($? >> 8!=0){
			printf "\t$node_name: child exited with value %d\n", $? >> 8;
		};
	}
}

sub multinomial_test{
	my ($node,$ref_name,$rh_group_counts)=@_;
	my $node_name=$node->get_name;
	my $ra_ref_counts=$rh_group_counts->{$ref_name}->{node};
	my $ra_counts=$rh_group_counts->{$node_name}->{node};
	my $tmp_fname=$node_name;
	$tmp_fname.="." if defined $tmp_fname;
	$tmp_fname.=gen_tempname(10);
	my $n=0;
	my $ng=0;
	for(my $i=0;$i<$ngroups;$i++){
		$ng++ if $ra_ref_counts->[$i]>0;
		$n+=$ra_counts->[$i];
	}
	return undef unless $n>10*$ng;
	return undef unless $ng>1;
	$n=0;
	for(my $i=0;$i<$ngroups;$i++){
		$n+=$ra_ref_counts->[$i];
	}
	return undef unless $n;
	my $pvalue;
	open OUTF, ">$tmp_fname" or die "\nUnable to open output file: $tmp_fname!";
	print OUTF "P\tX";
	for(my $i=0;$i<$ngroups;$i++){
		my $str=$ra_ref_counts->[$i]/$n;
		next unless $str>0;
		#$str=sprintf("%.4f",$str);
		print OUTF "\n".$str."\t".$ra_counts->[$i];
	}
	print OUTF "\n";
	close OUTF;
	my $tmp_out=$tmp_fname.".multinom_test.R.out";
	my $str=$multinom_r." $tmp_fname >$tmp_out";
	system($str);
	print_child_termination_status($node_name);
	open INPF, "<$tmp_out" or die "\nUnable to open input file: $tmp_out!";
	while(<INPF>){
		chomp;
		if(/\S+/){
			s/^\s+//;
			s/\s+$//;
			if(/P value.+?=\s*(\S+)/){
				$pvalue=$1;
				last;
			}
		}
	}
	close INPF;
	unless($? >> 8!=0){
		unlink $tmp_out;
		unlink $tmp_fname;
	}
	return $pvalue;
}

sub crosstab_test{
	my ($node,$ref_name,$rh_group_counts)=@_;
	my $node_name=$node->get_name;
	my @ref_counts=@{$rh_group_counts->{$ref_name}->{node}};
	my @counts=@{$rh_group_counts->{$node_name}->{node}};
	my $n=0;
	my $ng=0;
	for(my $i=0;$i<$ngroups;$i++){
		$counts[$i]+=$rh_group_counts->{$node_name}->{branch}->[$i];
		$ref_counts[$i]-=$counts[$i];
	}
	for(my $i=0;$i<$ngroups;$i++){
		if($ref_counts[$i]+$counts[$i]>0){
			$ng++;
			$n+=$counts[$i];
		}
	}
	return undef unless $ng>1;
	return undef unless $n>=$ng;
	$n=0;
	for(my $i=0;$i<$ngroups;$i++){
		$n+=$ref_counts[$i];
	}
	return undef unless $n>=$ng;
	my $pvalue;
	my $tmp_fname=$node_name;
	$tmp_fname.="." if defined $tmp_fname;
	$tmp_fname.=gen_tempname(10);
	open OUTF, ">$tmp_fname" or die "\nUnable to open output file: $tmp_fname!";
	my @tmp_counts=([],[]);
	for(my $i=0;$i<$ngroups;$i++){
		next unless $ref_counts[$i]+$counts[$i]>0;
		push @{$tmp_counts[0]},$ref_counts[$i];
		push @{$tmp_counts[1]},$counts[$i];
	}
	my $str=join "\t", @{$tmp_counts[0]};
	print OUTF $str;
	$str=join "\t", @{$tmp_counts[1]};
	print OUTF "\n$str\n";
	close OUTF;
	my $tmp_out=$tmp_fname.".crosstab_test.R.out";
	my $str=$crosstab_r." $tmp_fname >$tmp_out";
	system($str);
	print_child_termination_status($node_name);
	open INPF, "<$tmp_out" or die "\nUnable to open input file: $tmp_out!";
	while(<INPF>){
		chomp;
		if(/\S+/){
			s/^\s+//;
			s/\s+$//;
			if(/p-value.+?[=<]\s*(\S+)/){
				$pvalue=$1;
				last;
			}
		}
	}
	close INPF;
	unless($? >> 8!=0){
		unlink $tmp_out;
		unlink $tmp_fname;
	}
	return $pvalue;
}

sub calc_KL_distance{
	my ($ra_model_freqs,$ra_obs_freqs,$vsize)=@_;
	my $d;
	my $nf=0;
	my $no=0;
	my $epsilon=1e-3;
	for(my $i=0;$i<$vsize;$i++){
		if($ra_model_freqs->[$i]>0){
			if($ra_obs_freqs->[$i]>0){
				$d+=$ra_obs_freqs->[$i]*log($ra_obs_freqs->[$i]/$ra_model_freqs->[$i]);
			}
			$nf+=$ra_model_freqs->[$i];
			$no+=$ra_obs_freqs->[$i];
		}else{
			$d=undef;
			last;
		}
	}
	if(defined $d){
		die "\ncalc_KL_distance() Error: Value of model frequencies are not normalized: $nf!" if abs($nf-1.0)>$epsilon;
		die "\ncalc_KL_distance() Error: Value of observed frequencies are not normalized: $no!" if abs($no-1.0)>$epsilon;
	}
	return $d;
}

sub calc_KLD_for_prop_models{
	my ($ra_model_freqs,$ra_obs_freqs,$vsize)=@_;
	my @kld=(0) x ($vsize+1);
	$kld[0]=calc_KL_distance($ra_model_freqs,$ra_obs_freqs,$vsize);
	if($vsize>2){
		my @prop_models;
		for(my $i=0;$i<$vsize;$i++){
			push @prop_models, [(0) x $vsize];
			if($ra_obs_freqs->[$i]>0&&$ra_obs_freqs->[$i]<1){
				$prop_models[$i]->[$i]=$ra_obs_freqs->[$i];
				for(my $j=0;$j<$vsize;$j++){
					$prop_models[$i]->[$j]=(1-$ra_obs_freqs->[$i])*$ra_model_freqs->[$j]/(1-$ra_model_freqs->[$i]) if $i!=$j;
				}
			}else{
				$kld[$i+1]=undef;
			}
		}
		for(my $i=0;$i<$vsize;$i++){
			if(defined $kld[$i+1]){
				$kld[$i+1]=calc_KL_distance($prop_models[$i],$ra_obs_freqs,$vsize);
			}
		}
	}
	return @kld;
}
#Start script
open INPF, "<$groups_fn" or die "\nUnable to open input file: $groups_fn!";
my %groups;
while(<INPF>){
	chomp;
	s/^\s+//;
	s/\s+$//;
	my @line=split '\s+';
	next unless ($line[0]=~/^\d+$/)&&($line[1]=~/^\d+$/);
	if(($line[0]>0)&&($line[1]>0)){
		$site2group{$line[0]}=$line[1]-1;
		$groups{$line[1]}=1;
	}else{
		die "\nError in the input file $groups_fn: the site and group numbering must start from one!";
	}
}
close INPF;
$ngroups=keys %groups;
if($args{b}){
	init_subst_map($xparr_fn,$tree,\%subst_map,undef);
}elsif($args{t}){
	init_subst_map($xparr_fn,$tree,undef,\%subst_map);
}
calc_groups_counts($tree,$ngroups,\%site2group,\%subst_map,\%groups_counts);
my $tmp_file_list=gen_tempname(10);
my %nodes_tests;
struct NodeTestInfo => {
	ref_node_name => '$',
	test_pvalue => '$',
	test_ref_node_name => '$'
};
my @sign_tests;
if($branch_test){
	print STDERR "\nStart tests on the tree branches:";
	my $p = Time::Progress->new(min => 1, max => $ntests);
	my $nti=NodeTestInfo->new();
	my $name=$tree->get_root->get_name;
	$nti->ref_node_name($name);
	$nodes_tests{$name}=$nti;
	my $i=1;
	$tree->visit_breadth_first(
		-in => sub{
				my $node=shift;
				my $name=$node->get_name;
				unless($node->is_root||$node->is_terminal){
					my $pnode=$node->get_parent;
					my $pname=$pnode->get_name;
					my $ref_node_name=$nodes_tests{$pname}->ref_node_name();
					my $test_pvalue;
					if($branch_test==1){
						$test_pvalue=multinomial_test($node,$ref_node_name,\%groups_counts);
					}elsif($branch_test==2){
						$test_pvalue=crosstab_test($node,$ref_node_name,\%groups_counts);
					}else{
						die "\nUnknown branch test!";
					}
					my $nti=NodeTestInfo->new();
					if(defined $test_pvalue){
						$nti->test_ref_node_name($ref_node_name);
						$nti->test_pvalue($test_pvalue);
						if($test_pvalue*$ntests<=$alpha){
							$nti->ref_node_name($name);
							push @sign_tests,$nti;
						}else{
							$nti->ref_node_name($ref_node_name);
						}
					}else{
						$nti->ref_node_name($ref_node_name);
					}
					$nodes_tests{$name}=$nti;
					print STDERR $p->report("\r%20b  ETA: %E", $i++);
				}
			}
	);
}
#print results
print "Branch test: $branch_tests[$branch_test]";
print "\nTotal number of tests: $ntests";
print "\nNumber of sign. tests: ";
print scalar @sign_tests;
#my @ord_nodes_tests=sort {$nodes_tests{$a}->test_pvalue <=> $nodes_tests{$b}->test_pvalue} keys %nodes_tests;
print "\nNodeID\tPNodeID";
if($branch_test){
	print "\tp-value\tRefNodeID\tKLD\tMaxSkewedGroupID\tKLD_MSG_Model";
}
print "\tGroupCounts:1";
for(my $i=2;$i<=$ngroups;$i++){
	print " $i";
}
if($branch_test){
	print "\tRefGropsFreqs:1";
	for(my $i=2;$i<=$ngroups;$i++){
		print " $i";
	}
}
print "\tBranchGropsFreqs:1";
for(my $i=2;$i<=$ngroups;$i++){
	print " $i";
}
if($branch_test){
	my %groups_freqs;
	struct BranchInfo => {
		test_pvalue => '$',
		max_skewed_group_id => '$',
		ref_groups_freqs => '@',
		obs_groups_freqs => '@',
		ref_KLD => '$',
		max_skewed_group_KLD => '$'
	};
	my %sign_changes;
	$tree->visit_breadth_first(
			-in => sub{
					my $node=shift;
					my $name=$node->get_name;
					unless($node->is_terminal){
						my $ra_counts=$groups_counts{$name}->{node};
						my $pname;
						unless($node->is_root){
							my $pnode=$node->get_parent;
							$pname=$pnode->get_name;
						}else{
							my $str=join " ",@{$ra_counts};
							print "\n$name\t\t\t\t\t\t\t$str\t\t";
						}
						my $nti=$nodes_tests{$name};
						if(defined $nti->test_pvalue){
							my $test_ref_node_name=$nti->test_ref_node_name;
							print "\n$name\t$pname\t".$nti->test_pvalue."\t$test_ref_node_name\t";
							unless(defined $groups_freqs{$test_ref_node_name}){
								my $ra_counts=$groups_counts{$test_ref_node_name}->{node};
								my $n1=0;
								foreach my $n(@{$ra_counts}){$n1+=$n;}
								$groups_freqs{$test_ref_node_name}=[@{$ra_counts}];
								for(my $i=0;$i<$ngroups;$i++){
									$groups_freqs{$test_ref_node_name}->[$i]/=$n1;
								}
							}
							my @fobs=(0) x $ngroups;
							my $max_skewed_group_id;
							my $min_kld;
							my $kld;
							my $n1=0;
							foreach my $n(@{$ra_counts}){$n1+=$n;}
							if($n1>10){
								for(my $i=0;$i<$ngroups;$i++){
									$fobs[$i]=$ra_counts->[$i]/$n1;
								}
								my @fobs_pcount;
								my @fref;
								my @gridx;
								for(my $i=0;$i<$ngroups;$i++){
									#remove empty groups
									if($groups_counts{$test_ref_node_name}->{node}->[$i]>0){
										push @fref, $groups_freqs{$test_ref_node_name}->[$i];
										push @gridx,$i+1;
									}
								}
								my $ngr=@gridx;
								for(my $i=0;$i<$ngroups;$i++){
									push @fobs_pcount,($ra_counts->[$i]+1)/($n1+$ngr) if $groups_counts{$test_ref_node_name}->{node}->[$i]>0;
								}
								if($ngr>2){
									my @kld=calc_KLD_for_prop_models(\@fref,\@fobs_pcount,$ngr);
									$min_kld=$kld[0];
									$kld=$kld[0];
									for(my $i=1;$i<=$ngr;$i++){
										if(defined $kld[$i]){
											if($kld[$i]<$min_kld){
												$min_kld=$kld[$i];
												$max_skewed_group_id=$i;
											}
										}else{
											$max_skewed_group_id=undef;
											last;
										}
									}
									if(defined $max_skewed_group_id){
										$max_skewed_group_id=$gridx[$max_skewed_group_id-1];
									}
								}
							}
							my $str=sprintf("%.6f",$kld);
							$str.="\t$max_skewed_group_id\t";
							$str.=sprintf("%.6f",$min_kld) if defined $max_skewed_group_id;
							$str.="\t".join " ",@{$ra_counts};
							$str.="\t".join " ",@{$groups_counts{$test_ref_node_name}->{node}};
							$str.="\t".join " ",@{$groups_counts{$name}->{branch}};
							#$str.="\t";
							#for(my $i=0;$i<$ngroups;$i++){
							#	my $fi=$groups_freqs{$test_ref_node_name}->[$i];
							#	$fi=sprintf("%.4f",$fi);
							#	$str.=$fi;
							#	$str.=" " if $i<$ngroups-1;
							#}
							print $str;
							if($nti->test_pvalue*$ntests<=$alpha){
								my $brinf=BranchInfo->new();
								$brinf->test_pvalue($nti->test_pvalue);
								$brinf->max_skewed_group_id($max_skewed_group_id);
								@{$brinf->ref_groups_freqs}=@{$groups_freqs{$test_ref_node_name}};
								@{$brinf->obs_groups_freqs}=@fobs;
								$brinf->ref_KLD($kld);
								$brinf->max_skewed_group_KLD($min_kld);
								$sign_changes{$name}=$brinf;
							}
						}
					}
				}
		);
	#print in FigTree format
	my %group2color;
	my %branch_color;
	if(defined $colors_fn){
		open INPF, "<$colors_fn" or die "\nUnable to open input file: $colors_fn!";
		while(<INPF>){
			$_=$` if(/#/);
			chomp;
			my @line=split '\t';
			for(my $i=0;$i<2;$i++){
				$line[$i]=~s/^\s+//;
				$line[$i]=~s/\s+$//;
			}
			if($line[0]>0&&$line[1]=~/\S+/){
				$group2color{$line[0]}=$line[1];
			}
		}
		close INPF;
		my $ncolors=keys %group2color;
		foreach my $name(keys %sign_changes){
			my $grid=0;
			$grid=$sign_changes{$name}->max_skewed_group_id if $ngroups>2;
			if(defined $grid){
				my $color=$group2color{$grid};
				if(defined $color){
					$branch_color{$name}=$color;
				}elsif($ncolors<$ngroups){
					$branch_color{$name}=$group2color{1};
				}else{
					die "\nUndefined color for the group ID=$grid"
				}
			}
		}
		warn "\nNumber of groups $ngroups is greater than the number of colors $ncolors:\n\tthe color of the first group will be used for all branches!" if $ngroups>$ncolors;
	}
	my %branch_pvalue;
	my %max_skewed_group;
	my %branch_freq_change;
	my %ref_KLD;
	my %max_skewed_KLD;
	my %obs_groups_freqs;
	foreach my $name(keys %sign_changes){
		my $grid=$sign_changes{$name}->max_skewed_group_id;
		$branch_pvalue{$name}=$sign_changes{$name}->test_pvalue;
		$ref_KLD{$name}=sprintf("%.4f",$sign_changes{$name}->ref_KLD);
		my @fobs=@{$sign_changes{$name}->obs_groups_freqs};
		for(my $i=0;$i<$ngroups;$i++){
			$fobs[$i]*=100;
			$fobs[$i]=sprintf("%.1f",$fobs[$i]);
		}
		my $str=join " ",@fobs;
		$obs_groups_freqs{$name}=$str;
		if(defined $grid){
			$max_skewed_group{$name}=$grid;
			$max_skewed_KLD{$name}=sprintf("%.4f",$sign_changes{$name}->max_skewed_group_KLD);
			$branch_freq_change{$name}=$sign_changes{$name}->obs_groups_freqs($grid-1)/$sign_changes{$name}->ref_groups_freqs($grid-1)-1.0;
			$branch_freq_change{$name}=sprintf("%.4f",$branch_freq_change{$name});
		}
	}
	my %tree2str_args;
	$tree2str_args{groups_freqs}=\%obs_groups_freqs;
	$tree2str_args{ref_KLD}=\%ref_KLD;
	$tree2str_args{pvalue}=\%branch_pvalue;
	if(defined $colors_fn){
		$tree2str_args{color}=\%branch_color;
	}
	if($ngroups>2){
		$tree2str_args{MSG_ID}=\%max_skewed_group;
		$tree2str_args{MSG_KLD}=\%max_skewed_KLD;
		$tree2str_args{freq_change} = \%branch_freq_change;
	}
	my $str=tree2str($tree,%tree2str_args);
	my ($basename,$dir,$ext) = fileparse($xparr_fn,'\.[^\.]*$');
	my $out_figtree_fn=$dir.$basename.".$branch_tests[$branch_test].site_groups.tre";
	open OPF, ">$out_figtree_fn" or die "\nUnable to open output file $out_figtree_fn!";
	print OPF "#NEXUS\nbegin trees;";
	print OPF "\n\ttree tree_pairs = [&R] ";
	print OPF $str;
	print OPF "\nend;";
	close OPF;
}else{
	$tree->visit_breadth_first(
			-in => sub{
					my $node=shift;
					my $name=$node->get_name;
					unless($node->is_terminal){
						my $ra_counts=$groups_counts{$name}->{node};
						my $str.="\t".join " ",@{$ra_counts};
						my $pname;
						unless($node->is_root){
							my $pnode=$node->get_parent;
							$pname=$pnode->get_name;
							$str.="\t".join " ",@{$groups_counts{$name}->{branch}};
							print "\n$name\t$pname";
						}else{
							$str.="\t";
							print "\n$name\t";
						}
						print $str;
					}
				}
		);
}