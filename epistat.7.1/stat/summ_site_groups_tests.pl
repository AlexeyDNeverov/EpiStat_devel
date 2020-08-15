#!/usr/bin/env perl
#This script calculates overlap statistics for branches identified by 'xparr_site_groups_test.pl'
#https://cran.r-project.org/web/packages/SuperExactTest/vignettes/set_html.html
#Usage: [options] <input_files_list>
#Options:
#	-a <parent|node> - Specifies a definition of a split node for overlap testing. If parent is accounted then test is applied for parental nodes.
#	[-x] <xparr_fn> - A file name storing a tree in XPAR format. If parameter was specified the pvalues for overlaps are calculated. 
#		The null-model assumes that branch is selected with probability proportional to its' length
#<input_files_list> - list of 'xparr_site_groups_test.pl' output files
use strict;
use lib "$ENV{EPISTAT_LIB}";
#my $multinom_R="$ENV{EPISTAT_HOME}/stat/multisetIntersect.R";
my $multinom_R;
use Class::Struct;
use File::Basename;
use Getopt::Std;
use Bio::Phylo::IO;
use List::BinarySearch qw(binsearch_pos);
use Time::Progress;

my %args;
if(!getopts('x:a:',\%args)){
	die "\nError in option string!";
}
my $splits_type;
if(defined $args{a}){
	if($args{a}=~/^parent/i){
		$splits_type="parent";
	}elsif($args{a}=~/^node/i){
		$splits_type="node";
	}
}else{
	die "\nPlease specify, which nodes the overlap test should be applied (-a)?";
}
my $xpar_fn=$args{x};
my $tree;
if(defined $xpar_fn){
	print STDERR "\nStart of reading the tree file...";
	$tree=Bio::Phylo::IO->parse(
		-file => $xpar_fn,
		-format => 'adjacency');
#!!!The return value of Adjacency parser not corresponded to interface of the Bio::Phylo::IO->parse
	$tree=$tree->[0]->first;
	print STDERR " done!";
}
my @infiles;
struct NodeInfo => {
	parent_id => '$',
	node_id => '$',
	branch_length => '$'
};

my $branch_num;
my $sign_thresh=0.05;
my $nsamples=10000;

my @node_info;
my @node_counts;
my @node_sign_counts;
my %node_name2idx;

open INPF, "<$ARGV[0]" or die "\nUnable to open output file: $ARGV[0]!";
while(<INPF>){
	chomp;
	s/^\s+//;
	s/\s+$//;
	if(/\S+/){
		push @infiles,$_;
	}
}
close INPF;
my $nfiles=@infiles;
for(my $i=0;$i<$nfiles;$i++){ 
	my $fn=$infiles[$i];
	open INPF,"<$fn" or die "\nUnable to open input file: $fn!";
	my $I=0;
	while(<INPF>){
		chomp;
		if(/Total number of tests:\s*(\d+)/){
			$branch_num=$1;
			$I++;
		}elsif(/Number of sign\. tests:\s*(\d+)/){
			$I++;
		}elsif(/Branch test:/){
			$I++;
		}
		last if $I==3;
	}
	die "\nWrong format of the input file $fn!" unless $I==3;
	while(<INPF>){
		chomp;
		if(/\S+/){
			my @line=split '\t';
			if($I>3){
				my $name=$line[0];
				my $pname=$line[1];
				my $pval=$line[2];
				if(($name ne "")&&($pname ne "")&&($pval ne "")){
					my $sidx;
					if(defined $node_name2idx{$name}){
						$sidx=$node_name2idx{$name};
					}else{
						$sidx=-1;
						my $si=NodeInfo->new();
						$si->parent_id($pname);
						$si->node_id($name);
						push @node_info,$si;
						push @node_counts,[(0) x $nfiles];
						push @node_sign_counts,[(0) x $nfiles];
						$node_name2idx{$name}=$#node_info;
					}
					$node_counts[$sidx]->[$i]++;
					if($pval*$branch_num<$sign_thresh){
						$node_sign_counts[$sidx]->[$i]++
					}
				}
			}
			$I++;
		}
	}
	close INPF;
}
{
#remove branches which were not tested for all input files
	my @tmp_ninfo;
	my @tmp_nscounts;
	foreach my $name (keys %node_name2idx){
		my $sidx=$node_name2idx{$name};
		my $i=0;
		foreach my $f(@{$node_counts[$sidx]}){
			$i++ if $f;
		}
		if($i==$nfiles){
			push @tmp_ninfo,$node_info[$sidx];
			push @tmp_nscounts,$node_sign_counts[$sidx];
		}
	}
	@node_info=@tmp_ninfo;
	@node_sign_counts=@tmp_nscounts;
	@node_counts=();
	%node_name2idx=();
}
sub get_parent_split_counts{
	my ($ra_node_info,$ra_sign_counts,$nfiles,$rh_splits)=@_;
	%{$rh_splits}=();
	for(my $i=0;$i<@{$ra_node_info};$i++){
		my $pname=$ra_node_info->[$i]->parent_id;
		if(!defined $rh_splits->{$pname}){
			$rh_splits->{$pname}=[];
			@{$rh_splits->{$pname}}=(0) x $nfiles;
		}
		for(my $j=0;$j<$nfiles;$j++){
			$rh_splits->{$pname}->[$j]++ if $ra_sign_counts->[$i]->[$j];
		}
	}
	return scalar(keys(%{$rh_splits}));
}

sub get_node_split_counts{
	my ($ra_node_info,$ra_sign_counts,$nfiles,$rh_splits)=@_;
	%{$rh_splits}=();
	for(my $i=0;$i<@{$ra_node_info};$i++){
		my $name=$ra_node_info->[$i]->node_id;
		if(!defined $rh_splits->{$name}){
			$rh_splits->{$name}=[];
			@{$rh_splits->{$name}}=(0) x $nfiles;
		}
		for(my $j=0;$j<$nfiles;$j++){
			$rh_splits->{$name}->[$j]++ if $ra_sign_counts->[$i]->[$j];
		}
	}
	return scalar(keys(%{$rh_splits}));
}

sub calc_overlaps_counts{
	my ($rh_splits,$nfiles)=@_;
	my @counts=(0) x ($nfiles+1);
	foreach my $name(keys %{$rh_splits}){
		my $si=$rh_splits->{$name};
		my $sn=0;
		for(my $i=0;$i<$nfiles;$i++){
			$sn++ if $si->[$i];
		}
		$counts[$sn]++;
	}
	return @counts;
}
#calculate statistics of subsets of significant splits overlaps
my %splits;
my $tot_tests;
if($splits_type eq "parent"){
	$tot_tests=get_parent_split_counts(\@node_info,\@node_sign_counts,$nfiles,\%splits);
}elsif($splits_type eq "node"){
	$tot_tests=get_node_split_counts(\@node_info,\@node_sign_counts,$nfiles,\%splits);
#test:
#my $n=0;
#foreach my $node(@{$tree->get_internals}){
#	my $name=$node->get_name;
#	next unless defined $splits{$name};
#	my $k=0;
#	foreach my $snode(@{$node->get_sisters}){
#		next if $snode->is_terminal;
#		$k++ if defined $splits{$snode->get_name};
#	}
#	$n++ if $k==1;
#}
#print "$n\t$tot_tests";
#exit;
#end test
}else{
	die "\nUnknown a split type!";
}
my @sign_tests=(0) x $nfiles;
my @overlap_counts=(0) x ($nfiles+1);
my @upper_pvals=(0) x ($nfiles+1);
my @exp_overlap_counts=(0) x ($nfiles+1);
my @overlaps;
my @sign_sets;
for(my $i=0;$i<$nfiles;$i++){
	push @sign_sets, [];
}

foreach my $name(keys %splits){
	my $si=$splits{$name};
	my $sn=0;
	for(my $i=0;$i<$nfiles;$i++){
		if($si->[$i]){
			$sn++;
			$sign_tests[$i]++;
		}
	}
	$overlap_counts[$sn]++;
	
	if($sn>0){
		$overlaps[$sn]={} unless defined $overlaps[$sn];
		$overlaps[$sn]->{$name}=[] unless defined $overlaps[$sn]->{$name};
		for(my $i=0;$i<$nfiles;$i++){
			if($si->[$i]){
				push @{$sign_sets[$i]},$name;
				push @{$overlaps[$sn]->{$name}},$i+1;
			}
		}
	}
}
for(my $i=$nfiles-1;$i>=0;$i--){
	$overlap_counts[$i]+=$overlap_counts[$i+1];
}

#/////////////////////////////////////////////////////////
if(defined $tree){
	my @node_info_tmp;
	foreach my $node(@{$tree->get_internals}){
		my $name=$node->get_name;
		my $l=0;
		if($splits_type eq "node"){
			$l=$node->get_branch_length;
		}elsif($splits_type eq "parent"){
			foreach my $chnode(@{$node->get_children}){
				$l+=$chnode->get_branch_length;
			}
		}
		if(defined $splits{$name}){
			my $ni=NodeInfo->new();
			$ni->node_id($name);
			$ni->parent_id($node->get_parent->get_name) unless $node->is_root;
			$ni->branch_length($l);
			push @node_info_tmp,$ni;
		}
	}
	my @fake_sign_counts;
	my @clength;
	for(my $i=0;$i<@node_info_tmp;$i++){
		my $name=$node_info_tmp[$i]->node_id;
		my $l=$node_info_tmp[$i]->branch_length;
		push @clength,$l;
		$clength[-1]+=$clength[-2] if $i>0;
		push @fake_sign_counts,[(0) x $nfiles];
	}
	my $I=0;
	print STDERR "\nStart calculation of pvalues for the overlap counts:";
	my $p = Time::Progress->new(min => 1, max => $nsamples);
	while($I<$nsamples){
		for(my $i=0;$i<$nfiles;$i++){
			my $j=0;
			while($j<$sign_tests[$i]){
				my $smpl=rand $clength[-1];
				my $ind=binsearch_pos {$a <=> $b} $smpl,@clength;
				if($fake_sign_counts[$ind]->[$i]==0){
					$fake_sign_counts[$ind]->[$i]++;
					$j++;
				}
			}
		}
		my %fake_splits;
		#If $splits_type eq "parent" we sampled parental splits initially!
		get_node_split_counts(\@node_info_tmp,\@fake_sign_counts,$nfiles,\%fake_splits);
		my @fake_overlap_counts=calc_overlaps_counts(\%fake_splits,$nfiles);
		for(my $i=$nfiles-1;$i>=0;$i--){
			$fake_overlap_counts[$i]+=$fake_overlap_counts[$i+1];
		}
		for(my $i=0;$i<=$nfiles;$i++){
			$upper_pvals[$i]++ if $fake_overlap_counts[$i]>=$overlap_counts[$i];
			$exp_overlap_counts[$i]+=$fake_overlap_counts[$i];
		}
		for(my $i=0;$i<@fake_sign_counts;$i++){
			for(my $j=0;$j<$nfiles;$j++){
				$fake_sign_counts[$i]->[$j]=0;
			}
		}
		$I++;
		print STDERR $p->report("\r%20b  ETA: %E", $I);
	}
	for(my $i=0;$i<=$nfiles;$i++){
		$upper_pvals[$i]/=$nsamples;
		$exp_overlap_counts[$i]/=$nsamples;
	}
}

print "The type of splits: $splits_type";
print "\nTotal number of tests performed: $tot_tests";
print "\nNumbers of sign. tests in each input files:";
for(my $i=0;$i<$nfiles;$i++){
	print " $sign_tests[$i]";
}
print "\nMultiplicity of significant tests overlaps (x#files)";
for(my $i=0;$i<$nfiles+1;$i++){
	print " $i";
}
print "\nPvalues of significant tests overlaps" if defined $tree;
print "\n#Overlaps\t";
for(my $i=0;$i<$nfiles+1;$i++){
	print " $overlap_counts[$i]";
}
if(defined $tree){
	print "\n#ExpOverlaps\t";
	for(my $i=0;$i<$nfiles+1;$i++){
		print " $exp_overlap_counts[$i]";
	}
	print "\nPvalues\t";
	for(my $i=0;$i<$nfiles+1;$i++){
		print " $upper_pvals[$i]";
	}
}
print "\nOverlaps:\nMultiplicity\tSplitId\tInFilesIds";
for(my $i=1;$i<@infiles+1;$i++){
	my @tmp;
	foreach my $pname(keys %{$overlaps[$i]}){
		push @tmp,[($i,$pname)];
		push @{$tmp[-1]},@{$overlaps[$i]->{$pname}};
	}
	@tmp=sort 
			{
				for(my $i=0;$i<$a->[0];$i++){
					my $t=$a->[2+$i]<=>$b->[2+$i];
					return $t if $t!=0;
				}
				return 0;
			} @tmp;
	foreach my $line(@tmp){
		my $str=join "\t", @{$line};
		print "\n$str";
	}
}

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
		}
	}
}
if(defined $multinom_R){
	my $tmp=gen_tempname(5).".summ_site_groups_test";
	open OPF, ">$tmp.tmp" or die "\nUnable to open output file: $tmp.tmp!";
	print OPF "$tot_tests\n";
	for(my $i=0;$i<$nfiles;$i++){
		my $str=join "\t",@{$sign_sets[$i]};
		print OPF "$str\n";
	}
	close OPF;
	my ($basename,$dir,$ext) = fileparse($ARGV[0],'\.[^\.]*$');
	my $str="$multinom_R $tmp.tmp $basename>$basename.MultiSetIntersection.R.out";
	system($str);
	print_child_termination_status();
	unlink "$tmp.tmp"
}

