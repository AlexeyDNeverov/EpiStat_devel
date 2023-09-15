#!/usr/bin/env perl
#This script calculates pvalues for mean ranks for site pairs having sites from specified subset of sites and for a complementray subset of site pairs
#The script expects the results of application of 'xparr_enrich_test.pl' for a specified subset of sites as well as for random subsets of sites (with -r option)
#Run first example:
#	~/epistat/stat/xparr_enrich_test.pl -l or site_subset.xparr_enrich_test.prm>site_subset.or.out
#	seq 1 400|parallel ~/epistat/stat/xparr_enrich_test.pl "-l or -r site_subset.xparr_enrich_test.prm>site_subset.randomize/{}.or.out"
#Usage: options xparr_enrich_test_out_file
#options:
#	-d <DIR> - a path to folder with results of xparr_enrich_test for random site subsets
#	-s <STR> - a suffix of file names in the folder with random subsets results
use strict;
use File::Basename;
use Class::Struct;
use Getopt::Std;

struct TestInfo =>{
	obs_subset_rank => '$',
	obs_compl_rank => '$',
	fake_subset_rank => '$',
	fake_compl_rank => '$',
	upval => '$',
	lpval => '$',
};

my %args;
if(!getopts('d:s:',\%args)){
	die "\nError in option string!";
}
my $indir=$args{d};
die "\nThe folder with xparr_enrich_test results for random subsets of sites is not specified!" unless defined $indir;
my $sample_suff=$args{s};

sub parse_test_results{
	my ($in_fn,$out_test_info)=@_;
	open INPF,"<$in_fn" or die "\nUnable to open input file $in_fn!";
	while(<INPF>){
		if(/subset:/){
			if(/\(obs\)=\s*([0-9]*\.?[0-9]+?)/){
				$out_test_info->obs_subset_rank($1);
			}
			if(/\(<fake>\)=\s*([0-9]*\.?[0-9]+?)/){
				$out_test_info->fake_subset_rank($1);
			}
		}elsif(/complement:/){
			if(/\(obs\)=\s*([0-9]*\.?[0-9]+?)/){
				$out_test_info->obs_compl_rank($1);
			}
			if(/\(<fake>\)=\s*([0-9]*\.?[0-9]+?)/){
				$out_test_info->fake_compl_rank($1);
			}
		}elsif(/P\(s_exp<=s_obs\)=/){
			if(/P\(s_exp<=s_obs\)=\s*([0-9]*\.?[0-9]+?)/){
				$out_test_info->lpval($1);
			}
			if(/P\(s_exp>=s_obs\)=\s*([0-9]*\.?[0-9]+?)/){
				$out_test_info->upval($1);
			}
		}
	}
	close INPF;
}
my $ti=TestInfo->new();
parse_test_results($ARGV[0],$ti);
my @samples;
opendir(DIR, $indir) or die "can't opendir $indir: $!";
my $file;
while (defined($file = readdir(DIR))) {
	if(defined $sample_suff){
		my ($basename,$dir,$ext) = fileparse($file,'\..*$');
		if($ext eq $sample_suff){
			my $str=$indir."/".$file;
			if(-f $str){
				my $ti=TestInfo->new();
				parse_test_results($str,$ti);
				push @samples,$ti;
			}
		}
	}else{
		my $str=$indir."/".$file;
		if(-f $str){
			my $ti=TestInfo->new();
			parse_test_results($str,$ti);
			push @samples,$ti;
		}
	}
}
closedir(DIR);
my $n=@samples;
print "Samples num=$n";
print "\nSubset rank (obs)=".$ti->obs_subset_rank;
@samples=sort {$a->obs_subset_rank <=> $b->obs_subset_rank} @samples;
my $i=0;
my $j;
for(;$i<$n;$i++){
	last if $samples[$i]->obs_subset_rank>$ti->obs_subset_rank;
	$j=$i if !defined($j)&&($samples[$i]->obs_subset_rank==$ti->obs_subset_rank);
}
unless(defined $j){
	if($i<$n){
		$j=$i;
	}else{
		$j=$n;
	}
}
print "\nP(sample_subset_rank<=obs_subset_rank)=".sprintf("%.4f", $i/$n);
print "\nP(sample_subset_rank>=obs_subset_rank)=".sprintf("%.4f", ($n-$j)/$n);

print "\n\nCompl rank (obs)=".$ti->obs_compl_rank;
@samples=sort {$a->obs_compl_rank <=> $b->obs_compl_rank} @samples;
$i=0;
$j=undef;
for(;$i<$n;$i++){
	last if $samples[$i]->obs_compl_rank>$ti->obs_compl_rank;
	$j=$i if !defined($j)&&($samples[$i]->obs_compl_rank==$ti->obs_compl_rank);
}
unless(defined $j){
	if($i<$n){
		$j=$i;
	}else{
		$j=$n;
	}
}
print "\nP(sample_subset_rank<=obs_compl_rank)=".sprintf("%.4f", $i/$n);
print "\nP(sample_subset_rank>=obs_compl_rank)=".sprintf("%.4f", ($n-$j)/$n);
