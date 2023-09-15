#!/usr/bin/env perl
#This script makes a file required by 'run_epistat.pl' with associations of sites and phenotypes for conditioned search for epistasis
#The script changes format of representation of epistat pvalues project files of associations of sites and phenotypes
#Usage: options [pheno_labels]
#options:
#	-p <FN> - Specifies a *.site_pairs file
#	-r <phenotype_id(,phenotype_id)*> - Specifies the order of association statistics in the output file
#		<phenotype_id>=uint
#	-P <1|2> - Specifies where the phenotypes are: 1 - in the background slot of the site_pairs file, 2 - in the foreground slot
#	[-l] <FN> - Specifies a *.lower.pvalue file
#	[-u] <FN> - Specifies a *.upper.pvalue file
#		!!!At lest one of '-l' or '-u' options must be specified
#pheno_labels - a TAB delimited file with text labels for phenotype_ids which used to ease of results inspection. If it is omitted, then phenotype_ids are used instead
use strict;
use lib "$ENV{EPISTAT_LIB}";
use Getopt::Std;
use File::Basename;
use SitePairMatrix;

my %args;
if(!getopts('p:r:l:u:P:',\%args)){
	die "\nError in option string!";
}

my $pairs_fn=$args{p};
die "\nUndefined a required file *.site_pairs! Use -p option for specification!" unless defined $pairs_fn;
my $sp_matrix=SitePairMatrix->new();
$sp_matrix->init_sites($pairs_fn,0);
my $npairs=$sp_matrix->{NLINES};
my @phen_labels;
my %phen2idx;
my $f_phen_in_bgr;
if(defined $args{r}){
	my @line=split ",",$args{r};
	my $i=0;
	foreach my $phen_id(@line){
		$phen2idx{$phen_id}=$i++;
	}
	unless(defined $ARGV[0]){
		@phen_labels=@line;
	}else{
		open INPF, "<$ARGV[0]" or die "\nUnable to open input file:$ARGV[0]!";
		my $n=0;
		while(<INPF>){
			chomp;
			s/\s+$//;
			my @tmp=split "\t";
			my $phi=$phen2idx{$tmp[1]};
			if(defined $phi){
				$phen_labels[$phi]=$tmp[0];
				$n++;
			}
		}
		close INPF;
		die "\nUnable to find labels for all phenotype codes in the file $ARGV[0]!" unless $n==@line; 
	}
}else{
	die "\nUndefined what associations need to be extracted! Use -r option for specification!";
}
if(defined $args{P}){
	if($args{P}=~/^[12]$/){
		if($args{P}==1){
			$f_phen_in_bgr=1;
		}else{
			$f_phen_in_bgr=0;
		}
	}else{
		die "\nUnpropper value has been set in the '-P' option! Allowed values are 1 or 2";
	}
}
die "\nNo default value for the '-P' option please specify!" unless defined $f_phen_in_bgr;

my $lpval_fn=$args{l};
my $upval_fn=$args{u};
die "\nAt least one files with pvalues is required!" unless defined($lpval_fn)||defined($upval_fn);
my @lpvals;
my @upvals;
if(defined $lpval_fn){
	open INPF,"<$lpval_fn" or die "\nUnable to open input file: $lpval_fn!";
	while(<INPF>){
		chomp;
		s/^\s+//;
		s/\s+$//;
		if(/\S+/){
			push @lpvals,$_;
		}
	}
	close INPF;
	die "\nUnexpected number of rows in the file: $lpval_fn!" unless $npairs==@lpvals;
}
if(defined $upval_fn){
	open INPF,"<$upval_fn" or die "\nUnable to open input file: $upval_fn!";
	while(<INPF>){
		chomp;
		s/^\s+//;
		s/\s+$//;
		if(/\S+/){
			push @upvals,$_;
		}
	}
	close INPF;
	die "\nUnexpected number of rows in the file: $upval_fn!" unless $npairs==@upvals;
}
unless(defined $lpval_fn){
	for(my $i=0;$i<$npairs;$i++){
		$lpvals[$i]=1.0-$upvals[$i];
	}
}
unless(defined $upval_fn){
	for(my $i=0;$i<$npairs;$i++){
		$upvals[$i]=1.0-$lpvals[$i];
	}
}

my @sites;
my @phens;
if($f_phen_in_bgr){
	$sp_matrix->get_sites(\@phens,\@sites);
}else{
	$sp_matrix->get_sites(\@sites,\@phens);
}
my %site2idx;
for(my $i=0;$i<@sites;$i++){
	$site2idx{$sites[$i]}=$i;
}
my @stats;
for(my $i=0;$i<$npairs;$i++){
	my ($bgs,$fgs)=$sp_matrix->line2site_pair($i);
	($bgs,$fgs)=($fgs,$bgs) if $f_phen_in_bgr;
	my $phi=$phen2idx{$fgs};
	if(defined $phi){
		my $si=$site2idx{$bgs};
		$stats[$si]=[] unless defined $stats[$si];
		$stats[$si]->[$phi]=$lpvals[$i]<$upvals[$i]?$lpvals[$i]:$upvals[$i];
	}
}
die "\nSome sites missed for  some of specified phenotypes!" unless @stats==@sites;
print "site\t".join("\t",@phen_labels);
for(my $i=0;$i<@sites;$i++){
	print "\n$sites[$i]\t".join("\t",@{$stats[$i]});
}