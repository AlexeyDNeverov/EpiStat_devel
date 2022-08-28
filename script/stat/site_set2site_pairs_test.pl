#!/usr/bin/env perl
#This script makes contingency table for overlap of two sets of pairs: the set of selected pairs (e. g., best predicted pairs) and a set of pairs with one or both sites from the specified site set.
#The script uses the xparr_enrich_test.pl parameter file.
#Usage: [options] parameters_file selected_pairs_file
#Params:
#[options] - Options have priorities over the parameter file
#	-l --link_logic <and|or> - Sets the 'PairsLinkLogic' parameter.
#<parameters> - The file with script params.
#		Contents:
#		[MutationNumbersFilter]="sites"|"pairs" - Apply filter of mutations' numbers on sites or site pairs. No default value.
#		[PairsOrdering]="0"|"1" - The type of site pairs: 0 - unordered or 1 - ordered
#			Deault="1"
#		[EpiStatistics]="FN" - file with epistatic statistics for pairs. Is require to apply the MutationNumbersFilter="pairs" filter.
#		[BGR_SiteMinMutations]="<uint>" - Minimal number of mutations in background sites
#		[FGR_SiteMinMutations]="<uint>" - Minimal number of mutations in foreground sites
#		[BgrSelectedSites="FN"] - If the option is accounted, only sites in the background which are listed in the specified file are used for FDR estimation.
#		[FgrSelectedSites="FN"] - If the option is accounted, only sites in the foreground which are listed in the specified file are used for FDR estimation.
#		[BothSitesSelected="0|1"] Defines the way to apply selected site filter for intragenic site pairs. No default value.
#				!!!The parameter is required if PairsType="intragene"
#				0 - at least one site in a pair has to be in one of selected site subsets: 'BgrSelectedSites' or 'FgrSelectedSites'
#				1 - both sites in a pair have to be in one of selected site subsets: 'BgrSelectedSites' or 'FgrSelectedSites'
#		PairsType=("intra[gene]"|"inter[gene]")
#		SitePairs="FN" - list of site pairs for that the epistatic statistics were calculated
#		TestSiteSet="FN" - the set of sites for enrichment testing 
#		LinkSiteSetWith="<backgr[ound]|target|pairs>" Search sites from the site set through the background or the target sites.
#		[ConvertSiteNumbering="FN"] - the table converting site coordinates from the XPARR file into coordinates in the site set
#		[PairsLinkLogic="<and|or>"] No defaul value. Used only with LinkSiteSetWith="pairs". 
#			Logic function that links a specified site set with a site pair: 'and' - with both sites, 'or' - at least one site in a pair.
#<selected_pairs>
#		<header_string>?
#		^<bgr_site>\t<fgr_site>
#			header_string=<column_name>\t<column_name>(\t<column_name>)*
#			<bgr_site>=uint
#			<fgr_site>=uint

use strict;
die "Set up the system variable \$EPISTAT_HOME" unless defined $ENV{EPISTAT_HOME};
use lib "$ENV{EPISTAT_LIB}";
#use File::Basename;
#my ($name,$dir,$ext) = fileparse($ARGV[0],'\.?[^\.]*$');
use Getopt::Std;
use SitePairMatrix;
use Class::Struct;
use AssociationStatistics::SiteMutationNumberFilter;
use AssociationStatistics::PairMutationNumberFilter;
use AssociationStatistics::SelectedSitesFilter;

my $fishers_test_r="$ENV{EPISTAT_HOME}/stat/fishers_exact_test.R";

my %args;
if(!getopts('l:',\%args)){
	die "\nError in option string!";
}

struct SiteFilterParams=>{
	mnumber_filter_type => '$', #'sites'|'pairs'
	pairs_ordering => '$', #'ord' - 1|'unord' - 0
	f_intragene => '$',
	site_pairs_fn => '$',
	epistat_ext => '$',
	obs_epistat_fn => '$',
	bgr_site_min_muts => '$',
	fgr_site_min_muts => '$',
	bg_filtered_sites_fn => '$',
	fg_filtered_sites_fn => '$',
	f_both_sites_selected => '$'
};

my $filters_params;
my $site_pair_mtx;
my $selected_sites_filter;

my ($pairs_fn,$stat_fn,$ordstat_fn,$fake_fn,$sites_fn,$cds_fn,$convert_sites_fn);
my $f_intragene;
my $pairs_ordering=1;
my $link_idx;
my $link_logic;
my $stat_type;
my $locus_name;

open INFILE, "$ARGV[0]" or die "\nUnable to open input file: $ARGV[0]!";
while(<INFILE>){
	$_=$` if(/#/);
	chomp;
	if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
		my ($key,$value)=($1,$2);
		if($key eq "SitePairs"){
			$pairs_fn=$value;
		}elsif($key eq "PairsType"){
			if($value=~m/^intra(gene)*/i){
				$f_intragene=1;
			}elsif($value=~m/^inter(gene)*/i){
				$f_intragene=0;
			}else{
				die "\nUnknown value ($value) for the PairsType parameter in file $ARGV[0]!"; 
			}
		}elsif($key eq "EpiStatistics"){
			$stat_fn=$value;
		}elsif($key eq "TestSiteSet"){
			$sites_fn=$value;
		}elsif($key eq "ConvertSiteNumbering"){
			$convert_sites_fn=$value;
		}elsif($key eq	"LinkSiteSetWith"){
			if($value=~/^backgr/i){
				$link_idx=0;
			}elsif($value=~/^target/i){
				$link_idx=1;
			}elsif($value=~/^pairs/i){
				#intragenic epistasis
				$link_idx=2;
			}else{
				die "\nOnly 'background','target' or 'pairs' values allowed for the 'LinkSiteSetWith' parameter!";
			}
		}elsif($key eq	"PairsLinkLogic"){
			if($value=~/^and/i){
				$link_logic=1;
			}elsif($value=~/^or/i){
				$link_logic=0;
			}else{
				die "\nOnly 'AND' or 'OR' values allowed for the 'PairsLinkLogic' parameter!";
			}
		}elsif($key eq "MutationNumbersFilter"){
			die "The parameter 'MutationNumbersFilter' must be 'sites' OR 'pairs'!" unless $value=~m/^sites|pairs$/i;
			$filters_params=SiteFilterParams->new unless defined $filters_params;
			$filters_params->mnumber_filter_type(lc($value));
		}elsif($key eq "PairsOrdering"){
			if($value==0||$value==1){
				$pairs_ordering=$value;
			}else{
				die "\nUnknown value ($value) for the 'PairsOrdering' parameter in file $ARGV[0]!"; 
			}
		}elsif($key eq "BGR_SiteMinMutations"){
			$filters_params=SiteFilterParams->new unless defined $filters_params;
			die "\nUINT value for the parameter 'BGR_SiteMinMutations' is expected!" unless $value=~/^\d+$/;
			$filters_params->bgr_site_min_muts($value);
		}elsif($key eq "FGR_SiteMinMutations"){
			$filters_params=SiteFilterParams->new unless defined $filters_params;
			die "\nUINT value for the parameter 'FGR_SiteMinMutations' is expected!" unless $value=~/^\d+$/;
			$filters_params->fgr_site_min_muts($value);
		}elsif($key eq "BgrSelectedSites"){
			$filters_params=SiteFilterParams->new unless defined $filters_params;
			$filters_params->bg_filtered_sites_fn($value);
		}elsif($key eq "FgrSelectedSites"){
			$filters_params=SiteFilterParams->new unless defined $filters_params;
			$filters_params->fg_filtered_sites_fn($value);
		}elsif($key eq "BothSitesSelected"){
			die "\nOnly 0 or 1 values are allowed for the 'BothSitesSelected' parameter!" unless $value=~m/^[01]$/;
			$filters_params=SiteFilterParams->new unless defined $filters_params;
			$filters_params->f_both_sites_selected($value);
		}
	}
}
close INFILE;
if(defined($args{l})){
	if($args{l}=~/and/i){
		$link_logic=1;
	}elsif($args{l}=~/or/i){
		$link_logic=0;
	}else{
		die "\nOnly 'AND' or 'OR' values allowed for the 'PairsLinkLogic' parameter!";
	}
}

die "\nNo default value for the 'LinkSiteSetWith' parameter!" unless defined $link_idx;
die "\nNo default value for the 'PairsLinkLogic' parameter!" unless defined($link_logic)||($link_idx<2);
die "\nThe pairs type was not defined intra- or intergenic!" unless defined $f_intragene;
#die "\nUndefined whether site pairs are ordered or not!" unless defined $pairs_ordering;
die "\nUnordered site pairs could be defined only for intragenic pairs!" unless ($pairs_ordering)||($f_intragene);
die "\nLinkSiteSetWith='pairs' is allowed only for intragenic pairs!" if (!$f_intragene)&&($link_idx==2);

if(defined $filters_params){
	$filters_params->obs_epistat_fn($stat_fn);
	$filters_params->site_pairs_fn($pairs_fn);
	$filters_params->f_intragene($f_intragene);
	$filters_params->pairs_ordering($pairs_ordering);
	die "\nUnable to initialize a site pair filter!" unless defined($filters_params->site_pairs_fn)&&defined($filters_params->f_intragene);
	my $infile_reading_mode=0;
	if(defined $filters_params->mnumber_filter_type){
		die "\nUnable to initialize a site pair filter!" unless $filters_params->fgr_site_min_muts>0||$filters_params->bgr_site_min_muts>0;
		if($filters_params->mnumber_filter_type() eq 'sites'){
			if($filters_params->bgr_site_min_muts>0){
				$infile_reading_mode++ unless $infile_reading_mode==3;
			}
			if($filters_params->fgr_site_min_muts>0){
				$infile_reading_mode+=2 unless $infile_reading_mode>=2;
			}
		}else{
			die "\nUnknown value for the \"MutationNumbersFilter\" parameter!" unless($filters_params->mnumber_filter_type() ne 'pairs');
		}
	}elsif(defined($filters_params->fgr_site_min_muts)||defined($filters_params->bgr_site_min_muts)){
		die "\nUnable to initialize mutatin number filter: the \"MutationNumbersFilter\" parameter is undefined!"
	}
	$site_pair_mtx=SitePairMatrix->new($filters_params->site_pairs_fn,$filters_params->f_intragene,$infile_reading_mode);
	$site_pair_mtx->init_sites;
	if(defined($filters_params->bg_filtered_sites_fn)||defined($filters_params->fg_filtered_sites_fn)){
		die "\nThe way how to select site pairs is not completely defined. Please specify 'BothSitesSelected'!" unless defined($filters_params->f_both_sites_selected)||(!$filters_params->pairs_ordering);
		$selected_sites_filter=AssociationStatistics::SelectedSitesFilter->new($site_pair_mtx,0,$filters_params->bg_filtered_sites_fn,$filters_params->fg_filtered_sites_fn,$filters_params->pairs_ordering,$filters_params->f_both_sites_selected);
	}
}

sub symmetrize_filter{
	my ($sp_matrix,$ra_filter)=@_;
	for(my $i=0;$i<@{$ra_filter};$i++){
		my ($bgs,$tgs)=$sp_matrix->line2sites_idxs($i);
		my $isym=$sp_matrix->site_idx_pair2line($tgs,$bgs);
		$ra_filter->[$i]=1 if($ra_filter->[$i]||$ra_filter->[$isym]);
	}
}

my %site_pairs;
my %selected_pairs;
my %sites;
my %site_ext2int;

if($convert_sites_fn){
	open INPF, "<$convert_sites_fn" or die "\nUnable to open input file $convert_sites_fn!";
	while(<INPF>){
		chomp;
		s/^\s+//;
		s/\s+$//;
		my @line=split '\t';
		if($line[0]=~/\d+/){
			$site_ext2int{$line[1]}=$line[0] if $line[1]=~/\S+/;
		}
	}
	close INPF;
}
#reading site set
my %site_subset;
open INFILE, "<$sites_fn" or die "\nUnable to open input file: $sites_fn!";
while(<INFILE>){
	$_=~s/\s//g;
	my @tmp=split ",";
	foreach my $str(@tmp){
		if($str=~/(\d+)-(\d+)/){
			for(my $i=$1;$i<=$2;$i++){
				my $a=$i;
				$a=$site_ext2int{$i} if $convert_sites_fn;
				$site_subset{$a}=1 if defined $a;
			}
		}elsif($str=~/(\d+)/){
			my $a=$1;
			$a=$site_ext2int{$1} if $convert_sites_fn;
			$site_subset{$a}=1 if defined $a;
		}
	}
}
close INFILE;
die "\nThe testing subset of sites is empty!" unless scalar(keys %site_subset);

my @site_pair_filter=(1) x $site_pair_mtx->{NLINES};
if(defined $filters_params){
	my $site_pair_filter;
	if($filters_params->mnumber_filter_type eq 'pairs'){
		$site_pair_filter=AssociationStatistics::PairMutationNumberFilter->new($site_pair_mtx,$filters_params->obs_epistat_fn,
			$filters_params->bgr_site_min_muts,$filters_params->fgr_site_min_muts);
	}elsif($filters_params->mnumber_filter_type eq 'sites'){
		$site_pair_filter=AssociationStatistics::SiteMutationNumberFilter->new($site_pair_mtx,
			$filters_params->bgr_site_min_muts,$filters_params->fgr_site_min_muts);
	}
	@site_pair_filter=$site_pair_filter->apply(\@site_pair_filter);
	@site_pair_filter=$selected_sites_filter->apply(\@site_pair_filter) if defined $selected_sites_filter;
	symmetrize_filter($site_pair_mtx,\@site_pair_filter) unless $filters_params->pairs_ordering;
}

{
	open PAIRS, "<$pairs_fn" or die "\nUnable to open input file $pairs_fn!";
	<PAIRS>; #skip header line
	my $ntg=0;
	my $nbg=0;
	my $I=0;
	while(<PAIRS>){
		chomp;
		if(/^(\d+)\t(\d+)\t(\d+)\t(\d+)/){
			my ($bgs,$tgs)=($1,$2);
			next unless (!defined($filters_params))||$site_pair_filter[$I++];
			if((!$pairs_ordering)&&($f_intragene)){
				next unless $bgs<$tgs;
			}
			$site_pairs{"$bgs,$tgs"}=[($bgs,$tgs)];
		}
	}
	close PAIRS;
}
open INFILE, "<$ARGV[1]" or die "\nUnable to open a list of selected site pairs $ARGV[1]!";
while(<INFILE>){
	if(/^(\d+)\t(\d+)/){
		my ($bgs,$tgs)=($1,$2);
		$selected_pairs{"$bgs,$tgs"}=[($bgs,$tgs)] if defined $site_pairs{"$bgs,$tgs"};
	}
}
close INFILE;
die "\nThe subset of selected site pairs is empty! Check file format and applied site pair filters" unless scalar(keys %selected_pairs);

sub is_in_site_set{
	my ($bgs,$tgs,$link_idx,$link_logic,$rh_site_set)=@_;
	my $t=0;
	if($link_idx==0){
		$t=1 if defined $rh_site_set->{$bgs};
	}elsif($link_idx==1){
		$t=1 if defined $rh_site_set->{$tgs};
	}elsif($link_idx==2){
		if($link_logic==1){
			$t=1 if defined($rh_site_set->{$bgs})&&defined($rh_site_set->{$tgs});
		}else{
			$t=1 if defined($rh_site_set->{$bgs})||defined($rh_site_set->{$tgs});
		}
	}
	return $t;
}

sub gen_tempname{
	my $nchar=shift;
	my @chars = ( "A" .. "Z", "a" .. "z", 0 .. 9 );
	return join("", @chars[ map { rand @chars } ( 1 .. $nchar ) ]);
}

#calculate contingency table
my @ctable=(0,0,0,0);
foreach my $sp(keys %site_pairs){
	my ($bgs,$tgs)=@{$site_pairs{$sp}};
	if(defined $selected_pairs{$sp}){
		if(is_in_site_set($bgs,$tgs,$link_idx,$link_logic,\%site_subset)){
			$ctable[0]++;
		}else{
			$ctable[1]++;
		}
	}else{
		if(is_in_site_set($bgs,$tgs,$link_idx,$link_logic,\%site_subset)){
			$ctable[2]++;
		}else{
			$ctable[3]++;
		}
	}
}
print "\nrows: selected_pairs\t1\t0\ncols: in_site_set\t1\t0\n";
my $tmp_fname=gen_tempname(10);
open OPF, ">$tmp_fname" or die "\nUnable to open output file: $tmp_fname!";
#print OPF "selected_pairs_1_0\tin_site_set_1_0\n";
print OPF $ctable[0]."\t".$ctable[1]."\n".$ctable[2]."\t".$ctable[3]."\n";
close OPF;
my $str="$fishers_test_r $tmp_fname";
my $out_fn=$tmp_fname.".fishers_test.out";
$str.=">$out_fn";

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

system($str);
print_child_termination_status();
system("cat $out_fn");
unlink $tmp_fname;
unlink $out_fn;