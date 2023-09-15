#!/usr/bin/env perl
#This script calculates enrichment statistics for the provided set of sites among leading and trailing epistatically interacting sites.
#Params:
#[options] - Options have priorities over the parameter file
#	-l --link_logic <and|or> - Sets the 'PairsLinkLogic' parameter.
#<parameters> - The file with script params.
#		Contents:
#		[MutationNumbersFilter]="sites"|"pairs" - Applies filter of mutations' numbers on sites or site pairs. No default value.
#		[PairsOrdering]="0"|"1" - The type of site pairs: 0 - unordered or 1 - ordered
#			Deault="1"
#		[Epistat_Ext]="STR" - An extation of fake files with epistatic statistics, e.g. ".stat.exp.22"
#		[EpiStatistics]="FN" - file with epistatic statistics for pairs
#		[FDR_Dir]="STR" - A path to a directory with random samples for a FDR estimate
#		[BGR_SiteMinMutations]="<uint>" - Minimal number of mutations in background sites
#		[FGR_SiteMinMutations]="<uint>" - Minimal number of mutations in foreground sites
#		[BgrSelectedSites="FN"] - If the option is accounted, only sites in the background which are listed in the specified file are used for FDR estimation.
#		[FgrSelectedSites="FN"] - If the option is accounted, only sites in the foreground which are listed in the specified file are used for FDR estimation.
#		[BothSitesSelected="0|1"] Defines the way to apply selected site filter for intragenic site pairs. No default value.
#				!!!The parameter is required if PairsType="intragene"
#				0 - at least one site in a pair has to be in one of selected site subsets: 'BgrSelectedSites' or 'FgrSelectedSites'
#				1 - both sites in a pair have to be in one of selected site subsets: 'BgrSelectedSites' or 'FgrSelectedSites'
#		PairsType=("intra[gene]"|"inter[gene]")
#		StatisticsType="<median|mean>" - The measure of central tendency used for comparison of distributions in the specified subset with the complement
#		SitePairs="FN" - list of site pairs for that the epistatic statistics were calculated
#		ObsOrderingStat="FN" - A file name with statistics for an observed dataset used for ranking of site pairs 
#		FakeOrderingStat="FN" - A file name with statistics for fake datasets used for ranking of site pairs 
#		SortingOrder="<(asc)ending|(desc)ending>" - A sorditing order of statistics from the best to the worst
#		TestSiteSet="FN" - the set of sites for enrichment testing 
#		LinkSiteSetWith="<backgr[ound]|target|pairs>" Search sites from the site set through the background or the target sites.
#		ConvertSiteNumbering="FN" - the table converting site coordinates from the XPARR file into coordinates in the site set
#		[PairsLinkLogic="<and|or>"] No defaul value. Used only with LinkSiteSetWith="pairs". 
#			Logic function that links a specified site set with a site pair: 'and' - with both sites, 'or' - at least one site in a pair.
#		[CDS="FN"] - The XML file with CDS statement.
#		[BackgroundSegmID]= "<1|2>" - Default=1. The CDS for background segment.
#		[LocusName="STR"] - The locus name containing sites from the testing site set.
#WARN: Results of analisis after application of any data fiter that selects different subsets of site pairs for observed or sampled datasets may be misleading!
use strict;
use lib "$ENV{EPISTAT_LIB}";
#use File::Basename;
#my ($name,$dir,$ext) = fileparse($ARGV[0],'\.?[^\.]*$');
use Getopt::Std;
use Class::Struct;
use DnaUtilities::CDS::annotation;
use SitePairMatrix;
use AssociationStatistics::SiteMutationNumberFilter;
use AssociationStatistics::PairMutationNumberFilter;
use AssociationStatistics::SelectedSitesFilter;
use List::Util qw/shuffle/;

my %args;
if(!getopts('l:r',\%args)){
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
	f_both_sites_selected => '$',
};

my $site_filter_params;
my $site_pair_mtx;
my $selected_sites_filter;

my $bgr_cds=1;
my $tgt_cds=1;

my ($pairs_fn,$stat_fn,$ordstat_fn,$fake_fn,$sites_fn,$cds_fn,$convert_sites_fn,$out_csv_fn);
my $f_pval_filters;
my $f_intragene;
my $pairs_ordering=1;
my $link_idx;
my $link_logic;
my $stat_type;
my $locus_name;
my $sort_order;
my $fdr_samples_dir;
my @sort_order=("ASCENDING","DESCENDING");

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
		}elsif($key eq "ObsOrderingStat"){
			$ordstat_fn=$value;
		}elsif($key eq "FakeOrderingStat"){
			$fake_fn=$value;
		}elsif($key eq "SortingOrder"){
			if($value=~/^desc(ending)*/i){
				$sort_order=1;
			}elsif($value=~/^asc(ending)*/i){
				$sort_order=0;
			}
		}elsif($key eq "StatisticsType"){
			if($value=~/^median/i){
				$stat_type=0;
			}elsif($value=~/^mean/i){
				$stat_type=1;
			}
		}elsif($key eq "TestSiteSet"){
			$sites_fn=$value;
		}elsif($key eq "BackgroundSegmID"){
			$bgr_cds=$value;
			die "\nWrong CDS id for background segment was accounted: $bgr_cds! Expects <1|2>." unless($bgr_cds==1||$bgr_cds==2);
		}elsif($key eq "ConvertSiteNumbering"){
			$convert_sites_fn=$value;
		}elsif($key eq "CDS"){
			$cds_fn=$value;
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
		}elsif($key eq	"LocusName"){
			$locus_name=$value;
		}elsif($key eq "MutationNumbersFilter"){
			die "The parameter 'MutationNumbersFilter' must be 'sites' OR 'pairs'!" unless $value=~m/^sites|pairs$/i;
			$site_filter_params=SiteFilterParams->new unless defined $site_filter_params;
			$site_filter_params->mnumber_filter_type(lc($value));
		}elsif($key eq "FDR_Dir"){
			$value.="/" unless $value=~/\/$/;
			$fdr_samples_dir=$value;
		}elsif($key eq "Epistat_Ext"){
			$site_filter_params=SiteFilterParams->new unless defined $site_filter_params;
			$site_filter_params->epistat_ext($value);
		}elsif($key eq "PairsOrdering"){
			if($value==0||$value==1){
				$pairs_ordering=$value;
			}else{
				die "\nUnknown value ($value) for the 'PairsOrdering' parameter in file $ARGV[0]!"; 
			}
		}elsif($key eq "BGR_SiteMinMutations"){
			$site_filter_params=SiteFilterParams->new unless defined $site_filter_params;
			die "\nUINT value for the parameter 'BGR_SiteMinMutations' is expected!" unless $value=~/^\d+$/;
			$site_filter_params->bgr_site_min_muts($value);
		}elsif($key eq "FGR_SiteMinMutations"){
			$site_filter_params=SiteFilterParams->new unless defined $site_filter_params;
			die "\nUINT value for the parameter 'FGR_SiteMinMutations' is expected!" unless $value=~/^\d+$/;
			$site_filter_params->fgr_site_min_muts($value);
		}elsif($key eq "BgrSelectedSites"){
			$site_filter_params=SiteFilterParams->new unless defined $site_filter_params;
			$site_filter_params->bg_filtered_sites_fn($value);
		}elsif($key eq "FgrSelectedSites"){
			$site_filter_params=SiteFilterParams->new unless defined $site_filter_params;
			$site_filter_params->fg_filtered_sites_fn($value);
		}elsif($key eq "BothSitesSelected"){
			die "\nOnly 0 or 1 values are allowed for the 'BothSitesSelected' parameter!" unless $value=~m/^[01]$/;
			$site_filter_params=SiteFilterParams->new unless defined $site_filter_params;
			$site_filter_params->f_both_sites_selected($value);
		}else{
			die "\nWrong parameter: $key!";
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
die "\nAn order of sorting statistisc is not accounted!" unless defined $sort_order;
die "\nNo default value for the StatisticsType parameter! Please specify 'median or 'mean'." unless defined $stat_type;
die "\nThe pairs' type was not defined intra- or intergenic!" unless defined $f_intragene;
#die "\nUndefined whether site pairs are ordered or not!" unless defined $pairs_ordering;

if(!$f_intragene){
	if($bgr_cds==2){
		$tgt_cds=1;
	}else{
		$tgt_cds=2;
	}
	if($link_idx<2){
		$bgr_cds--;
		$tgt_cds--;
	}else{
		die "\nThe same site set $sites_fn is linked to both sites in intergenic pairs!";
	}
}else{
	$tgt_cds=$bgr_cds;
}

my $nfake_datasets;
my @fdr_samples;
open FAKE, "<$fake_fn" or die "\nUnable to open input file $fake_fn!";
$_=<FAKE>;
if(/^#SAMPLE_ID:/){
	$_=~s/^#SAMPLE_ID://;
	chomp;
	@fdr_samples=split '\t';
	$nfake_datasets=@fdr_samples;
}else{
	die "\nThe header line with sample IDs is required for fake data samples in $fake_fn!";
}
close FAKE;
die "\nNo fake datasets!" if $nfake_datasets==0;
my $ra_site_pair_filter;
{
	my $infile_reading_mode=0;
	if(defined $site_filter_params){
		$site_filter_params->obs_epistat_fn($stat_fn);
		$site_filter_params->site_pairs_fn($pairs_fn);
		$site_filter_params->f_intragene($f_intragene);
		$site_filter_params->pairs_ordering($pairs_ordering);
		die "\nUnable to initialize a site pair filter!" unless defined($site_filter_params->site_pairs_fn)&&defined($site_filter_params->f_intragene);
		if(defined $site_filter_params->mnumber_filter_type){
			die "\nUnable to initialize a site pair filter!" unless $site_filter_params->fgr_site_min_muts>0||$site_filter_params->bgr_site_min_muts>0;
			die "\nUnordered site pairs could be defined only for intragenic pairs!" unless ($site_filter_params->pairs_ordering)||($site_filter_params->f_intragene);
			if($site_filter_params->mnumber_filter_type() eq 'pairs'){
				die "\nSuffix of fake files with epistatic statistics is not defined!" unless defined $site_filter_params->epistat_ext;
				die "\nThe 'FDR_Dir' parameter is required if 'MutationNumbersFilter'='pairs'!" unless defined $fdr_samples_dir;
			}elsif($site_filter_params->mnumber_filter_type() eq 'sites'){
				if($site_filter_params->bgr_site_min_muts>0){
					$infile_reading_mode++ unless $infile_reading_mode==3;
				}
				if($site_filter_params->fgr_site_min_muts>0){
					$infile_reading_mode+=2 unless $infile_reading_mode>=2;
				}
			}else{
				die "\nUnknown value for the \"MutationNumbersFilter\" parameter!";
			}
		}elsif(defined($site_filter_params->fgr_site_min_muts)||defined($site_filter_params->bgr_site_min_muts)){
			die "\nUnable to initialize mutatin number filter: the \"MutationNumbersFilter\" parameter is undefined!"
		}
	}
	if(defined($site_filter_params)){
		$site_pair_mtx=SitePairMatrix->new($pairs_fn,$f_intragene,$infile_reading_mode);
		$site_pair_mtx->init_sites;
		$ra_site_pair_filter=[];
		my $n=1;
		if((defined($site_filter_params)&&($site_filter_params->mnumber_filter_type eq 'pairs'))){
			$n+=$nfake_datasets;
		}
		for(my $i=0;$i<$n;$i++){
			$ra_site_pair_filter->[$i]=[];
			@{$ra_site_pair_filter->[$i]}=(1) x $site_pair_mtx->{NLINES};
		}
		if(defined $site_filter_params){
			if(defined($site_filter_params->bg_filtered_sites_fn)||defined($site_filter_params->fg_filtered_sites_fn)){
				die "\nThe way how to select site pairs is not completely defined. Please specify 'BothSitesSelected'!" unless defined($site_filter_params->f_both_sites_selected)||(!$site_filter_params->pairs_ordering);
				$selected_sites_filter=AssociationStatistics::SelectedSitesFilter->new($site_pair_mtx,0,$site_filter_params->bg_filtered_sites_fn,$site_filter_params->fg_filtered_sites_fn,$site_filter_params->pairs_ordering,$site_filter_params->f_both_sites_selected);
			}
		}
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

sub get_best_stat{
	my ($sort_order,$s1,$s2)=@_;
	if(defined($s1)&&defined($s2)){
		if($sort_order==0){
			return ($s1<=$s2)?$s1:$s2;
		}
		return ($s1>=$s2)?$s1:$s2;
	}
	return defined($s1)?$s1:$s2;
}

#information about detected pairs
struct SitePairInfo => {
	bg_site => '$',
	tg_site => '$',
	bg_nsubst => '$',
	tg_nsubst => '$',
	ordering_stats => '@'
};

struct SiteLocation =>{
	locus_name => '$',
	position => '$'
};

struct SiteInfo =>{
	location => '$',
	ordering_stats => '@',
	best_site_pair => '$'
};

my @site_pairs;
my %sites;

my @annotation;
my %site_ext2int;
if($cds_fn){
	#Reading the XML annotation file
	@annotation=DnaUtilities::CDS::annotation::read_annotation_xml($cds_fn);
	die "\nError: the subset: $sites_fn\n\t has linked with both background and foreground sites, thus common CDS for both is required!" if(@annotation>1&&$link_idx==2);
}elsif($convert_sites_fn){
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
if(exists $args{r}){
	#randomize subset of sites
	$site_pair_mtx=SitePairMatrix->new($pairs_fn,$f_intragene,0);
	$site_pair_mtx->init_sites;
	my $n=0;
	my @bgr_sites;
	my @fgr_sites;
	$site_pair_mtx->get_sites(\@bgr_sites,\@fgr_sites);
	foreach my $site(@bgr_sites){
		$n++ if defined $site_subset{$site};
	}
	@bgr_sites=shuffle @bgr_sites;
	@bgr_sites=@bgr_sites[0..$n-1];
	%site_subset=();
	for(my $i=0;$i<$n;$i++){
		$site_subset{$bgr_sites[$i]}=1;
	}
}

if(defined $site_filter_params){
	my $site_pair_filter;
	my $n=@{$ra_site_pair_filter};
	if($site_filter_params->mnumber_filter_type eq 'pairs'){
		$site_pair_filter=AssociationStatistics::PairMutationNumberFilter->new($site_pair_mtx,$site_filter_params->obs_epistat_fn,
			$site_filter_params->bgr_site_min_muts,$site_filter_params->fgr_site_min_muts);
	}elsif($site_filter_params->mnumber_filter_type eq 'sites'){
		$site_pair_filter=AssociationStatistics::SiteMutationNumberFilter->new($site_pair_mtx,
			$site_filter_params->bgr_site_min_muts,$site_filter_params->fgr_site_min_muts);
	}
	for(my $i=0;$i<$n;$i++){
		@{$ra_site_pair_filter->[$i]}=$site_pair_filter->apply($ra_site_pair_filter->[$i]);
		@{$ra_site_pair_filter->[$i]}=$selected_sites_filter->apply($ra_site_pair_filter->[$i]) if defined $selected_sites_filter;
		symmetrize_filter($site_pair_mtx,$ra_site_pair_filter->[$i]) unless $site_filter_params->pairs_ordering;
		if($site_filter_params->mnumber_filter_type eq 'pairs'){
			#update filter for the next fake sample
			if($i<$n-1){
				my $fake_fn=$fdr_samples_dir.$fdr_samples[$i].$site_filter_params->epistat_ext;
				$site_pair_filter=AssociationStatistics::PairMutationNumberFilter->new($site_pair_mtx,$fake_fn,
					$site_filter_params->bgr_site_min_muts,$site_filter_params->fgr_site_min_muts);
			}
		}
	}
}

sub is_filter_ok{
	my ($ra_filter,$site_pair_idx,$sample_idx)=@_;
	#'$sample_idx' = [0(observation),$nfake_datasets]
	#die "\nError is_filter_ok(): Undefined site pair index!" unless defined $site_pair_idx;
	if(defined $ra_filter){
		if(@{$ra_filter}==1){
			return $ra_filter->[0]->[$site_pair_idx];
		}
		return $ra_filter->[$sample_idx]->[$site_pair_idx] if defined $sample_idx;
		for(my $i=0;$i<@{$ra_filter};$i++){
			return 1 if $ra_filter->[$i]->[$site_pair_idx];
		}
		return 0;
	}
	return 1;
}

{
	open ORDSTAT, "<$ordstat_fn" or die "\nUnable to open input file $ordstat_fn!";
	open FAKE, "<$fake_fn" or die "\nUnable to open input file $fake_fn!";
	my $ntg=0;
	my $nbg=0;
	open PAIRS, "<$pairs_fn" or die "\nUnable to open input file $pairs_fn!";
	<PAIRS>; #skip header line
	for(my $I=0;!(eof(PAIRS)||eof(ORDSTAT)||eof(FAKE));$I++){
		my $str=<PAIRS>;
		my $spi=SitePairInfo->new();
		$spi->bg_site(0);
		$spi->tg_site(0);
		$spi->bg_nsubst(0);
		$spi->tg_nsubst(0);
		if($str=~/^(\d+)\t(\d+)\t(\d+)\t(\d+)/){
			$spi->bg_site($1);
			$spi->tg_site($2);
			$spi->bg_nsubst($3);
			$spi->tg_nsubst($4);
		}
		my $str=<ORDSTAT>;
		chomp $str;
		if($str=~/([-+]?\d*\.?\d+([eE]?[-+]?\d*))/){
			$spi->ordering_stats->[0]=$1 if is_filter_ok($ra_site_pair_filter,$I,0);
		}else{
			die "\nA value is expected for the line $I in $ordstat_fn!" unless $str=~/na/i;
		}
		$str=<FAKE>;
		while($str=~/^\s*#/){
			$str=<FAKE>;
		}
		chomp $str;
		my @fake_stats=split '\t', $str;
		die "\nThe numbers of columns in some rows of the '$fake_fn' file are not equal to the value in 1st row  = $nfake_datasets!" unless $nfake_datasets==@fake_stats;
		next unless defined($spi->bg_site)&&defined($spi->tg_site);
		next unless is_filter_ok($ra_site_pair_filter,$I);
		#die "\nUndefined ordering statistics\n\t$str\nin the file: $ordstat_fn!" unless defined $spi->ordering_stat;
		if($link_idx==2){
			for(my $i=0;$i<$nfake_datasets;$i++){
				if($fake_stats[$i]=~/([-+]?\d*\.?\d+([eE]?[-+]?\d*))/){
					$spi->ordering_stats->[$i+1]=$fake_stats[$i] if is_filter_ok($ra_site_pair_filter,$I,$i+1);
				}
			}
			my $bg_loc=SiteLocation::new();
			my $tg_loc=SiteLocation::new();
			if(defined $cds_fn){
				my ($site,$locus)=DnaUtilities::CDS::annotation::alnpos2annotation(0,$spi->bg_site,\@annotation);
				$bg_loc->position($site);
				$bg_loc->locus_name($locus) if(defined $locus_name);
				($site,$locus)=DnaUtilities::CDS::annotation::alnpos2annotation(0,$spi->tg_site,\@annotation);
				$tg_loc->position($site);
				$tg_loc->locus_name($locus) if(defined $locus_name);
			}else{
				$bg_loc->position($spi->bg_site);
				$tg_loc->position($spi->tg_site);
			}
			push @site_pairs,[($bg_loc,$tg_loc,$spi)];
		}else{
			my @key_sites;
			if($link_idx==0){
				$key_sites[0]=$spi->bg_site;
				$key_sites[1]=$spi->tg_site if (!$pairs_ordering)&&($f_intragene);
			}else{
				$key_sites[0]=$spi->tg_site;
				$key_sites[1]=$spi->bg_site if (!$pairs_ordering)&&($f_intragene);
			}
			foreach my $key_site(@key_sites){
				if(!defined $sites{$key_site}){
					my $site_info=SiteInfo::new();
					my $site_loc=SiteLocation::new();
					if(defined $cds_fn){
						my ($site,$locus)=DnaUtilities::CDS::annotation::alnpos2annotation($link_idx==0?$bgr_cds:$tgt_cds,$key_site,\@annotation);
						$site_loc->position($site);
						$site_loc->locus_name($locus) if(defined $locus_name);
					}else{
						$site_loc->position($key_site);
					}
					$site_info->location($site_loc);
					if(is_filter_ok($ra_site_pair_filter,$I,0)){
						$site_info->ordering_stats->[0]=$spi->ordering_stats(0);
						$site_info->best_site_pair($spi) if defined $spi->ordering_stats(0);
					}
					$sites{$key_site}=$site_info;
				}elsif(is_filter_ok($ra_site_pair_filter,$I,0)){
					my $bs=get_best_stat($sort_order,$sites{$key_site}->ordering_stats(0),$spi->ordering_stats(0));
					if(defined($bs)&&($bs!=$sites{$key_site}->ordering_stats(0))){
						$sites{$key_site}->ordering_stats->[0]=$bs;
						$sites{$key_site}->best_site_pair($spi);
					}
				}
				for(my $i=0;$i<$nfake_datasets;$i++){
					next unless is_filter_ok($ra_site_pair_filter,$I,$i+1);
					if($fake_stats[$i]=~/([-+]?\d*\.?\d+([eE]?[-+]?\d*))/){
						if(!defined $sites{$key_site}->ordering_stats->[$i+1]){
							$sites{$key_site}->ordering_stats->[$i+1]=$fake_stats[$i];
						}else{
							$sites{$key_site}->ordering_stats->[$i+1]=get_best_stat($sort_order,$sites{$key_site}->ordering_stats($i+1),$fake_stats[$i]);
						}
					}
				}
			}
		}
	}
	close PAIRS;
	close FAKE;
	close ORDSTAT;
}

sub is_site_in_subset{
	my ($rsite_loc,$rh_site_set)=@_;
	if(!defined($rsite_loc->locus_name)|| $rsite_loc->locus_name eq $locus_name){
		return 1 if defined($rh_site_set->{$rsite_loc->position});
	}
	return 0;
}

sub get_sites_keys{
	my ($hr_sites,$sample_idx)=@_;
	my @keys;
	foreach my $site(keys %{$hr_sites}){
		push @keys,$site if defined $hr_sites->{$site}->ordering_stats($sample_idx);
	}
	@keys=sort {
		return $hr_sites->{$b}->ordering_stats($sample_idx) <=> $hr_sites->{$a}->ordering_stats($sample_idx) if $sort_order;
		$hr_sites->{$a}->ordering_stats($sample_idx) <=> $hr_sites->{$b}->ordering_stats($sample_idx);
		} @keys;
	return @keys;
}

sub calc_ranks{
	my $ra_ord_stats_values=shift;
	my @ranks;
	my $val=$ra_ord_stats_values->[0];
	my $R=1;
	my $n=0;
	my $N=@{$ra_ord_stats_values};
	for(my $i=1;$i<$N;$i++){
		my $str=$ra_ord_stats_values->[$i];
		if($str==$val){
			$n++;
			$R+=$i+1;
		}else{
			$R/=($n+1) if $n;
			for(my $j=$i-$n-1;$j<$i;$j++){
				$ranks[$j]=$R;
			}
			$R=$i+1;
			$n=0;
			$val=$str;
		}
	}
	$R/=($n+1) if $n;
	for(my $j=$N-$n-1;$j<$N;$j++){
		$ranks[$j]=$R;
	}
	return @ranks;
}

sub calc_subset_statistics{
	my ($ra_ordstat,$ra_keys_subset,$ra_keys_compl,$stat_type)=@_;
	my $msubset;
	my $mcompl;
	my $n_subset=@{$ra_keys_subset};
	my $n_compl=@{$ra_keys_compl};
	if($stat_type==0){ #median
		if($n_subset){
			my $i=int($n_subset/2);
			$msubset=$ra_ordstat->[$ra_keys_subset->[$i]];
			if($n_subset%2==0){
				$msubset+=$ra_ordstat->[$ra_keys_subset->[$i-1]];
				$msubset/=2;
			}
		}
		if($n_compl){
			my $i=int($n_compl/2);
			$mcompl=$ra_ordstat->[$ra_keys_compl->[$i]];
			if($n_compl%2==0){
				$mcompl+=$ra_ordstat->[$ra_keys_compl->[$i-1]];
				$mcompl/=2;
			}
		}
	}elsif($stat_type==1){ #mean
		if($n_subset){
			$msubset=0;
			foreach my $i(@{$ra_keys_subset}){
				$msubset+=$ra_ordstat->[$i];
			}
			$msubset/=$n_subset;
		}
		if($n_compl){
			$mcompl=0;
			foreach my $i(@{$ra_keys_compl}){
				$mcompl+=$ra_ordstat->[$i];
			}
			$mcompl/=$n_compl;
		}
	}
	return ($msubset,$mcompl,$n_subset,$n_compl);
}

sub get_statistics_for_sites{
	my ($hr_sites,$idx)=@_;
	my @keys=get_sites_keys($hr_sites,$idx);
	return (undef,undef,0,0) unless scalar(@keys);
	my @ord_stats_values;
	foreach my $key(@keys){
		push @ord_stats_values,$hr_sites->{$key}->[2]->ordering_stats($idx);
	}
	#calculate ranks
	my @ranks=calc_ranks(\@ord_stats_values);
	my @keys_subset;
	my @keys_compl;
	for(my $i=0;$i<@keys;$i++){
		my $key=$keys[$i];
		if(is_site_in_subset $hr_sites->{$key}->location,\%site_subset){
			push @keys_subset,$i;
		}else{
			push @keys_compl,$i;
		}
	}
	warn "\nUnable to compare two sets for the sample $idx: one is empty!" unless @keys_subset>0&&@keys_compl>0;
	return calc_subset_statistics(\@ranks,\@keys_subset,\@keys_compl,$stat_type);
}

#returns indices of site pairs in array ordered according increase of ordering statistics for the observed (idx==0) or for a fake data
#for two complementary pairs the one with the best statistics is presented
sub get_pairs_keys{
	my ($ha_site_pairs,$idx)=@_;
	my %hkeys;
	my %pair_keys;
	for(my $i=0;$i<@{$ha_site_pairs};$i++){
		my $bgs=$ha_site_pairs->[$i]->[2]->bg_site;
		my $tgs=$ha_site_pairs->[$i]->[2]->tg_site;
		$pair_keys{$bgs.",".$tgs}=$i;
	}
	for(my $i=0;$i<@{$ha_site_pairs};$i++){
		my $val=$ha_site_pairs->[$i]->[2]->ordering_stats($idx);
		if(defined $val){
			my $bgs=$ha_site_pairs->[$i]->[2]->bg_site;
			my $tgs=$ha_site_pairs->[$i]->[2]->tg_site;
			my $j=$pair_keys{$tgs.",".$bgs}; #complementary pair
			if(defined $j){
				my $rval=$ha_site_pairs->[$j]->[2]->ordering_stats($idx);
				my $best_val=get_best_stat($sort_order,$val,$rval);
				$hkeys{$val==$best_val?$i:$j}=1 unless defined($hkeys{$i})||defined($hkeys{$j});
			}else{
				$hkeys{$i}=1;
			}
		}
	}
	return sort {
		return $ha_site_pairs->[$b]->[2]->ordering_stats($idx) <=> $ha_site_pairs->[$a]->[2]->ordering_stats($idx) if($sort_order);
		$ha_site_pairs->[$a]->[2]->ordering_stats($idx) <=> $ha_site_pairs->[$b]->[2]->ordering_stats($idx);
		} keys(%hkeys);
}

sub is_pair_in_subset{
	my ($bg_loc,$fg_loc,$rh_site_subset,$logand)=@_;
	if($logand){
		return is_site_in_subset($bg_loc,$rh_site_subset)&&
			is_site_in_subset($fg_loc,$rh_site_subset);
	}
	return is_site_in_subset($bg_loc,$rh_site_subset)||
			is_site_in_subset($fg_loc,$rh_site_subset);
}

sub get_statistics_for_pairs{
	my ($ha_site_pairs,$idx)=@_;
	my @keys=get_pairs_keys($ha_site_pairs,$idx);
	return (undef,undef,0,0) unless scalar(@keys);
	my @ord_stats_values;
	foreach my $key(@keys){
		push @ord_stats_values,$ha_site_pairs->[$key]->[2]->ordering_stats($idx);
	}
	#calculate ranks
	my @ranks=calc_ranks(\@ord_stats_values);
	my @keys_subset;
	my @keys_compl;
	for(my $i=0;$i<@keys;$i++){
		my $key=$keys[$i];
		if(is_pair_in_subset($ha_site_pairs->[$key]->[0],$ha_site_pairs->[$key]->[1],\%site_subset,$link_logic)){
			push @keys_subset,$i;
		}else{
			push @keys_compl,$i;
		}
	}
	warn "\nUnable to compare two sets for the sample $idx: one is empty!" unless @keys_subset>0&&@keys_compl>0;
	return calc_subset_statistics(\@ranks,\@keys_subset,\@keys_compl,$stat_type);
}

my @wresults;
my ($msubset,$mcompl,$nsubset,$ncompl);
if($link_idx<2){
	($msubset,$mcompl,$nsubset,$ncompl)=get_statistics_for_sites(\%sites,0);
}else{
	($msubset,$mcompl,$nsubset,$ncompl)=get_statistics_for_pairs(\@site_pairs,0);
}
my ($msubset_fake,$mcompl_fake)=(0,0);
for(my $i=1;$i<=$nfake_datasets;$i++){
	my @vals;
	if($link_idx<2){
		@vals=get_statistics_for_sites(\%sites,$i);
	}else{
		@vals=get_statistics_for_pairs(\@site_pairs,$i);
	}
	if(defined($vals[0])&&defined($vals[1])){
		$msubset_fake+=$vals[0];
		$mcompl_fake+=$vals[1];
		push @wresults,$vals[0]-$vals[1];
	}
}
my $nfakes=@wresults;
$msubset_fake/=$nfakes;
$mcompl_fake/=$nfakes;
@wresults=sort {$a <=> $b} @wresults;

my $diff;
$diff=$msubset-$mcompl if defined($msubset)&&defined($mcompl);
my $p_eless=-1;
my $p_egreater=-1;
if(defined $diff){
	for($p_eless=0;$p_eless<$nfakes;$p_eless++){
		$p_egreater=$p_eless if ($p_egreater==-1)&&($wresults[$p_eless]==$diff);
		last if($wresults[$p_eless]>$diff);
	}
	$p_egreater=$p_eless if $p_egreater==-1;
	$p_eless/=$nfakes;
	$p_eless=sprintf("%.4f",$p_eless);
	$p_egreater=$nfakes-$p_egreater;
	$p_egreater/=$nfakes;
	$p_egreater=sprintf("%.4f",$p_egreater);
}
print "Number of observations in the subset=$nsubset\n";
print "Number of observations in the complement=$ncompl\n";
print "SortingOrder=\"".$sort_order[$sort_order]."\"\n";
if($stat_type==0){
	$_="Median ";
}elsif($stat_type==1){
	$_="Mean ";
}
print $_;
print "subset: (obs)=$msubset\t(<fake>)=$msubset_fake\n";
print $_;
print "complement: (obs)=$mcompl\t(<fake>)=$mcompl_fake\n";
print "s=".$_." nominal statistics(subset) - Median nominal statistics(complement)";
print "\nP(s_exp<=s_obs)=$p_eless\tP(s_exp>=s_obs)=$p_egreater";
print "\nNSamples=".scalar($nfakes);


if($link_idx<2){
	my @keys=get_sites_keys(\%sites,0);
	if(@keys){
		print "\nObserved sites from the subset{";
		print "\nsite\tstat\tbest_pair:bgr_site\tfgr_site";
		for(my $i=0;$i<@keys;$i++){
			if(is_site_in_subset $sites{$keys[$i]}->location,\%site_subset){
				print "\n";
				my $stat=sprintf("%.4f",$sites{$keys[$i]}->ordering_stats(0));
				print "$keys[$i]\t$stat";
				my $bspi=$sites{$keys[$i]}->best_site_pair;
				if(defined $bspi){
					print "\t".$bspi->bg_site."\t".$bspi->tg_site;#."\t".$bspi->ordering_stats(0);
				}
			}
		}
		print "}\n";
	}
}else{
	my @keys=get_pairs_keys(\@site_pairs,0);
	if(@keys){
		print "\nObserved site pairs having ";
		if($link_logic){
			print "both sites";
		}else{
			print "at least one site";
		};
		print " from the subset{";
		print "\nbgr_site\tfgr_site\tstat";
		foreach my $i(@keys){
			if(is_pair_in_subset($site_pairs[$i]->[0],$site_pairs[$i]->[1],\%site_subset,$link_logic)){
				print "\n";
				my $stat=sprintf("%.4f",$site_pairs[$i]->[2]->ordering_stats(0));
				my $bgs=$site_pairs[$i]->[2]->bg_site;
				my $tgs=$site_pairs[$i]->[2]->tg_site;
				print "$bgs\t$tgs\t$stat";
			}
		}
		print "}\n";
	}
}
