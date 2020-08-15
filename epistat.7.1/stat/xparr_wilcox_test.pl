#!/usr/bin/env perl
#This script calculates enrichment statistics for the provided set of sites among leading and trailing epistatically interacting sites.
#A subset of sites or site pairs is compared against its complement by nonparametric U-test.
#Usage: [options] <parameters>
#[options] - Options have priorities over the parameter file
#	[-a] alternative=("t(wo.sided)"|"l(ess)"|"g(reater)") - Alternative hypothesis for R wilcox.test function
#	Default=two.sided
#	-l --link_logic <and|or> - Sets the 'PairsLinkLogic' parameter.
#Params:
#<parameters> - The file with script params.
#		Contents:
#		[MutationNumbersFilter]="sites"|"pairs" - Apply filter of mutations' numbers on sites or site pairs. No default value.
#		[PairsOrdering]="0"|"1" - The type of site pairs: 0 - unordered or 1 - ordered
#		[BGR_SiteMinMutations]="<uint>" - Minimal number of mutations in background sites
#		[FGR_SiteMinMutations]="<uint>" - Minimal number of mutations in foreground sites
#		PairsType=("intra[gene]"|"inter[gene]")
#		SitePairs="FN" - list of site pairs for that the epistatic statistics were calculated
#		EpiStatistics="FN" - file with epistatic statistics for pairs
#		ObsOrderingStat="FN" - A file name with statistics for an observed dataset used for ranking of site pairs 
#		SortingOrder="<(asc)ending|(desc)ending>" - A sorditing order of statistics from the best to the worst
#		TestSiteSet="FN" - the set of sites for enrichment testing 
#		LinkSiteSetWith="<backgr[ound]|target|pairs>" Default="background". Search sites from the site set through the background or the target sites.
#		[PairsLinkLogic="<and|or>"] No defaul value. Used only with LinkSiteSetWith="pairs". 
#			Logic function that links a specified site set with a site pair: 'and' - with both sites, 'or' - at least one site in a pair.
#		[CDS="FN"] - The XML file with CDS statement.
#		[BackgroundSegmID]= "<1|2>" - Default=1. The CDS for background segment.
#		[LocusName="STR"] - The locus name containing sites from the testing site set.

use strict;
use lib "$ENV{EPISTAT_LIB}";
use lib "$ENV{EPISTAT_HOME}";
my $wilcox_r="$ENV{EPISTAT_HOME}/stat/wilcox_test.R";

use Getopt::Std;
use Class::Struct;
use DnaUtilities::CDS::annotation;
use SitePairMatrix;
use AssociationStatistics::SiteMutationNumberFilter;
use AssociationStatistics::PairMutationNumberFilter;
use AssociationStatistics::SelectedSitesFilter;

my %args;
if(!getopts('a:l:',\%args)){
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
	fg_filtered_sites_fn => '$'
};

my $wilcox_alternstive;
my $filters_params;
my $site_pair_mtx;
my $selected_sites_filter;

my $bgr_cds=1;
my $tgt_cds=1;

my ($pairs_fn,$stat_fn,$ordstat_fn,$sites_fn,$cds_fn,$out_csv_fn);
my $f_intragene;
my $link_idx;
my $link_logic;
my $locus_name;
my $sort_order;
my $pairs_ordering;
my @sort_order=("ASCENDING","DESCENDING");
my @logic_functs=("OR","AND");

if(defined $args{a}){
	if($args{a}=~/t(wo\.sided)?|l(ess)?|g(reater)?/i){
		$wilcox_alternstive=$args{a};
	}else{
		die "\nUnpropper value for the -a option!";
	}
}

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
		}elsif($key eq "SortingOrder"){
			if($value=~/^desc(ending)*/i){
				$sort_order=1;
			}elsif($value=~/^asc(ending)*/i){
				$sort_order=0;
			}
		}elsif($key eq "TestSiteSet"){
			$sites_fn=$value;
		}elsif($key eq "BackgroundSegmID"){
			$bgr_cds=$value;
			die "\nWrong CDS id for background segment was accounted: $bgr_cds! Expects <1|2>." unless($bgr_cds==1||$bgr_cds==2);
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
		}else{
			warn "\nUnknown parameter: $key!";
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
die "\nThe pairs' type was not defined intra- or intergenic!" unless defined $f_intragene;

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
if(defined $filters_params){
	$filters_params->obs_epistat_fn($stat_fn);
	$filters_params->site_pairs_fn($pairs_fn);
	$filters_params->f_intragene($f_intragene);
	$filters_params->pairs_ordering($pairs_ordering);
	die "\nUnable to initialize a site pair filter!" unless defined($filters_params->site_pairs_fn)&&defined($filters_params->f_intragene);
	my $infile_reading_mode=0;
	if(defined $filters_params->mnumber_filter_type){
		die "\nUndefined whether site pairs are ordered or not!" unless defined $filters_params->pairs_ordering;
		die "\nUnable to initialize a site pair filter!" unless $filters_params->fgr_site_min_muts>0||$filters_params->bgr_site_min_muts>0;
		die "\nUnordered site pairs could be defined only for intragenic pairs!" unless ($filters_params->pairs_ordering)||($filters_params->f_intragene);
		if($filters_params->mnumber_filter_type() eq 'sites'){
			if($filters_params->bgr_site_min_muts>0){
				$infile_reading_mode++ unless $infile_reading_mode==3;
			}
			if($filters_params->fgr_site_min_muts>0){
				$infile_reading_mode+=2 unless $infile_reading_mode>=2;
			}
		}elsif($filters_params->mnumber_filter_type() ne 'pairs'){
			die "\nUnknown value for the \"MutationNumbersFilter\" parameter!";
		}
	}elsif(defined($filters_params->fgr_site_min_muts)||defined($filters_params->bgr_site_min_muts)){
		die "\nUnable to initialize mutatin number filter: the \"MutationNumbersFilter\" parameter is undefined!"
	}
	$site_pair_mtx=SitePairMatrix->new($filters_params->site_pairs_fn,$filters_params->f_intragene,$infile_reading_mode);
	$site_pair_mtx->init_sites;
	if(defined($filters_params->bg_filtered_sites_fn)||defined($filters_params->fg_filtered_sites_fn)){
		$selected_sites_filter=AssociationStatistics::SelectedSitesFilter->new($site_pair_mtx,0,$filters_params->bg_filtered_sites_fn,$filters_params->fg_filtered_sites_fn);
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
#reading site set
my %site_subset;
open INFILE, "<$sites_fn" or die "\nUnable to open input file: $sites_fn!";
while(<INFILE>){
	$_=~s/\s//g;
	my @tmp=split ",";
	foreach my $str(@tmp){
		if($str=~/(\d+)-(\d+)/){
			for(my $i=$1;$i<=$2;$i++){
				$site_subset{$i}=1;
			}
		}elsif($str=~/(\d+)/){
			$site_subset{$1}=1;
		}
	}
}
close INFILE;

#information about detected pairs
struct SitePairInfo => {
	bg_site => '$',
	tg_site => '$',
	bg_nsubst => '$',
	tg_nsubst => '$',
	epistat => '$',
	npair_subst => '$',
	ordering_stats => '$'
};

struct SiteLocation =>{
	locus_name => '$',
	position => '$'
};

struct SiteInfo =>{
	location => '$',
	ordering_stats => '$',
	best_site_pair => '$'
};

my @site_pairs;
my %sites;

my @annotation;
if($cds_fn){
	#Reading the XML annotation file
	@annotation=DnaUtilities::CDS::annotation::read_annotation_xml($cds_fn);
	die "\nError: the subset: $sites_fn\n\t has linked with both background and foreground sites, thus common CDS for both is required!" if(@annotation>1&&$link_idx==2);
}

my $ra_site_pair_filter;
if(defined $filters_params){
	$ra_site_pair_filter=[];
	my $site_pair_filter;
	if($filters_params->mnumber_filter_type eq 'pairs'){
		$site_pair_filter=AssociationStatistics::PairMutationNumberFilter->new($site_pair_mtx,$filters_params->obs_epistat_fn,
			$filters_params->bgr_site_min_muts,$filters_params->fgr_site_min_muts);
	}elsif($filters_params->mnumber_filter_type eq 'sites'){
		$site_pair_filter=AssociationStatistics::SiteMutationNumberFilter->new($site_pair_mtx,
			$filters_params->bgr_site_min_muts,$filters_params->fgr_site_min_muts);
	}
	@{$ra_site_pair_filter}=(1) x $site_pair_mtx->{NLINES};
	@{$ra_site_pair_filter}=$site_pair_filter->apply($ra_site_pair_filter);
	@{$ra_site_pair_filter}=$selected_sites_filter->apply($ra_site_pair_filter) if defined $selected_sites_filter;
	symmetrize_filter($site_pair_mtx,$ra_site_pair_filter) unless $filters_params->pairs_ordering;
}
sub is_filter_ok{
	my ($ra_filter,$site_pair_idx)=@_;
	die "\nError is_filter_ok(): Undefined site pair index!" unless defined $site_pair_idx;
	if(defined $ra_filter){
		return $ra_filter->[$site_pair_idx];
	}
	return 1;
}

{
	open PAIRS, "<$pairs_fn" or die "\nUnable to open input file $pairs_fn!";
	<PAIRS>; #skip header line
	open STAT, "<$stat_fn" or die "\nUnable to open input file $stat_fn!";
	open ORDSTAT, "<$ordstat_fn" or die "\nUnable to open input file $ordstat_fn!";
	my $ntg=0;
	my $nbg=0;
	for(my $I=0;!(eof(PAIRS)||eof(STAT)||eof(ORDSTAT));$I++){
		my $str=<PAIRS>;
		my $spi=SitePairInfo->new();
		$spi->bg_site(0);
		$spi->tg_site(0);
		$spi->bg_nsubst(0);
		$spi->tg_nsubst(0);
		$spi->epistat(0);
		$spi->npair_subst(0);
		if($str=~/^(\d+)\t(\d+)\t(\d+)\t(\d+)/){
			$spi->bg_site($1);
			$spi->tg_site($2);
			$spi->bg_nsubst($3);
			$spi->tg_nsubst($4);
		};
		$str=<STAT>;
		if($str=~/^([\d.eE-]+)\t(\d+\.\d+)/){
			$spi->epistat($1);
			$spi->npair_subst($2);
		};
		$str=<ORDSTAT>;
		chomp $str;
		if($str=~/([-+]?\d*\.?\d+([eE]?[-+]?\d*))/){
			$spi->ordering_stats($1) if is_filter_ok($ra_site_pair_filter,$I);
		}
		
		next unless defined($spi->bg_site)&&defined($spi->tg_site);
		next unless is_filter_ok($ra_site_pair_filter,$I);
		if($link_idx==2){
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
					$site_info->ordering_stats($spi->ordering_stats);
					$site_info->best_site_pair($spi) if defined $spi->ordering_stats;
					$sites{$key_site}=$site_info;
				}else{
					my $bs=get_best_stat($sort_order,$sites{$key_site}->ordering_stats,$spi->ordering_stats);
					if(defined($bs)&&($bs!=$sites{$key_site}->ordering_stats)){
						$sites{$key_site}->ordering_stats($bs);
						$sites{$key_site}->best_site_pair($spi);
					}
				}
			}
		}
	};
	close PAIRS;
	close STAT;
	close ORDSTAT;
}

sub is_site_in_subset{
	my ($rsite_loc,$rh_site_set)=@_;
	if(!defined($rsite_loc->locus_name)|| $rsite_loc->locus_name eq $locus_name){
		return 1 if defined($rh_site_set->{$rsite_loc->position});
	};
	return 0;
}

sub get_sites_keys{
	my ($hr_sites)=@_;
	my @keys;
	foreach my $site(keys %{$hr_sites}){
		push @keys,$site if defined $hr_sites->{$site}->ordering_stats;
	}
	@keys=sort {$hr_sites->{$a}->ordering_stats <=> $hr_sites->{$b}->ordering_stats} @keys;
	return @keys;
}

sub get_pairs_keys{
	my ($ha_site_pairs)=@_;
	my %hkeys;
	my %pair_keys;
	for(my $i=0;$i<@{$ha_site_pairs};$i++){
		my $bgs=$ha_site_pairs->[$i]->[2]->bg_site;
		my $tgs=$ha_site_pairs->[$i]->[2]->tg_site;
		$pair_keys{$bgs.",".$tgs}=$i;
	}
	for(my $i=0;$i<@{$ha_site_pairs};$i++){
		my $val=$ha_site_pairs->[$i]->[2]->ordering_stats;
		if(defined $val){
			my $bgs=$ha_site_pairs->[$i]->[2]->bg_site;
			my $tgs=$ha_site_pairs->[$i]->[2]->tg_site;
			my $j=$pair_keys{$tgs.",".$bgs};
			if(defined $j){
				my $rval=$ha_site_pairs->[$j]->[2]->ordering_stats;
				my $best_val=get_best_stat($sort_order,$val,$rval);
				$hkeys{$val==$best_val?$i:$j}=1 unless defined($hkeys{$i})||defined($hkeys{$j});
			}else{
				$hkeys{$i}=1;
			}
		}
	}
	return sort {$ha_site_pairs->[$a]->[2]->ordering_stats <=> $ha_site_pairs->[$b]->[2]->ordering_stats} keys(%hkeys);
}

sub gen_tempname{
	my $nchar=shift;
	my @chars = ( "A" .. "Z", "a" .. "z", 0 .. 9 );
	return join("", @chars[ map { rand @chars } ( 1 .. $nchar ) ]);
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

my $tmp_fname=gen_tempname(10);
open OPF, ">$tmp_fname" or die "\nUnable to open output file: $tmp_fname!";
if($link_idx<2){
	print OPF "site\tstat\tin_subset";
	my @keys=get_sites_keys(\%sites);
	for(my $i=0;$i<@keys;$i++){
		my $stat=sprintf("%.4f",$sites{$keys[$i]}->ordering_stats);
		print OPF "\n".$keys[$i]."\t".$stat."\t";
		if(is_site_in_subset $sites{$keys[$i]}->location,\%site_subset){
			print OPF 1;
		}else{
			print OPF 0;
		}
		my $bspi=$sites{$keys[$i]}->best_site_pair;
		#if(defined $bspi){
		#	print OPF "\t".$bspi->bg_site."\t".$bspi->tg_site;#."\t".$bspi->ordering_stats(0);
		#}
	}
}else{
	print OPF "bg_site\tfg_site\tstat\tin_subset";
	my @keys=get_pairs_keys(\@site_pairs);
	foreach my $i(@keys){
		my $stat=sprintf("%.4f",$site_pairs[$i]->[2]->ordering_stats);
		my $bgs=$site_pairs[$i]->[2]->bg_site;
		my $tgs=$site_pairs[$i]->[2]->tg_site;
		print OPF "\n$bgs\t$tgs\t$stat\t";
		if(is_pair_in_subset($site_pairs[$i]->[0],$site_pairs[$i]->[1],\%site_subset,$link_logic)){
			print OPF 1;
		}else{
			print OPF 0;
		}
	}
}
close OPF;
my $str="$wilcox_r $tmp_fname ";
$str.="$wilcox_alternstive " if defined $wilcox_alternstive;
my $out_fn=$ARGV[0];
$out_fn=~s/\.[^.]+$//;
$out_fn.=".$logic_functs[$link_logic]" if $link_idx==2;
$out_fn.=".$wilcox_alternstive" if defined $wilcox_alternstive;
$out_fn.=".wilcox_test.out";
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

print STDERR "$\n$str";
system($str);
print_child_termination_status();
unlink $tmp_fname;
