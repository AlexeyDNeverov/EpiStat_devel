#!/usr/bin/env perl
#This script calculates Q-values for pvalues of epistatic statistics
#Usage: <parameters>
#<parameters> - configuration file
#		Contents:
#		[MutationNumbersFilter]="sites"|"pairs" - Apply filter of mutations' numbers on sites or site pairs. No default value.
#		[PairsOrdering]="0"|"1" - The type of site pairs: 0 - unordered or 1 - ordered
#			Default="1"
#		[PairsType]=("intra[gene]"|"inter[gene]")
#		[Pairs]="FN" - A file with site pairs
#		[Epistat_Ext="STR"] - An extation of files with epistatic statistics
#		[ObsEpistat="FN"] - A file name with observed epistatic statistics
#		[BGR_SiteMinMutations]="<uint>" - Minimal number of mutations in background sites
#		[FGR_SiteMinMutations]="<uint>" - Minimal number of mutations in foreground sites
#		[BgrSelectedSites="FN"] - If the option is accounted, only sites in the background which are listed in the specified file are used for FDR estimation.
#		[FgrSelectedSites="FN"] - If the option is accounted, only sites in the foreground which are listed in the specified file are used for FDR estimation.
#		[BothSitesSelected="0|1"] Defines the way to apply selected site filter for intragenic site pairs. No default value.
#				!!!The parameter is required if PairsType="intragene"
#				0 - at least one site in a pair has to be in one of selected site subsets: 'BgrSelectedSites' or 'FgrSelectedSites'
#				1 - both sites in a pair have to be in one of selected site subsets: 'BgrSelectedSites' or 'FgrSelectedSites'
#		FDRSampleNum="<uint>" - Number of samples used to calculate FDR
#		ObsPvalues="FN" - A file with pvalues
#		FDR_Dir="STR" - A path to a directory with random samples for a FDR estimate
#		FDRPvalue_Ext="STR" - An extaition of files storing pvalues on fakes


use strict;
use Class::Struct;
use File::Basename;

use lib "$ENV{EPISTAT_LIB}";

use SitePairMatrix;
use AssociationStatistics::SiteMutationNumberFilter;
use AssociationStatistics::PairMutationNumberFilter;
use AssociationStatistics::SelectedSitesFilter;


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
my ($pvalues_fn,$fdr_dir,$nfdr_samples,$pvalue_ext,$fdr_basename_suffix);
my $selected_sites_filter;
open INFILE, "$ARGV[0]" or die "\nUnable to open input file: $ARGV[0]!";
while(<INFILE>){
	$_=$` if(/#/);
	if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
		my ($key,$value)=($1,$2);
		if($key eq "FDRSampleNum"){
			die "\nUINT value for the parameter 'FDRSampleNum' is expected!" unless $value=~/^\d+$/;
			$nfdr_samples=$value;
		}elsif($key eq "ObsPvalues"){
			$pvalues_fn=$value;
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
		}elsif($key eq "Pairs"){
			$filters_params=SiteFilterParams->new unless defined $filters_params;
			$filters_params->site_pairs_fn($value);
		}elsif($key eq "Epistat_Ext"){
			$filters_params=SiteFilterParams->new unless defined $filters_params;
			$filters_params->epistat_ext($value);
		}elsif($key eq "FDRPvalue_Ext"){
			$value=~s/\s+$//;
			if($value=~/\..*$/){
				if($-[0]<length($value)){
					$fdr_basename_suffix=substr $value,0,$-[0];
					$value=substr $value,$-[0],length($value)-$-[0];
				}
			}
			$pvalue_ext=$value;
		}elsif($key eq "FDR_Dir"){
			$value.="/" unless $value=~/\/$/;
			$fdr_dir=$value;
		}elsif($key eq "MutationNumbersFilter"){
			die "The parameter 'MutationNumbersFilter' must be 'sites' OR 'pairs'!" unless $value=~m/^sites|pairs$/i;
			$filters_params=SiteFilterParams->new unless defined $filters_params;
			$filters_params->mnumber_filter_type(lc($value));
		}elsif($key eq "ObsEpistat"){
			$filters_params=SiteFilterParams->new unless defined $filters_params;
			$filters_params->obs_epistat_fn($value);
		}elsif($key eq "BGR_SiteMinMutations"){
			$filters_params=SiteFilterParams->new unless defined $filters_params;
			die "\nUINT value for the parameter 'BGR_SiteMinMutations' is expected!" unless $value=~/^\d+$/;
			$filters_params->bgr_site_min_muts($value);
		}elsif($key eq "FGR_SiteMinMutations"){
			$filters_params=SiteFilterParams->new unless defined $filters_params;
			die "\nUINT value for the parameter 'FGR_SiteMinMutations' is expected!" unless $value=~/^\d+$/;
			$filters_params->fgr_site_min_muts($value);
		}elsif($key eq "PairsType"){
			$filters_params=SiteFilterParams->new unless defined $filters_params;
			if($value=~m/^intra(gene)*/i){
				$filters_params->f_intragene(1);
			}elsif($value=~m/^inter(gene)*/i){
				$filters_params->f_intragene(0);
			}else{
				die "\nUnknown value ($value) for the PairsType parameter in file $ARGV[0]!"; 
			}
		}elsif($key eq "PairsOrdering"){
			$filters_params=SiteFilterParams->new unless defined $filters_params;
			if($value==0||$value==1){
				$filters_params->pairs_ordering($value);
			}else{
				die "\nUnknown value ($value) for the 'PairsOrdering' parameter in file $ARGV[0]!"; 
			}
		}else{
			die "\nUnknown parameter in $ARGV[0]: '$key'!";
		}
	}
}
close INFILE;

die "\nThe 'FDR_Dir' parameter hasn't been set!" unless defined $fdr_dir;
die "\nThe 'Pvalue_Ext' parameter hasn't been set!" unless defined $pvalue_ext;
die "\nThe 'FDRSampleNum' parameter hasn't been set!" unless defined $nfdr_samples;
if(defined $filters_params){
	die "\nUnable to initialize a site pair filter!" unless defined($filters_params->site_pairs_fn)&&defined($filters_params->f_intragene);
	my $infile_reading_mode=0;
	$filters_params->pairs_ordering(1) unless defined $filters_params->pairs_ordering;
	if(defined $filters_params->mnumber_filter_type){
		die "\nUnable to initialize a site pair filter!" unless $filters_params->fgr_site_min_muts>0||$filters_params->bgr_site_min_muts>0;
		die "\nUnordered pairs could be defined only for intragenic pairs!" unless ($filters_params->pairs_ordering)||($filters_params->f_intragene);
		if($filters_params->mnumber_filter_type() eq 'pairs'){
			die "\nUndefined required parameter \"ObsEpistat\" for initialization of the 'pairs' filter!" unless defined $filters_params->obs_epistat_fn;
			die "\nSuffix of fake files with epistatic statistics is not defined!" unless defined $filters_params->epistat_ext;
		}elsif($filters_params->mnumber_filter_type() eq 'sites'){
			if($filters_params->bgr_site_min_muts>0){
				$infile_reading_mode++ unless $infile_reading_mode==3;
			}
			if($filters_params->fgr_site_min_muts>0){
				$infile_reading_mode+=2 unless $infile_reading_mode>=2;
			}
		}else{
			die "\nUnknown \"MutationNumbersFilter\" parameter's value!";
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


sub read_statistics{
	my ($stats_fn,$ra_stat)=@_;
	open INPF, "<$stats_fn" or die "\nUnable to open input file: $stats_fn!";
	@{$ra_stat}=();
	while(<INPF>){
		chomp;
		$_=~s/\s*$//;
		my @line=split "\t";
		if($line[0] ne ""){
			push @{$ra_stat},$line[0];
		}
	}
	close INPF;
}

sub calc_fake_pvalue_counts{
	my ($ra_obs_pvals,$ra_fake_pvalues,$ra_out_counts,$ra_site_pair_filter)=@_;
	my $n=@{$ra_fake_pvalues};
	if(defined $ra_site_pair_filter){
		die "\nWrong number of elements in input arrays!" unless @{$ra_site_pair_filter}==$n;
	}
	my @obs_pvalues=sort {$a<=>$b} @{$ra_obs_pvals};
	my @fake_pvalues;
	if(!defined $ra_site_pair_filter){
		@fake_pvalues=sort {$a<=>$b} @{$ra_fake_pvalues};
	}else{
		for(my $i=0;$i<$n;$i++){
			push @fake_pvalues, $ra_fake_pvalues->[$i] if($ra_site_pair_filter->[$i]);
		}
		@fake_pvalues=sort {$a<=>$b} @fake_pvalues;
	}
	my @counts;
	#my @pvalues;
	my $j=0;
	my $i;
	for($i=0;$i<@obs_pvalues;$i++){
		my $pval=$obs_pvalues[$i];
		next unless defined $pval;
		#last if $pval == 1;
#print "\n$i\t$pval\t$j\t$fake_pvalues[$j]\t$fake_pvalues[$j-1]";
		if($i>0 && $obs_pvalues[$i-1]<$pval){
			for(;$j<@fake_pvalues;$j++){
				last if $fake_pvalues[$j]>$obs_pvalues[$i-1];
			}
			#push @pvalues,$obs_pvalues[$i-1];
			push @counts, $j;
		}
	}
	if($i>0){
		for(;$j<@fake_pvalues;$j++){
			last if $fake_pvalues[$j]>$obs_pvalues[$i-1];
		}
		#push @pvalues,$obs_pvalues[$i-1];
		push @counts, $j;
	}
	@{$ra_out_counts}=@counts;
	#@{$ra_out_pvalues}=@pvalues;
}

sub calc_obs_pvalue_counts{
	my ($ra_obs_pvalues,$ra_out_counts,$ra_out_pvals,$ra_site_pair_filter)=@_;
	my $n=@{$ra_obs_pvalues};
	if(defined $ra_site_pair_filter){
		die "\nWrong number of elements in input arrays!" unless @{$ra_site_pair_filter}==$n;
	}
	my @obs_pvalues;
	if(!defined $ra_site_pair_filter){
		@obs_pvalues=sort {$a<=>$b} @{$ra_obs_pvalues};
	}else{
		for(my $i=0;$i<$n;$i++){
			push @obs_pvalues, $ra_obs_pvalues->[$i] if($ra_site_pair_filter->[$i]);
		}
		@obs_pvalues=sort {$a<=>$b} @obs_pvalues;
	}
	my @counts;
	my @pvalues;
	my $i;
	for($i=0;$i<@obs_pvalues;$i++){
		my $pval=$obs_pvalues[$i];
		next unless defined $pval;
		#last if $pval == 1;
		if($i>0 && $obs_pvalues[$i-1]<$pval){
			push @pvalues,$obs_pvalues[$i-1];
			push @counts, $i;
		}
	}
	if($i>0){
		push @pvalues,$obs_pvalues[$i-1];
		push @counts, $i;
	}
	@{$ra_out_counts}=@counts;
	@{$ra_out_pvals}=@pvalues;
}

sub symmetrize_filter{
	my ($sp_matrix,$ra_filter)=@_;
	for(my $i=0;$i<@{$ra_filter};$i++){
		my ($bgs,$tgs)=$sp_matrix->line2sites_idxs($i);
		my $isym=$sp_matrix->site_idx_pair2line($tgs,$bgs);
		$ra_filter->[$i]=1 if($ra_filter->[$i]||$ra_filter->[$isym]);
	}
}

#begin script
my @obs_pvalues;
my @obs_pvals;
my @obs_counts;
my $site_pair_filter;
my $ra_site_pair_filter;
read_statistics($pvalues_fn,\@obs_pvalues);
if(defined $filters_params){
	if($filters_params->mnumber_filter_type eq 'pairs'){
		$site_pair_filter=AssociationStatistics::PairMutationNumberFilter->new($site_pair_mtx,$filters_params->obs_epistat_fn,
			$filters_params->bgr_site_min_muts,$filters_params->fgr_site_min_muts);
	}elsif($filters_params->mnumber_filter_type eq 'sites'){
		$site_pair_filter=AssociationStatistics::SiteMutationNumberFilter->new($site_pair_mtx,
			$filters_params->bgr_site_min_muts,$filters_params->fgr_site_min_muts);
	}
	$ra_site_pair_filter=[];
	@{$ra_site_pair_filter}=(1) x @obs_pvalues;
	@{$ra_site_pair_filter}=$site_pair_filter->apply($ra_site_pair_filter);
	@{$ra_site_pair_filter}=$selected_sites_filter->apply($ra_site_pair_filter) if defined $selected_sites_filter;
	symmetrize_filter($site_pair_mtx,$ra_site_pair_filter) unless $filters_params->pairs_ordering;
}
calc_obs_pvalue_counts(\@obs_pvalues,\@obs_counts,\@obs_pvals,$ra_site_pair_filter);
opendir(DIR, $fdr_dir) or die "can't opendir $fdr_dir: $!";
my $file;
my @fdr_samples;
my $nfdr_stat=0;
my $n;
$n=length($fdr_basename_suffix) if defined $fdr_basename_suffix;
while (defined($file = readdir(DIR))) {
	my ($basename,$dir,$ext) = fileparse($file,'\..*$');
	if($ext eq $pvalue_ext){
		$basename=substr $basename,0,-$n if(defined $fdr_basename_suffix);
		push @fdr_samples,$basename;
	}elsif(defined($filters_params)&&$filters_params->mnumber_filter_type eq 'pairs'){
		$nfdr_stat++ if $ext eq $filters_params->epistat_ext;
	}
}
closedir(DIR);
die "\nNumber of '*$pvalue_ext' files in $fdr_dir is not equal to FDRSampleNum=$nfdr_samples!" unless $nfdr_samples==@fdr_samples;
if(defined($filters_params)&&$filters_params->mnumber_filter_type eq 'pairs'){
	my $str=$filters_params->epistat_ext;
	die "\nNumber of '*$str' files in $fdr_dir is not equal to FDRSampleNum=$nfdr_samples!" unless $nfdr_samples==$nfdr_stat;
}
my @avg_fake_discovs;
my @upper_fake_counts;
for(my $i=0;$i<$nfdr_samples;$i++){
	my $fake_fn=$fdr_dir.$fdr_samples[$i];
	my @fake_pvalues;
	read_statistics($fake_fn.$fdr_basename_suffix.$pvalue_ext,\@fake_pvalues);
	my @fake_counts;
	if(defined($filters_params)&&$filters_params->mnumber_filter_type eq 'pairs'){
		$site_pair_filter=AssociationStatistics::PairMutationNumberFilter->new($site_pair_mtx,$fake_fn.$filters_params->epistat_ext,
			$filters_params->bgr_site_min_muts,$filters_params->fgr_site_min_muts);
		@{$ra_site_pair_filter}=(1) x @obs_pvalues;
		@{$ra_site_pair_filter}=$site_pair_filter->apply($ra_site_pair_filter);
		@{$ra_site_pair_filter}=$selected_sites_filter->apply($ra_site_pair_filter) if defined $selected_sites_filter;
		symmetrize_filter($site_pair_mtx,$ra_site_pair_filter) unless $filters_params->pairs_ordering;
	}
	calc_fake_pvalue_counts(\@obs_pvals,\@fake_pvalues,\@fake_counts,$ra_site_pair_filter);
	die "\nUnexpected error!" unless @fake_counts==@obs_counts;
	for(my $k=0;$k<@obs_counts;$k++){
		$avg_fake_discovs[$k]+=$fake_counts[$k];
		$upper_fake_counts[$k]+=1 if $fake_counts[$k]>=$obs_counts[$k];
	}
}
for(my $j=0;$j<@obs_counts;$j++){
	$avg_fake_discovs[$j]/=$nfdr_samples;
	$upper_fake_counts[$j]/=$nfdr_samples;
}
print "pvalue\t#obs\t#exp\tP";
for(my $j=0;$j<@obs_pvals;$j++){
	my $pval=$obs_pvals[$j];
	#last if $pval==1;
	print "\n$pval";
	my $str=$obs_counts[$j];
	print "\t$str";
	$str=sprintf("%.1f", $avg_fake_discovs[$j]);
	print "\t$str";
	$str=sprintf("%.3f", $upper_fake_counts[$j]);
	print "\t$str";
}