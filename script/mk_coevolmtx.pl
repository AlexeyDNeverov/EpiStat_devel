#!/usr/bin/env perl
#The script generates matrix of coevolution measures for sites in one or two different proteins
#Usage: options <parameters>
#options:
#	-m <method> - Method used to calculate matrix values
#		<method>=Z-score|NZ-score|covariance-like|correlation-like
#	-b - Output bipartite graph matrix
#Params:
#	<parameters> - file with parameters
#		NumberOfDatasets="<uint>" - Number of datasets in the parameter file. Have to be declarated before dataset sections
#		BlockMtxName="STR" - The name of output block matrix
#			Default=undef - Print to STDOUT
#		BlockMtxExt="STR" - Extension of output file that contains block matrix
#			Default=undef - Print to STDOUT
#		SubMtxExt="STR" - Extension of output files that contain submatrices of a block matrix
#			Default=".sub.mtx";
#		OutDir="STR" - Folder for writing of output files
#			Default="."
#		NonNegativeElementFilter="0|1" - Set negative elements in the matrix to zeros
#			Default="0"
#		MutationNumbersFilter="sites"|"pairs" - Apply filter of mutations' numbers on sites or site pairs. No default value.
#		BGR_SiteMinMutations="<uint>" - Minimal number of mutations in background sites
#			Default="0"
#		FGR_SiteMinMutations="<uint>" - Minimal number of mutations in foreground sites
#			Default="0"
#		DistatnceType="SYM"|"ASYM"|"ALLELE" - The type of site pairs: symetric or asymetric
#			Default="ASYM"
#			ALLELE - the epistatic statistic is naturally symmetric and acconts for allele pairs. Files with statistics are upper triangular.
#		[Dataset] - Section statement
#			Name="STR" - Name of dataset
#			PairsType=("intra[gene]"|"inter[gene]")
#			EpiStat="FN" - A file with epistatic statistics
#			Mean="FN" - A file with mean epistatic statistics values
#			Pairs="FN" - A file with site pairs
#			BgrGeneName="STR" - An alias of background gene. Optional
#			FgrGeneName="STR" - An alias of foreground gene. Optional
#			Variance="FN" - A file with variances of epistatic statistics values. Optional
#			OrdPairsCovariance="FN" - A file with covariances of epistatic statistics of two orderes site pairs. Optional. Covariance is used only for symmetric intragenic matrices!
#			LowerPvalue="FN" - A file with pvalues. Optional
#			LowerPvalueThreshold="(0,1]"
#			UpperPvalue="FN" - A file with pvalues. Optional
#			UpperPvalueThreshold="(0,1]"
#			BgrSelectedSites="FN" - If the option is accounted, only sites in the background which are listed in the specified file are used for matrix construction. Optional.
#			FgrSelectedSites="FN" - If the option is accounted, only sites in the foreground which are listed in the specified file are used for matrix construction. Optional.
#			BothSitesSelected="0|1" Defines the way to apply selected site filter for intragenic site pairs. No default value.
#				!!!The parameter is required if PairsType="intragene"
#				0 - at least one site in a pair has to be in one of selected site subsets: 'BgrSelectedSites' or 'FgrSelectedSites'
#				1 - both sites in a pair have to be in one of selected site subsets: 'BgrSelectedSites' or 'FgrSelectedSites'
#		#End of description. Number of datasets have to be setup in the 'NumberOfDatasets' parameter

use strict;
use Getopt::Std;
use IO::File;
use Class::Struct;
use Cwd;

use lib "$ENV{EPISTAT_LIB}";

use AssociationStatistics::CorrelationLike;
use AssociationStatistics::CovarianceLike;
use AssociationStatistics::ZScore;
use AssociationStatistics::NormZScore;
use AssociationStatistics::SymCorrelationLike;
use AssociationStatistics::SymCovarianceLike;
use AssociationStatistics::SymZScore;
use AssociationStatistics::SymNormZScore;

struct DatasetDesc =>{
	name => '$',
	f_transposed => '$',
	f_intragene_pairs => '$',
	mean_fn => '$',
	variance_fn => '$',
	ord_pairs_covariance_fn => '$',
	pairs_fn => '$',
	epistat_fn => '$',
	lower_pvalue_fn => '$',
	upper_pvalue_fn => '$',
	bgr_site_moments12_fn => '$',
	fgr_site_moments12_fn => '$',
	lower_pvalue_threshold => '$',
	upper_pvalue_threshold => '$',
	bgr_gene_name => '$',
	fgr_gene_name => '$',
	bgr_selected_sites_fn => '$',
	fgr_selected_sites_fn => '$',
	f_both_sites_selected => '$'
};

struct GeneralSettings =>{
	no_negative_elements => '$',
	f_mutation_numbers => '$',
	bgr_min_mutations => '$',
	fgr_min_mutations => '$',
	distance_type => '$'
};

my %args;
if(!getopts('m:b',\%args)){
	die "\nError in option string!";
}

my @dataset_desc;
my $datasets_num;
my $method;
my $matrix_name;
my $matrix_ext;
my $outdir;
my $f_bipartite=0;
my $submatrix_ext=".sub.mtx";
my $general_settings=GeneralSettings->new;

die "\nNo default value for the option '-m' (method)!" unless defined $args{m};
if($args{m} eq "Z-score" || $args{m} eq "NZ-score" || $args{m} eq "covariance-like" || $args{m} eq "correlation-like"){
	$method=$args{m};
}else{
	die "\nError: Unknown method: $args{m}!";
}

$f_bipartite=1 if($args{b});
$general_settings->distance_type("asym");

#begin script
open INFILE, "$ARGV[0]" or die "\nUnable to open input file: $ARGV[0]!";
while(<INFILE>){
	$_=$` if(/#/);
	if(/\[\s*(.+)\s*\]/){
		my $name=$1;
		if($name eq "Dataset"){
			my $ds=DatasetDesc->new();
			$ds->f_transposed(0);
			push @dataset_desc,$ds;
		}else{
			die "\nUnknown section name: $name!";
		}
	}elsif(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
		my ($key,$value)=($1,$2);
		if($key eq "NumberOfDatasets"){
			$datasets_num=$value;
		}elsif($key eq "OutDir"){
			$outdir=$value;
			$outdir=~s/^\s+//;
			$outdir=~s/\s+$//;
			$outdir=Cwd::abs_path($outdir) if $outdir=~/\.\.?/;
			$outdir.="/" unless $outdir=~/\/$/;
		}elsif($key eq "BlockMtxName"){
			$matrix_name=$value;
		}elsif($key eq "BlockMtxExt"){
			$matrix_ext=$value;
		}elsif($key eq "SubMtxExt"){
			$submatrix_ext=$value;
		}elsif($key eq "NonNegativeElementFilter"){
			die "The parameter 'NonNegativeElementFilter' must take 0 or 1 values!" unless $value==0 || $value==1;
			$general_settings->no_negative_elements($value) if $value==1;
		}elsif($key eq "BGR_SiteMinMutations"){
			die "The parameter 'BGR_SiteMinMutations' must be UINT!" unless $value=~/^\d+$/;
			$general_settings->bgr_min_mutations($value);
		}elsif($key eq "FGR_SiteMinMutations"){
			die "The parameter 'FGR_SiteMinMutations' must be UINT!" unless $value=~/^\d+$/;
			$general_settings->fgr_min_mutations($value);
		}elsif($key eq "MutationNumbersFilter"){
			die "The parameter 'MutationNumbersFilter' must be 'sites' OR 'pairs'!" unless $value=~m/^sites|pairs$/i;
			$general_settings->f_mutation_numbers(lc($value));
		}elsif($key eq "DistatnceType"){
			die "The parameter 'DistatnceType' must be 'SYM','ASYM' or 'ALLELE'!" unless $value=~m/^sym|asym|allele/i;
			$general_settings->distance_type(lc($&));
			if(($general_settings->distance_type eq "sym")&&$f_bipartite){
				warn "\nBipartite graph reprezentation is inpossible for symmetric assosiations between sites!\nThe -b option will be ignored!";
				$f_bipartite=0;
			}
		}else{
			die "\nThe parameter 'NumberOfDatasets' must be set up befor datasets!" unless defined $datasets_num;
			die "\nThe dataset section hasn't been declarated!" unless defined $dataset_desc[-1];
			if($key eq "PairsType"){
				if($value=~/^intra(gene)*/){
					$dataset_desc[-1]->f_intragene_pairs(1);
				}elsif($value=~/^inter(gene)*/){
					$dataset_desc[-1]->f_intragene_pairs(0);
				}else{
					die "\nError: Unknown type of site pairs!";
				}
			}elsif($key eq "Mean"){
				$dataset_desc[-1]->mean_fn($value);
			}elsif($key eq "Pairs"){
				$dataset_desc[-1]->pairs_fn($value);
			}elsif($key eq "EpiStat"){
				$dataset_desc[-1]->epistat_fn($value);
			}elsif($key eq "Variance"){
				$dataset_desc[-1]->variance_fn($value);
			}elsif($key eq "OrdPairsCovariance"){
				$dataset_desc[-1]->ord_pairs_covariance_fn($value);
			}elsif($key eq "BgrGeneName"){
				$dataset_desc[-1]->bgr_gene_name($value);
			}elsif($key eq "FgrGeneName"){
				$dataset_desc[-1]->fgr_gene_name($value);
			}elsif($key eq "LowerPvalue"){
				$dataset_desc[-1]->lower_pvalue_fn($value);
			}elsif($key eq "UpperPvalue"){
				$dataset_desc[-1]->upper_pvalue_fn($value);
			}elsif($key eq "LowerPvalueThreshold"){
				$dataset_desc[-1]->lower_pvalue_threshold($value);
			}elsif($key eq "UpperPvalueThreshold"){
				$dataset_desc[-1]->upper_pvalue_threshold($value);
			}elsif($key eq "BgrSiteMarginStat"){
				$dataset_desc[-1]->bgr_site_moments12_fn($value);
			}elsif($key eq "FgrSiteMarginStat"){
				$dataset_desc[-1]->fgr_site_moments12_fn($value);
			}elsif($key eq "Name"){
				$dataset_desc[-1]->name($value);
			}elsif($key eq "BgrSelectedSites"){
				$dataset_desc[-1]->bgr_selected_sites_fn($value);
			}elsif($key eq "FgrSelectedSites"){
				$dataset_desc[-1]->fgr_selected_sites_fn($value);
			}elsif($key eq "BothSitesSelected"){
				die "\nOnly 0 or 1 values are allowed for the 'BothSitesSelected' parameter!" unless $value=~m/^[01]$/;
				$dataset_desc[-1]->f_both_sites_selected($value);
			}elsif($key eq "Tranformation"){
				if($value =~/^transp/i){
					$dataset_desc[-1]->f_transposed(1);
				}else{
					die "\nUnknown transformation type: $value!";
				}
			}else{
				die "\nUnknown parameter: $key in the input file $ARGV[0]!";
			}
		}
	}
}
close INFILE;
die "\nThe value of 'NumberOfDatasets' parameter is not equal to number of datasets!" unless @dataset_desc==$datasets_num;
die "\nUnable to define objects (sites or site pairs) for which the filter on mutations numbers have to be applied on" if 
	(!defined($general_settings->f_mutation_numbers))&&($general_settings->bgr_min_mutations>0||$general_settings->fgr_min_mutations>0);
foreach my $ds(@dataset_desc){
	if(defined($ds->bgr_selected_sites_fn)||defined($ds->fgr_selected_sites_fn)){
		die "\nIn the intragenic dataset ".$ds->name." the way how to select site pairs is not completely defined. Please specify 'BothSitesSelected'!" unless defined($ds->f_both_sites_selected)||(!$ds->f_intragene_pairs);
	}
}
sub print_submatrix{
	my ($out_fh,$f_transp,$name)=@_;
	if($f_transp){
		print $out_fh "\nt\t$name";
	}else{
		print $out_fh "\n\t$name";
	}
}

sub union_ofSiteSets{
	my ($hr_set1,$hr_set2)=@_;
	my %uset;
	foreach my $str(keys %{$hr_set1}){
		$uset{$str}=1;
	}
	foreach my $str(keys %{$hr_set2}){
		$uset{$str}=1;
	}
	return %uset;
}

sub get_matrix_norm_constant{
	my ($ra_datasets)=@_;
	my $M=0;
	foreach my $dts(@{$ra_datasets}){
		my @stat=$dts->get_statistics(0);
		foreach my $val(@stat){
			$M=abs($val) if abs($val)>$M;
		}
	}
	return $M;
}

sub mk_asym_matrix{
	my ($method,$descr,$settings)=@_;
	my $dts;
	if($method eq 'Z-score'){
		$dts=AssociationStatistics::ZScore->new($descr,$settings);
	}elsif($method eq 'NZ-score'){
		$dts=AssociationStatistics::NormZScore->new($descr,$settings);
	}elsif($method eq 'covariance-like'){
		$dts=AssociationStatistics::CovarianceLike->new($descr,$settings);
	}elsif($method eq 'correlation-like'){
		$dts=AssociationStatistics::CorrelationLike->new($descr,$settings);
	}
	return $dts;
}

#Make configuration of the covariance matrix
if($datasets_num>0&&$datasets_num<=4){
	my @intragenes;
	my @intergenes;
	my @datasets;
	my %sites;
	my %gene_counts;
	for(my $i=0;$i<$datasets_num;$i++){
		if($dataset_desc[$i]->f_intragene_pairs){
			push @intragenes,$i;
			die "\nTransposition of intragenic matrix is not allowed!" if $dataset_desc[$i]->f_transposed;
		}else{
			push @intergenes,$i;
		}
	}
	my $conf_ngenes;
	if(@intragenes>=0 && @intragenes<=2 && @intergenes==2 && 
		$dataset_desc[$intergenes[0]]->bgr_gene_name eq $dataset_desc[$intergenes[1]]->fgr_gene_name &&
		$dataset_desc[$intergenes[0]]->fgr_gene_name eq $dataset_desc[$intergenes[1]]->bgr_gene_name){
		$conf_ngenes=2;
		die "\nBoth intergenic matrices are transposed!" if $dataset_desc[$intergenes[0]]->f_transposed()&&$dataset_desc[$intergenes[1]]->f_transposed();
	}elsif(@intragenes==1 && @intergenes==0){
		$conf_ngenes=1;
	}else{
		die "\nError: Wrong configuration!";
	}
	my $dmtx=$conf_ngenes;
	$dmtx*=2 if $f_bipartite;
	if($general_settings->distance_type eq "asym"||$general_settings->distance_type eq "allele"){
		for(my $i=0;$i<$datasets_num;$i++){
			my $dts=mk_asym_matrix($method,$dataset_desc[$i],$general_settings);
			if($dataset_desc[$i]->f_transposed){
				push @datasets,$dts->transpose();
			}else{
				push @datasets,$dts;
			}
		}
	}else{
		if(@intergenes==2){
			my $dts;
			my @idx=($intergenes[0],$intergenes[1]);
			unless($dataset_desc[$intergenes[0]]->f_transposed()||$dataset_desc[$intergenes[1]]->f_transposed()){
				if($method eq 'Z-score'){
					$dts=AssociationStatistics::SymZScore->new($dataset_desc[$intergenes[0]],$dataset_desc[$intergenes[1]],$general_settings);
				}elsif($method eq 'NZ-score'){
					$dts=AssociationStatistics::SymNormZScore->new($dataset_desc[$intergenes[0]],$dataset_desc[$intergenes[1]],$general_settings);
				}elsif($method eq 'covariance-like'){
					$dts=AssociationStatistics::SymCovarianceLike->new($dataset_desc[$intergenes[0]],$dataset_desc[$intergenes[1]],$general_settings);
				}elsif($method eq 'correlation-like'){
					$dts=AssociationStatistics::SymCorrelationLike->new($dataset_desc[$intergenes[0]],$dataset_desc[$intergenes[1]],$general_settings);
				}
			}else{
				@idx=($intergenes[1],$intergenes[0]) unless $dataset_desc[$intergenes[1]]->f_transposed();
				$dts=mk_asym_matrix($method,$dataset_desc[$idx[0]],$general_settings);
			}
			if(defined $dts){
				$datasets[$idx[0]]=$dts;
				$datasets[$idx[1]]=$dts->transpose();
			}
		}
		for(my $i=0;$i<@intragenes;$i++){
			my $dts;
			my $idx=$intragenes[$i];
			if($method eq 'Z-score'){
				$dts=AssociationStatistics::SymZScore->new($dataset_desc[$idx],undef,$general_settings);
			}elsif($method eq 'NZ-score'){
				$dts=AssociationStatistics::SymNormZScore->new($dataset_desc[$idx],undef,$general_settings);
			}elsif($method eq 'covariance-like'){
				$dts=AssociationStatistics::SymCovarianceLike->new($dataset_desc[$idx],undef,$general_settings);
			}elsif($method eq 'correlation-like'){
				$dts=AssociationStatistics::SymCorrelationLike->new($dataset_desc[$idx],undef,$general_settings);
			}
			$datasets[$idx]=$dts if(defined $dts);
		}
	}
	for(my $i=0;$i<$datasets_num;$i++){
		my $dts=$datasets[$i];
		if(defined $dts){
			my %bgr_sites;
			my %fgr_sites;
			$dts->get_sites(\%bgr_sites,\%fgr_sites);
			%bgr_sites=union_ofSiteSets(\%bgr_sites,\%fgr_sites) if($dataset_desc[$i]->f_intragene_pairs);
			my $str=$dataset_desc[$i]->bgr_gene_name;
			$gene_counts{$str}++;
			if(!defined $sites{$str}){
				$sites{$str}={};
			}
			foreach my $site(keys %bgr_sites){
				$sites{$str}->{$site}++;
			}
			if(!$dataset_desc[$i]->f_intragene_pairs){
				$str=$dataset_desc[$i]->fgr_gene_name;
				$gene_counts{$str}++;
				if(!defined $sites{$str}){
					$sites{$str}={};
				}
				foreach my $site(keys %fgr_sites){
					$sites{$str}->{$site}++;
				}
			}
		}else{last;}
	}
	my %gene_nsites;
	foreach my $gene_name(keys %gene_counts){
		if($gene_counts{$gene_name}>1){
			foreach my $site(keys %{$sites{$gene_name}}){
				delete $sites{$gene_name}->{$site} unless $sites{$gene_name}->{$site} == $gene_counts{$gene_name};
			}
		}
		$gene_nsites{$gene_name}=scalar(keys %{$sites{$gene_name}});
	}
	my $out_fh=*STDOUT;
	if(defined $matrix_name){
		my $out_fn=$outdir.$matrix_name.$matrix_ext;
		$out_fh=IO::File->new(">$out_fn");
		die "\nUnable to open output file $out_fn" unless $out_fh;
	}
	print $out_fh "Block Matrix: $matrix_name";
	print $out_fh "[$dmtx,$dmtx]";
	print $out_fh "\nOutDir=\"$outdir\"";
	print $out_fh "\nSubMtxExt=\"$submatrix_ext\"";
	if($conf_ngenes==2){
		#my $norm_const=0; #no normalization
		my $norm_const; #default normalization
		#if(($method eq 'Z-score')||($method eq 'correlation-like')){
		#	$norm_const=get_matrix_norm_constant(\@datasets);
		#}
		my @submatrices;
		my $gene1_name=$dataset_desc[$intergenes[0]]->bgr_gene_name;
		my $n1=$gene_nsites{$gene1_name};
		my $gene2_name=$dataset_desc[$intergenes[0]]->fgr_gene_name;
		my $n2=$gene_nsites{$gene2_name};
		print $out_fh "\nGenes: ($gene1_name";
		if($f_bipartite){
			print $out_fh "[2 x $n1],";
		}else{
			print $out_fh "[$n1],";
		}
		print $out_fh $gene2_name;
		if($f_bipartite){
			print $out_fh "[2 x $n2])";
		}else{
			print $out_fh "[$n2])";
		}
		print $out_fh "\nSubmatrices:\noperator\tname\tnrow\tncol";
		#row 1
		my $i=0;
		for(;$i<@intragenes;$i++){
			last if $dataset_desc[$intragenes[$i]]->bgr_gene_name eq $gene1_name;
		}
		my $dts;
		if($i<@intragenes){
			$dts=$datasets[$intragenes[$i]];
			$_=$dts->get_name;
			$submatrices[0]=[("$_\t$n1\t$n1",$dts)];
		}else{
			$submatrices[0]=[("\t$n1\t$n1",undef)];
		}
		$dts=$datasets[$intergenes[0]];
		$_=$dts->get_name;
		$submatrices[1]=[("$_\t$n1\t$n2",$dts)];
		#row 2
		$dts=$datasets[$intergenes[1]];
		$_=$dts->get_name;
		$submatrices[2]=[("$_\t$n2\t$n1",$dts)];
		$i=0;
		for(;$i<@intragenes;$i++){
			last if $dataset_desc[$intragenes[$i]]->bgr_gene_name eq $gene2_name;
		}
		if($i<@intragenes){
			$dts=$datasets[$intragenes[$i]];
			$_=$dts->get_name;
			$submatrices[3]=[("$_\t$n2\t$n2",$dts)]; 
		}else{
			$submatrices[3]=[("\t$n2\t$n2",undef)];
		}
		if($f_bipartite){
			#print row 1
			print_submatrix($out_fh,0,"0\t$n1\t$n1");
			my $dts=$submatrices[0]->[1];
			if(defined $dts){
				$_=$dts->get_name;
				print_submatrix($out_fh,0,$submatrices[0]->[0]);
				open OPF, ">$outdir$_$submatrix_ext" or die "\nUnable to open output file $outdir$_$submatrix_ext!";
				$dts->print_matrix(*OPF,$sites{$gene1_name},$sites{$gene1_name},'-norm'=>$norm_const);
				close OPF; 
			}else{
				print_submatrix($out_fh,0,"0".$submatrices[0]->[0]);
			}
			print_submatrix($out_fh,0,"0\t$n1\t$n2");
			$dts=$submatrices[1]->[1];
			$_=$dts->get_name;
			print_submatrix($out_fh,0,$submatrices[1]->[0]);
			open OPF, ">$outdir$_$submatrix_ext" or die "\nUnable to open output file $outdir$_$submatrix_ext!";
			$dts->print_matrix(*OPF,$sites{$gene1_name},$sites{$gene2_name},'-norm'=>$norm_const);
			close OPF;
			#print row 2
			if(defined $submatrices[0]->[1]){
				print_submatrix($out_fh,1,$submatrices[0]->[0]);
			}else{
				print_submatrix($out_fh,0,"0".$submatrices[0]->[0]);
			}
			print_submatrix($out_fh,0,"0\t$n1\t$n1");
			print_submatrix($out_fh,1,$submatrices[2]->[0]);
			print_submatrix($out_fh,0,"0\t$n1\t$n2");
			#print row 3
			print_submatrix($out_fh,0,"0\t$n2\t$n1");
			$dts=$submatrices[2]->[1];
			$_=$dts->get_name;
			print_submatrix($out_fh,0,$submatrices[2]->[0]);
			open OPF, ">$outdir$_$submatrix_ext" or die "\nUnable to open output file $outdir$_$submatrix_ext!";
			$dts->print_matrix(*OPF,$sites{$gene2_name},$sites{$gene1_name},'-norm'=>$norm_const);
			close OPF;
			print_submatrix($out_fh,0,"0\t$n2\t$n2");
			$dts=$submatrices[3]->[1];
			if(defined $dts){
				$_=$dts->get_name;
				print_submatrix($out_fh,0,$submatrices[3]->[0]);
				open OPF, ">$outdir$_$submatrix_ext" or die "\nUnable to open output file $outdir$_$submatrix_ext!";
				$dts->print_matrix(*OPF,$sites{$gene2_name},$sites{$gene2_name},'-norm'=>$norm_const);
				close OPF; 
			}else{
				print_submatrix($out_fh,0,"0".$submatrices[3]->[0]);
			}
			#print row 4
			print_submatrix($out_fh,1,$submatrices[1]->[0]);
			print_submatrix($out_fh,0,"0\t$n2\t$n1");
			if(defined $submatrices[3]->[1]){
				print_submatrix($out_fh,1,$submatrices[3]->[0]);
			}else{
				print_submatrix($out_fh,0,"0".$submatrices[3]->[0]);
			}
			print_submatrix($out_fh,0,"0\t$n2\t$n2");
		}else{
			#print row 1
			my $dts=$submatrices[0]->[1];
			if(defined $dts){
				$_=$dts->get_name;
				print_submatrix($out_fh,0,$submatrices[0]->[0]);
				open OPF, ">$outdir$_$submatrix_ext" or die "\nUnable to open output file $outdir$_$submatrix_ext!";
				$dts->print_matrix(*OPF,$sites{$gene1_name},$sites{$gene1_name},'-norm'=>$norm_const);
				close OPF; 
			}else{
				print_submatrix($out_fh,0,"I".$submatrices[0]->[0]);
			}
			$dts=$submatrices[1]->[1];
			$_=$dts->get_name;
			print_submatrix($out_fh,0,$submatrices[1]->[0]);
			open OPF, ">$outdir$_$submatrix_ext" or die "\nUnable to open output file $outdir$_$submatrix_ext!";
			$dts->print_matrix(*OPF,$sites{$gene1_name},$sites{$gene2_name},'-norm'=>$norm_const);
			close OPF;
			#print row 2
			$dts=$submatrices[2]->[1];
			$_=$dts->get_name;
			print_submatrix($out_fh,0,$submatrices[2]->[0]);
			open OPF, ">$outdir$_$submatrix_ext" or die "\nUnable to open output file $outdir$_$submatrix_ext!";
			$dts->print_matrix(*OPF,$sites{$gene2_name},$sites{$gene1_name},'-norm'=>$norm_const);
			close OPF;
			$dts=$submatrices[3]->[1];
			if(defined $dts){
				$_=$dts->get_name;
				print_submatrix($out_fh,0,$submatrices[3]->[0]);
				open OPF, ">$outdir$_$submatrix_ext" or die "\nUnable to open output file $outdir$_$submatrix_ext!";
				$dts->print_matrix(*OPF,$sites{$gene2_name},$sites{$gene2_name},'-norm'=>$norm_const);
				close OPF; 
			}else{
				print_submatrix($out_fh,0,"I".$submatrices[3]->[0]);
			}
		}
	}elsif($conf_ngenes==1){
		my $gene1_name=$dataset_desc[$intragenes[0]]->bgr_gene_name;
		my $n1=$gene_nsites{$gene1_name};
		my $dts=$datasets[$intragenes[0]];
		$_=$dts->get_name;
		print $out_fh "\nGenes: ($gene1_name";
		if($f_bipartite){
			print $out_fh "[2 x $n1])";
		}else{
			print $out_fh "[$n1])";
		}
		print $out_fh "\nSubmatrices:\noperator\tname\tnrow\tncol";
		if($f_bipartite){
			#print row 1
			print_submatrix($out_fh,0,"0\t$n1\t$n1");
			print_submatrix($out_fh,0,$_."\t$n1\t$n1");
			#print row 2
			print_submatrix($out_fh,1,$_."\t$n1\t$n1");
			print_submatrix($out_fh,0,"0\t$n1\t$n1");
		}else{
			print_submatrix($out_fh,0,$_."\t$n1\t$n1");
		}
		open OPF, ">$outdir$_$submatrix_ext" or die "\nUnable to open output file $outdir$_$submatrix_ext!";
		$dts->print_matrix(*OPF,$sites{$gene1_name},$sites{$gene1_name});
		close OPF;
	}
	undef $out_fh if(defined $out_fh);
}
	
	