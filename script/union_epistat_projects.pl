#!/usr/bin/env perl
#This script concatenates files produced by 'run_epistat.pl' for specified projects
#The script generates in the specified folder for each type of results files a corresponding file with vertical concatenation
#	of corresponding files from projects folders
#Usage: options Params
#options:
#	-p <FN> - Specifies a file with a list of project folders for union
#	[-b] <FN> - Filter background sites. Expects a file in CSV format with a list of selected site IDs
#	[-t] <FN> - Filter foreground (target) sites. Expects a file in CSV format with a list of selected site IDs
#	[-o] <DIR> - Specifies an output folder. Default='.'
#Params:
#<epistat_prm> - epistat.pl configuration file
#<run_epistat_prm> - configuration file
use strict;
use Getopt::Std;
use File::Basename;

my $prj_list_fn;
my $out_dir="./";
my $bgr_filter_fn;
my $fgr_filter_fn;

my %args;
if(!getopts('p:b:t:o:',\%args)){
	die "\nError in option string!";
}
if(defined $args{o}){
	$out_dir=$args{o};
	$out_dir.="/" unless $out_dir=~/\/$/;
}
if(defined $args{p}){
	$prj_list_fn=$args{p};
}else{
	die "\nThe list with projects folders is not specified!";
}
$bgr_filter_fn=$args{b} if defined $args{b};
$fgr_filter_fn=$args{t} if defined $args{t};

my $pairs_ext=".site_pairs";
my $f_intragene;
my $f_stat_type;
my $epistat_prm_ext;
{
	my ($basename,$dir,$ext) = fileparse($ARGV[0],'\..*?$');
	$epistat_prm_ext=$ext;
}
open INFILE, "$ARGV[0]" or die "\nUnable to open input file: $ARGV[0]!";
while(<INFILE>){
	$_=$` if(/#/);
	if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
		my ($key,$value)=($1,$2);
		if($key eq "OutSitePairs"){
			$pairs_ext=$value;
		}elsif($key eq "PairsType"){
			if($value=~m/^intra(gene)*/i){
				$f_intragene=1;
			}elsif($value=~m/^inter(gene)*/i){
				$f_intragene=0;
			}else{
				die "\nUnknown value ($value) for the PairsType parameter in file $ARGV[0]!"; 
			}
		}elsif($key eq "StatType"){
			if($value=~m/^all_pairs/i){
				$f_stat_type=0;
			}elsif($value=~m/^independent_pairs/i){
				$f_stat_type=1;
			}elsif($value=~m/^branch_pairs/i){
				$f_stat_type=2;
			}elsif($value=~m/^all_allele_pairs/i){
				$f_stat_type=3;
			}elsif($value=~m/^indep_allele_pairs/i){
				$f_stat_type=3;
			}elsif($value=~m/^allele_branch_pairs/i){
				$f_stat_type=3;
			}else{
				die "\nUndefined type of statistics: StatType=$value!";
			}
		}
	}
}
close INFILE;
die "\nError: Parameter PairsType is undefined in the file $ARGV[0]!" unless defined $f_intragene;

sub get_epistat_ext{
	my ($epistat_prm)=@_;
	my ($stats_ext,$tau);
	open INFILE, "<$epistat_prm" or die "\nUnable to open input file: $epistat_prm!";
	while(<INFILE>){
		$_=$` if(/#/);
		if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
			my ($key,$value)=($1,$2);
			if($key eq "TAU"){
				if($value=~/\d+(\.\d+)*/){
					$tau=$value;
				}else{
					die "\nWrong value for the TAU parameter in file $epistat_prm!\n\tPositive number expected!";
				}
			}elsif($key eq "OutEpiStatistics"){
				$stats_ext=$value;
			}
		}
	}
	close INFILE;
	die "\nError get_epistat_ext(): Unable to define the epistat file extension!" unless defined($stats_ext)||defined($tau);
	unless(defined $stats_ext){
		$stats_ext=".stat.exp";
		$stats_ext.=".$tau" if defined $tau;
	}
	return $stats_ext;
}

my ($upper_pval_ext,$lower_pval_ext,$avg_ext,$var_ext);
my $site_stat_ext;
my $f_splited_bgtg;
my $cov_ext=".ord_site_pars.cov";
open INFILE, "$ARGV[1]" or die "\nUnable to open input file: $ARGV[1]!";
while(<INFILE>){
	$_=$` if(/#/);
	if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
		my ($key,$value)=($1,$2);
		if($key eq "OutUpperPValues"){
			$upper_pval_ext=$value;
		}elsif($key eq "OutLowerPValues"){
			$lower_pval_ext=$value;
		}elsif($key eq "OutAvgEpiStats"){
			$avg_ext=$value;
		}elsif($key eq "OutVarEpiStat"){
			$var_ext=$value;
		}elsif($key eq "OutCovOrdPairsEpiStat"){
			$cov_ext=$value;
		}elsif($key eq "OutSiteStats"){
			$site_stat_ext=$value;
		}elsif($key eq "SplitBgrAndFgr"){
			$value=~s/\s+$//;
			$value=~s/^\s+//;
			$f_splited_bgtg=$value if $value=~/^0|1$/;
		}
	}
}
close INFILE;
if(defined $f_splited_bgtg){
	unless($f_intragene){
		die "\nThe parameter 'SplitBgrAndFgr' could not be set up to the 0 value!" if $f_splited_bgtg==0;
	}
}else{
	if($f_intragene){
		$f_splited_bgtg=0;
	}else{
		$f_splited_bgtg=1;
	}
}
die "\nError: The OutUpperPValues parameter has no defaul value!" unless defined $upper_pval_ext;
die "\nError: The OutLowerPValues parameter has no defaul value!" unless defined $lower_pval_ext;
die "\nError: The OutAvgEpiStats parameter has no defaul value!" unless defined $avg_ext;
die "\nError: The OutVarEpiStat parameter has no defaul value!" unless defined $var_ext;
die "\nError: The OutBgrSiteStats parameter has no defaul value!" unless defined $site_stat_ext;

my @out_file_types;
my $pval_default_string="1";
{
	my $pvalue_mode_str="";
	if($f_intragene){
		$pvalue_mode_str=".intragene.unord_pairs";
		push @out_file_types,$pvalue_mode_str.$lower_pval_ext;
		push @out_file_types,$pvalue_mode_str.$upper_pval_ext;
		if($f_stat_type==3){
			$pvalue_mode_str="";
			push @out_file_types,$site_stat_ext;
		}else{
			$pvalue_mode_str=".intragene.ord_pairs";
			push @out_file_types,$cov_ext;
		}
	}
	push @out_file_types,$avg_ext;
	push @out_file_types,$var_ext;
	unless(($f_stat_type==3)&&$f_intragene&&($f_splited_bgtg==0)){
		push @out_file_types,$pvalue_mode_str.$lower_pval_ext;
		push @out_file_types,$pvalue_mode_str.$upper_pval_ext;
		push @out_file_types,".bgr".$site_stat_ext;
		push @out_file_types,".fgr".$site_stat_ext;
	}
}

sub check_outfiles{
	my $fn_prefix=shift;
	my @fn_suffixes=@_;
	foreach my $suff(@fn_suffixes){
		my $fn=$fn_prefix.$suff;
		return 0 unless (-e $fn);
	}
	return 1;
}

my @in_dirs;
open INFILE, "<$prj_list_fn" or die "\nUnable to open output file: $prj_list_fn!";
while(<INFILE>){
	chomp;
	s/^\s+//;
	s/\s+$//;
	$_.="/" unless /\/$/;
	if(/\S+/){
		push @in_dirs,$_;
	}
}
close INFILE;
my $ndir=@in_dirs;

sub read_site_subset{
	my ($sites_fn)=@_;
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
	return %site_subset;
}

my %bgr_sites;
my %fgr_sites;
%bgr_sites=read_site_subset($bgr_filter_fn) if(defined $bgr_filter_fn);
%fgr_sites=read_site_subset($fgr_filter_fn) if(defined $fgr_filter_fn);

#start script
my @in_project_names;
my @in_epistat_ext;
my $header_str;
my @site_pairs;
{
	my %h_site_pairs;
	my %bgr_site_counts;
	my %fgr_site_counts;

	for(my $I=0;$I<$ndir;$I++){
		my $in_dir=$in_dirs[$I];
		opendir(DIR, $in_dir) or die "can't opendir $in_dir: $!";
		my $file;
		my ($basename,$dir,$ext);
		my ($prj_name,$stat_ext);
		while (defined($file = readdir(DIR))) {
			($basename,$dir,$ext) = fileparse($file,'\..*$');
			if($ext eq $pairs_ext){
				$prj_name=$in_dir.$basename;
			}elsif($ext eq $epistat_prm_ext){
				$stat_ext=get_epistat_ext($in_dir.$basename.$ext);
			}
			last if defined($prj_name)&&($stat_ext);
		}
		closedir(DIR);
		die "\nError: Unable to find the *$pairs_ext and *$epistat_prm_ext files in the folder $in_dir!" unless defined($prj_name)&&($stat_ext);
		push @in_project_names,$prj_name;
		push @in_epistat_ext,$stat_ext;
		if(defined $prj_name){
			open INFILE, "<$in_dir$basename$pairs_ext" or die "\nUnable to open input file: $in_dir$basename$ext!";
			my $i=0;
			while(<INFILE>){
				chomp;
				s/^\s+//;
				s/\s+$//;
				if(/^\d+\t\d+/){
					my @lines=split '\t';
					my $t=1;
					if(defined $bgr_filter_fn){
						$t=0 unless defined $bgr_sites{$lines[0]};
					}
					if(defined $fgr_filter_fn){
						$t=0 unless defined $fgr_sites{$lines[1]};
					}
					if($t){
						$bgr_site_counts{$lines[0]}=[(0) x $ndir] unless defined $bgr_site_counts{$lines[0]};
						$fgr_site_counts{$lines[1]}=[(0) x $ndir] unless defined $fgr_site_counts{$lines[1]};
						$bgr_site_counts{$lines[0]}->[$I]=$lines[2] unless $bgr_site_counts{$lines[0]}->[$I]; 
						$fgr_site_counts{$lines[1]}->[$I]=$lines[3] unless $fgr_site_counts{$lines[1]}->[$I];
						$h_site_pairs{$lines[0].",".$lines[1]}=[($I,$i,$_)] if $t;
					}
					$i++;
				}elsif(/\t/){
					$header_str=$_ unless defined $header_str;
				}
			}
			close INFILE;
		}else{
			die "\nUnable to find a $pairs_ext file in the folder $in_dir! The file is required to initialize the project!";
		}
	}
	my @bgr_sites=sort {$a<=>$b} keys %bgr_site_counts;
	my @fgr_sites=sort {$a<=>$b} keys %fgr_site_counts;
	
	foreach my $bgs (@bgr_sites){
		my $nbg=0;
		for(my $I=0;$I<$ndir;$I++){
			$nbg++ if $bgr_site_counts{$bgs}->[$I];
		}
		my $ibg;
		if($nbg==1){
			#number of mutations in the site defined unambiguously
			for(my $I=0;$I<$ndir;$I++){
				if($bgr_site_counts{$bgs}->[$I]){
					$ibg=$I;
					last;
				}
			}
		}
		foreach my $fgs (@fgr_sites){
			next if $f_intragene&&($bgs==$fgs);
			push @site_pairs,[];
			if(defined $h_site_pairs{$bgs.",".$fgs}){
				@{$site_pairs[-1]}=@{$h_site_pairs{$bgs.",".$fgs}};
			}else{
				#The pair is absent in all datasets
				my $nfg=0;
				for(my $I=0;$I<$ndir;$I++){
					$nfg++ if $fgr_site_counts{$fgs}->[$I];
				}
				my $ifg;
				if($nfg==1){
					#number of mutations in the site defined unambiguously
					for(my $I=0;$I<$ndir;$I++){
						if($fgr_site_counts{$fgs}->[$I]){
							$ifg=$I;
							last;
						}
					}
				}
				my $str="$bgs\t$fgs\t";
				if(defined $ibg){
					$str.=$bgr_site_counts{$bgs}->[$ibg];
				}else{
					$str.="0";
				}
				$str.="\t";
				if(defined $ifg){
					$str.=$fgr_site_counts{$fgs}->[$ifg];
				}else{
					$str.="0";
				}
				$site_pairs[-1]->[2]=$str;
			}
		}
	}
}
my $nlines=@site_pairs;
#print site pairs
my ($basename,$dir,$ext) = fileparse($prj_list_fn,'\..*?$');
my $out_fn=$out_dir.$basename;
open OPF, ">$out_fn$pairs_ext" or die "\nUnable to open output file $out_fn$pairs_ext!";
print OPF "$header_str";
my %file_line2idx;
for(my $i=0;$i<$nlines;$i++){
	my $ra=$site_pairs[$i];
	print OPF "\n".$ra->[-1];
	$ra->[-1]=undef;
	if(defined($ra->[0])&&defined($ra->[1])){
		my $str=$ra->[0].",".$ra->[1];
		$file_line2idx{$str}=$i;
	}
}
close OPF;

sub get_lines{
	my ($infile,$data_idx,$rh_line2idx,$ra_lines)=@_;
	open INPF, "<$infile" or die "\nUnable to open input file: $infile!";
	my $i=0;
	while(<INPF>){
		s/^\s+//;
		s/\s+$//;
		if(/\S+/){
			my $idx=$rh_line2idx->{"$data_idx,$i"};
			$ra_lines->[$idx]=$_;
			$i++;
		}
	}
	close INPF;
}

sub get_epistat_lines{
	my ($infile,$data_idx,$rh_line2idx,$rh_lines,$rs_header)=@_;
	open INPF, "<$infile" or die "\nUnable to open input file: $infile!";
	while(<INPF>){
		s/^\s+//;
		s/\s+$//;
		if(/^(\d+)\t/){
			my $i=$1-1;
			my $idx=$rh_line2idx->{"$data_idx,$i"};
			$rh_lines->{$idx}=$';
		}elsif(/\S+/){
			${$rs_header}=$_ if(scalar(split "\t")>1);
		}
	}
	close INPF;
}

{
	#union of epistat files
	my %lines;
	my $header;
	for(my $I=0;$I<$ndir;$I++){
		my $prj=$in_project_names[$I];
		my $stats_ext=$in_epistat_ext[$I];
		die "\nIncomplete project: $prj! Some of required files are absent!" unless check_outfiles($prj,@out_file_types);
		get_epistat_lines($prj.$stats_ext,$I,\%file_line2idx,\%lines,\$header);
	}
	my $stats_ext=".stat.exp";
	open OPF, ">$out_fn$stats_ext" or die "\nUnable to open output file $out_fn$stats_ext!";
	print OPF "NMaxDataLines=\"$nlines\"\n$header";
	foreach my $i(sort {$a<=>$b} keys %lines){
		print OPF "\n".($i+1)."\t".$lines{$i};
	}
	close OPF;
}

foreach my $ext(@out_file_types){
	my @lines;
	for(my $I=0;$I<$ndir;$I++){
		my $prj=$in_project_names[$I];
		die "\nIncomplete project: $prj! Some of required files are absent!" unless check_outfiles($prj,@out_file_types);
		get_lines($prj.$ext,$I,\%file_line2idx,\@lines);
	}
	die "\Wrong number of total lines in *$ext files: ".scalar(@lines)." required $nlines!" unless $nlines==@lines;
	open OPF, ">$out_fn$ext" or die "\nUnable to open output file $out_fn$ext!";
	for(my $i=0;$i<$nlines;$i++){
		if(defined $lines[$i]){
			print OPF "$lines[$i]\n";
		}else{
			if(($ext eq $lower_pval_ext)||($ext eq $upper_pval_ext)){
				print OPF "$pval_default_string\n";
			}else{
				print OPF "0\n";
			}
		}
	}
	close OPF;
}
