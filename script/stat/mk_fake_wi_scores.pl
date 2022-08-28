#!/usr/bin/env perl
#This script runs mk_coevolmtx.pl for each of fake samples
#Usage: options <fdr_folder>
#options:
#	-m <method> - Method used to calculate matrix values
#		<method>=Z-score|NZ-score|covariance-like|correlation-like
#	-p <parameters> - file with parameters for the mk_coevolmtx.pl
#	-x <extension of fake epistat files> e.g. ".stat.exp.22"
#	[-P or  -C] - Call matrix inversion script
#		if both options omitted then result files contain only pseudocorrelations 
#		-P and -C - Prohibited!
#		[-P] <psicov_opt_fn> - Runs 'minvert/psicov.R' with options from the file
#		[-C] <cor2pcor_opt_fn> - Runs 'minvert/cor2pcor.R' with options from the file
#	[-o] <obs_wi_scores> - Output wi_scores for observations in the one column format. Expects a file with output of 'psicov.R' or 'cor2pcor.R' scripts
#	[-n] - No negative values. Set negative elements' values to zero
#	-e <extension for output files> 

use strict;
use Getopt::Std;
use File::Basename;

die "Set up the system variable \$EPISTAT_HOME" unless defined $ENV{EPISTAT_HOME};
my $mk_coevolmtx_cmd="$ENV{EPISTAT_HOME}/mk_coevolmtx.pl";
my $psicov_cmd="$ENV{EPISTAT_HOME}/minvert/psicov.R";
my $cor2pcor_cmd="$ENV{EPISTAT_HOME}/minvert/cor2pcor.R";

my %args;
if(!getopts('m:p:x:o:P:C:ne:',\%args)){
	die "\nError in option string!";
}
my $mk_coevolmtx_prm_ext=".mk_coevolmtx.prm";
my $psicov_out_ext=".psicov.R.out";
my $cor2pcor_out_ext=".cor2pcor.R.out";
my $fake_out_ext;
my $obs_out_ext;

my $fdrdir=$ARGV[0];
warn "\nWarning: A folder with fake samples is not accounted!" unless defined $fdrdir;
if(defined $fdrdir){
	$fdrdir=~s/\/\s*//;
	$fdrdir.="/";
}
if(defined $args{e}){
	$obs_out_ext=$args{e};
	$fake_out_ext=".fake".$obs_out_ext;
}else{
	die "\nThe extension for output files wasn't specified! Use -e option.";
}
my ($method,$prm_file,$fake_ext);
my $revmtx_script;
my $revmtx_optstr;
my $revmtx_out_ext;
my $f_nonegative=0;
die "\nNo default value for the option '-m' (method)!" unless defined $args{m};
#die "\nA method to revert of covariance matrices is undefined!" unless defined($args{P})||defined($args{C});
die "\nOptions '-C' and '-P' are mutual exclusive!" if defined($args{P})&&defined($args{C});
$f_nonegative=1 if defined $args{n};
if(defined $args{p}){
	$prm_file=$args{p};
}else{
	die "\nNo default value for the option '-p' (parameters)!";
}
if(defined $args{x}){
	$fake_ext=$args{x};
}else{
	die "\nNo default value for the option '-x' (fake_file_ext)!" unless defined $args{x};
}
if($args{m} eq "Z-score"){
	$method=$args{m};
}else{
	die "\nError: Unknown or unsupported method: $args{m}!";
}
if(defined $args{P}){
	open INPF, $args{P} or die "\nUnable to open input file: $args{P}!";
	$revmtx_script=$psicov_cmd;
	$revmtx_out_ext=$psicov_out_ext;
	while(<INPF>){
		if(/-\w/){
			$revmtx_optstr=$_;
			$revmtx_optstr=~s/^\s+//;
			$revmtx_optstr=~s/\s+$//;
			last;
		}
	}
	close INPF;
}
if(defined $args{C}){
	open INPF, $args{C} or die "\nUnable to open input file: $args{C}!";
	$revmtx_script=$cor2pcor_cmd;
	$revmtx_out_ext=$cor2pcor_out_ext;
	while(<INPF>){
		if(/-\w/){
			$revmtx_optstr=$_;
			$revmtx_optstr=~s/^\s+//;
			$revmtx_optstr=~s/\s+$//;
			last;
		}
	}
	close INPF;
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
		};
	}
}

sub conv_scores2column{
	my($fn,$rh_site_pair_idxs,$site_pairs_num,$ra_out_lines)=@_;
	open INFILE, "<$fn" or die "\nUnable to open output file: $fn!";
	@{$ra_out_lines}=();
	while(<INFILE>){
		if(/^(\d+\t\d+)/){
			chomp;
			my $sp_idx=$rh_site_pair_idxs->{$1};
			next unless defined $sp_idx;
			my @splitter=split '\t';
			my $str=$splitter[-1];
			if($str<0&&$f_nonegative){
				$str=0;
			}else{
				$str=sprintf("%.4f",$str);
			}
			$ra_out_lines->[$sp_idx]=$str;
		}
	}
	close INFILE;
	for(my $j=0;$j<$site_pairs_num;$j++){
		if(!defined $ra_out_lines->[$j]){
			$ra_out_lines->[$j]="na";
		}
	}
}

my @samples;
if(defined $fdrdir){
	opendir(DIR, $fdrdir) or die "can't opendir $fdrdir: $!";
	my $file;
	while (defined($file = readdir(DIR))) {
		my ($basename,$dir,$ext) = fileparse($file,'\..*$');
		if($ext eq $fake_ext){
			push @samples,$basename;
		}
	}
	closedir(DIR);
	die "\nNo fake samples *$fake_ext found in the $fdrdir" unless @samples>0;
}
#make parameter files for fake samples
my @lines;
my $block_mtx_idx;
my $block_mtx_ext;
my $sub_mtx_ext;
my $datasets_num=0;
my @block_name_idxs;
my @block_names;
my @site_pairs_fn;
open INFILE, "$prm_file" or die "\nUnable to open input file: $prm_file!";
my $I=0;
while(<INFILE>){
	$_=$` if(/#/);
	chomp;
	my $str=$_;
	if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
		my ($key,$value)=($1,$2);
		$datasets_num=$value if($key eq "NumberOfDatasets");
		if($key eq "BlockMtxName"){
			$str="BlockMtxName=";
			$block_mtx_idx=$I;
		}elsif($key eq "BlockMtxExt"){
			$block_mtx_ext=$value;
		}elsif($key eq "SubMtxExt"){
			$sub_mtx_ext=$value;
		}elsif($key eq "Name"){
			$str="\tName=";
			push @block_name_idxs,$I;
			push @block_names,$value;
		}elsif($key eq "Pairs"){
			push @site_pairs_fn,$value;
		}
	}
	push @lines,$str;
	$I++;
}
close INFILE;
die "\nUndefined number of datasets in the file: $prm_file!" unless $datasets_num;
die "\nMultiple datasets are not supported yet!" unless $datasets_num==1;
die "\nWrong number of datasets in the file: $prm_file!" unless $datasets_num==@block_names;
die "\nNumber of files declaring site pairs is not equal to dataset number: $prm_file!" unless $datasets_num==@site_pairs_fn;
my %change_params;
$change_params{EpiStat}=1;
$change_params{LowerPvalue}=1;
$change_params{UpperPvalue}=1;
if(@samples){
	my $tmp_fname=gen_tempname(10);
	open OPFILE, ">$tmp_fname" or die "\nUnable to open output file: $tmp_fname!";
	for($I=0;$I<@samples;$I++){
		my $sid=$samples[$I];
		my $fn=$sid;
		my $dataset_name=$site_pairs_fn[0];
		$dataset_name=~s/\.[^.]+$//;
		print OPFILE "$fn\n";
		open OPF, ">$fn".$mk_coevolmtx_prm_ext or die "\nUnable to open output file: $fn$mk_coevolmtx_prm_ext!";
		my $n=0;
		for(my $i=0;$i<@lines;$i++){
			$_=$lines[$i];
			if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
				my ($key,$value)=($1,$2);
				if(defined $change_params{$key}){
					if(index($value, $dataset_name) != -1) {
						my $ext=substr $value,length($dataset_name);
						$_=~s/\".+?$/\"/;
						if($key eq "EpiStat"){
							$_.=$fdrdir.$sid."$ext\"";
						}else{
							$_.=$fdrdir.$sid."_fake$ext\"";
						}
					}else{
						die "\nThe dataset name \"$dataset_name\" doesn't match to:\n\t$_";
					}
				}
			}
			print OPF;
			if($i==$block_mtx_idx){
				print OPF "\"$sid\"";
			}elsif($i==$block_name_idxs[$n]){
				#print OPF "\"$sid.$block_names[$n]\"";
				print OPF "\"$sid\"";
				$n++;
			}
			print OPF "\n";
		}
		close OPF;
	}
	close OPFILE;
	my $str="cat $tmp_fname|parallel ".$mk_coevolmtx_cmd." -m Z-score {}".$mk_coevolmtx_prm_ext;
	print STDERR "\n\nCalculating correlation matrices!\n$str";
	system($str);
	print_child_termination_status();
	if(defined $revmtx_script){
		my $str="cat $tmp_fname|parallel ".$revmtx_script." ".$revmtx_optstr." -f {}".$block_mtx_ext;
		print STDERR "\n\nInverting correlation matrices!\n$str";
		system($str);
		print_child_termination_status();
	}
	unlink $tmp_fname;
}
my %site_pair_idxs;
my $n=0;

for($I=0;$I<$datasets_num;$I++){
	my $block_name=$block_names[$I];
	open INFILE, "<$site_pairs_fn[$I]" or die "\nUnable to open input file: $site_pairs_fn[$I]";
	while(<INFILE>){
		$site_pair_idxs{$1}=$n++ if(/^(\d+\t\d+)/);
	}
	close INFILE;
	if(@samples){
		my @out_lines;
		print "\nJoin samples:";
		for(my $i=0;$i<@samples;$i++){
			print "\n\tsample: $i";
			my $fn=$samples[$i];
			my @tmp_lines;
			if(defined $revmtx_script){
				conv_scores2column($fn.$revmtx_out_ext,\%site_pair_idxs,$n,\@tmp_lines);
				unlink $fn.$revmtx_out_ext;
			}else{
				conv_scores2column($fn.$sub_mtx_ext,\%site_pair_idxs,$n,\@tmp_lines);
			}
			for(my $j=0;$j<$n;$j++){
				$out_lines[$j].="\t" if $i;
				$out_lines[$j].=$tmp_lines[$j];
			}
			#deleting temp files
			unlink $fn.$mk_coevolmtx_prm_ext;
			unlink $fn.$sub_mtx_ext;
			unlink $fn.$block_mtx_ext;
		}
		open OPF, ">$block_name".$fake_out_ext or die "\nUnable to open output file: $block_name$fake_out_ext!";
		print OPF "#SAMPLE_ID:";
		for(my $i=0;$i<@samples;$i++){
			print OPF $samples[$i];
			if($i<@samples-1){
				print OPF "\t";
			}else{
				print OPF "\n";
			}
		}
		for(my $j=0;$j<$n;$j++){
			print OPF "$out_lines[$j]\n";
		}
		close OPF;
	}
}

if($args{o}){
	my $wi_score_fn=$args{o};
	my @tmp_lines;
	conv_scores2column($wi_score_fn,\%site_pair_idxs,$n,\@tmp_lines);
	my ($basename,$dir,$ext) = fileparse($wi_score_fn,'\.[^\.\s]*$');
	my $fn=$dir.$basename.$obs_out_ext;
	open OPF, ">$fn" or die "\nUnable to open output file: $fn!";
	for(my $j=0;$j<$n;$j++){
		print OPF "$tmp_lines[$j]\n";
	}
	close OPF;
}
