#!/usr/bin/env perl
#This script starts processing of *.xparr batch by run_epistat.pl script
#	-i <DIR> - a folder containing input XPAR files
#   [-o] <DIR> - an output folder
#		Default="./"
#	[-e] <STR> - an extension string for XPAR files
#		Default=".xparr"
#	-E <run_epistat_opt_fn> - Runs 'run_epistat.pl' with options from the file
#	[-M] <mk_coevolmtx_opt_fn> - Runs 'mk_coevolmtx.pl' with options from the file
#	[-P] <psicov_opt_fn> - Runs 'minvert/psicov.R' with options from the file
#	[-C] <cor2pcor_opt_fn> - Runs 'minvert/cor2pcor.R' with options from the file
#	[-r] Rerun. Proceed previous unsuccessfully finished job
#Params:
#<parameters> - epistat.pl configuration file
#<parameters> - run_epistat.pl configuration file
#[<mk_coevolmtx_prm>] - The template of mk_coevolmtx.pl configuration file, -M option required! 
use strict;
use Getopt::Std;
use File::Basename;
use File::Path qw(make_path rmtree);

my $run_epistat_cmd="$ENV{EPISTAT_HOME}/run_epistat.pl";
my $mk_coevolmtx_cmd="$ENV{EPISTAT_HOME}/mk_coevolmtx.pl";
my $psicov_cmd="$ENV{EPISTAT_HOME}/minvert/psicov.R";
my $cor2pcor_cmd="$ENV{EPISTAT_HOME}/minvert/cor2pcor.R";

my %args;
if(!getopts('i:E:M:P:C:e:o:r',\%args)){
	die "\nError in option string!";
}
my $f_rerun=0;
$f_rerun=1 if defined $args{r};
my $xpar_dir=$args{i};
die "\nThe input folder with XPAR files is not specified! Use -i option!" unless defined $xpar_dir;
$xpar_dir=~s/\s+$//;
$xpar_dir=~s/\/$//;
$xpar_dir.="/";
my $xpar_ext=".xparr";
$xpar_ext=$args{e} if defined $args{e};
my $out_dir="./";
if(defined $args{o}){
	$out_dir=$args{o};
	$out_dir=~s/\s+$//;
	$out_dir=~s/\/$//;
	$out_dir.="/";
} 

my $run_epistat_opts;
if(defined $args{E}){
	$run_epistat_opts="";
	open INPF, $args{E} or die "\nUnable to open input file: $args{E}!";
	while(<INPF>){
		if(/-\w/){
			$run_epistat_opts=$_;
			$run_epistat_opts=~s/^\s+//;
			$run_epistat_opts=~s/\s+$//;
			$run_epistat_opts=~s/-r// if $f_rerun;
			last;
		}
	}
	close INPF;
}

my $mk_coevolmtx_opts;
my $mk_coevolmtx_prm;
if(defined $args{M}){
	open INPF, $args{M} or die "\nUnable to open input file: $args{M}!";
	while(<INPF>){
		if(/-\w/){
			$mk_coevolmtx_opts=$_;
			$mk_coevolmtx_opts=~s/^\s+//;
			$mk_coevolmtx_opts=~s/\s+$//;
			last;
		}
	}
	close INPF;
	$mk_coevolmtx_prm=$ARGV[2];
	die "The template configuration file for mk_coevolmtx.pl script is not accounted!" unless defined $mk_coevolmtx_prm;
}
my $psicov_opts;
if(defined $args{P}){
	open INPF, $args{P} or die "\nUnable to open input file: $args{P}!";
	$psicov_opts="";
	while(<INPF>){
		if(/-\w/){
			$psicov_opts=$_;
			$psicov_opts=~s/^\s+//;
			$psicov_opts=~s/\s+$//;
			last;
		}
	}
	close INPF;
}
my $cor2pcor_opts;
if(defined $args{C}){
	open INPF, $args{C} or die "\nUnable to open input file: $args{C}!";
	$cor2pcor_opts="";
	while(<INPF>){
		if(/-\w/){
			$cor2pcor_opts=$_;
			$cor2pcor_opts=~s/^\s+//;
			$cor2pcor_opts=~s/\s+$//;
			last;
		}
	}
	close INPF;
}

my @xpar_fn;
opendir(DIR, $xpar_dir) or die "can't opendir $xpar_dir: $!";
my $file;
while (defined($file = readdir(DIR))) {
	my ($basename,$dir,$ext) = fileparse($file,'\.[^\.]*$');
	if($ext eq $xpar_ext){
		push @xpar_fn,$basename;
	}
}
closedir(DIR);

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

sub gen_tempname{
	my $nchar=shift;
	my @chars = ( "A" .. "Z", "a" .. "z", 0 .. 9 );
	return join("", @chars[ map { rand @chars } ( 1 .. $nchar ) ]);
}

foreach my $fn(@xpar_fn){
	my $str=$out_dir.$fn."_dat/";
	rmtree([ "$str" ]) if((-d $str)&&(!$f_rerun));
	make_path($str) unless -d $str;
	unless(-e $str.$fn.$xpar_ext){
		my $str="cp $xpar_dir$fn$xpar_ext ".$str;
		system($str);
		print_child_termination_status();
	}
}
#start run_epistat.pl script
print STDERR "Starting run_epistat.pl:";
foreach my $fn(@xpar_fn){
	my $dat_dir=$out_dir.$fn."_dat/";
	my $str=$run_epistat_cmd." -x $dat_dir$fn$xpar_ext $run_epistat_opts";
	$str.=" -r" if $f_rerun;
	$str.=" $ARGV[0] $ARGV[1]";
	print STDERR "\n$str";
	system($str);
	print_child_termination_status();
}
#start mk_coevolmtx.pl script
my $block_mtx_ext;
my $tmp_fname=gen_tempname(10);
my $n=0;
if(defined $mk_coevolmtx_opts){
	my ($basename,$dir,$mk_coevolmtx_prm_ext) = fileparse($mk_coevolmtx_prm,'\..+?\s*$');
	my @lines;
	my $block_mtx_idx;
	my $sub_mtx_ext;
	my $datasets_num=0;
	my @block_name_idxs;
	my @block_names;
	open INFILE, "$mk_coevolmtx_prm" or die "\nUnable to open input file: $mk_coevolmtx_prm!";
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
			}
		}
		push @lines,$str;
		$I++;
	}
	close INFILE;
	die "\nUndefined number of datasets in the file: $mk_coevolmtx_prm!" unless $datasets_num;
	die "\nMultiple datasets are not supported yet!" unless $datasets_num==1;
	die "\nWrong number of datasets in the file: $mk_coevolmtx_prm!" unless $datasets_num==@block_names;
	open OPFILE, ">$tmp_fname" or die "\nUnable to open output file: $tmp_fname!";
	foreach my $sid(@xpar_fn){
		my $fn=$out_dir.$sid."_dat/".$sid;
		unless($f_rerun&&(-e $fn.$block_mtx_ext)){
			print OPFILE "$sid\n";
			open OPF, ">$sid".$mk_coevolmtx_prm_ext or die "\nUnable to open output file: $sid$mk_coevolmtx_prm_ext!";
			$I=0;
			for(my $i=0;$i<@lines;$i++){
				$_=$lines[$i];
				if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
					my ($key,$value)=($1,$2);
					if($key eq "EpiStat" ||
						$key eq "BgrSiteMarginStat"||
						$key eq "FgrSiteMarginStat"||
						$key eq "Mean"||
						$key eq "Pairs"||
						$key eq "Variance"||
						$key eq "LowerPvalue"||
						$key eq "UpperPvalue"||
						$key eq "OrdPairsCovariance"
						){
						$_=~s/\"[^\s.]*\./\"$fn\./;
					}
				}
				print OPF;
				if($i==$block_mtx_idx){
					print OPF "\"$fn\"";
				}elsif($i==$block_name_idxs[$I]){
					#print OPF "\"$sid.$block_names[$I]\"";
					print OPF "\"$fn\"";
					$I++;
				}
				print OPF "\n";
			}
			close OPF;
			$n++;
		}
	}
	close OPFILE;
	if($n){
		my $str="cat $tmp_fname|parallel $mk_coevolmtx_cmd $mk_coevolmtx_opts {}".$mk_coevolmtx_prm_ext;
		print STDERR "\n\nCalculating correlation matrices!\n$str";
		system($str);
		print_child_termination_status();
		foreach my $fn(@xpar_fn){
			unlink $fn.$mk_coevolmtx_prm_ext;
		}
	}
	unlink $tmp_fname;
}
#start minvert/psicov.R
if(defined($psicov_opts)){
	open OPFILE, ">$tmp_fname" or die "\nUnable to open output file: $tmp_fname!";
	foreach my $sid(@xpar_fn){
		print OPFILE "$sid\n";
	}
	close OPFILE;
	my $str="cat $tmp_fname|parallel $psicov_cmd -f $out_dir"."{}"."_dat/{}"."$block_mtx_ext $psicov_opts";
	print STDERR "\n\nCalculating inverse matrices using psicov.R!\n$str";
	system($str);
	print_child_termination_status();
	unlink $tmp_fname;
}
#start minvert/cor2pcor.R
if(defined($cor2pcor_opts)){
	open OPFILE, ">$tmp_fname" or die "\nUnable to open output file: $tmp_fname!";
	foreach my $sid(@xpar_fn){
		print OPFILE "$sid\n";
	}
	close OPFILE;
	my $str="cat $tmp_fname|parallel $cor2pcor_cmd -f $out_dir"."{}"."_dat/{}"."$block_mtx_ext $cor2pcor_opts";
	print STDERR "\n\nCalculating inverse matrices using cor2pcor.R!\n$str";
	system($str);
	print_child_termination_status();
	unlink $tmp_fname;
}

