#!/usr/bin/env perl
#The test script for Epistat distribution
use strict;
use File::Basename;
use Test::More tests => 6;

my $dirname = dirname(__FILE__);
$dirname.="/data";
$ENV{EPISTAT_HOME}="blib/script";
$ENV{EPISTAT_LIB}="blib/lib";

sub print_child_termination_status{
	my $cmd=shift;
	if ($? == -1) {
		diag("$cmd\n\tfailed to execute: $!\n");
	}elsif ($? & 127) {
		my $str=sprintf "\n\tchild died with signal %d, %s coredump\n",
					($? & 127),  ($? & 128) ? 'with' : 'without';
		diag($cmd.$str);
	}else {
		if($? >> 8!=0){
			my $str=sprintf "\n\tchild exited with value %d\n", $? >> 8;
			diag($cmd.$str);
		}
	}
}
sub gen_tempname{
	my $nchar=shift;
	my @chars = ( "A" .. "Z", "a" .. "z", 0 .. 9 );
	return join("", @chars[ map { rand @chars } ( 1 .. $nchar ) ]);
}
my $tmp_sample2xparr_prm="$dirname/".gen_tempname(5).".sample2xparr.prm";
open OPF, ">$tmp_sample2xparr_prm" or die "\nUnamble to open output file $tmp_sample2xparr_prm!";
open INPF, "<$dirname/mixed.sample2xparr.1.prm" or die "\nUnable to open input file $dirname/mixed.sample2xparr.prm!";
while(<INPF>){
	if(/XPAR=\"(.+?)\"/){
		print OPF "XPAR=\"$dirname/$1\"\n";
	}elsif(/FgrIdx2Site=\"(.+?)\"/){
		print OPF "FgrIdx2Site=\"$dirname/$1\"\n";
	}elsif(/FgrIdx2Branch=\"(.+?)\"/){
		print OPF "FgrIdx2Branch=\"$dirname/$1\"\n";
	}else{
		print OPF $_;
	}
}
close OPF;
close INPF;

my $tmp_block_mtx="$dirname/".gen_tempname(5);
open OPF, ">$tmp_block_mtx.block.mtx" or die "\nUnamble to open output file $tmp_block_mtx.block.mtx!";
open INPF,"<$dirname/mixed.block.mtx" or die "\nUnable to open input file $dirname/mixed.block.mtx!";
while(<INPF>){
	if(/OutDir=/){
		print OPF "OutDir=\"$dirname/\"\n";
	}else{
		print OPF;
	}
}
close OPF;
close INPF;

#start tests
my $t_run_epistat=1;
#test epistat.pl
my $cmd_str="$ENV{EPISTAT_HOME}/epistat.pl -x $dirname/mixed.xparr -m 2 $dirname/mixed.epistat.prm";
print STDERR "\n$cmd_str";
system($cmd_str);
print_child_termination_status($cmd_str);
unless(
	(-s "$dirname/mixed.stat.exp.0.804")&&
	(-s "$dirname/mixed.site_pairs")&&
	(-s "$dirname/mixed.intra.fgr.idx2branch")&&
	(-s "$dirname/mixed.intra.fgr.idx2site")&&
	(-s "$dirname/mixed.intra.fgr.mmtx")
	){
	fail("epistat.pl");
	$t_run_epistat=0;
}else{
	pass("epistat.pl");
}
unlink "$dirname/mixed.stat.exp.0.804" if -e "$dirname/mixed.stat.exp.0.804";
unlink "$dirname/mixed.site_pairs" if -e "$dirname/mixed.site_pairs";
unlink "$dirname/mixed.intra.fgr.idx2branch" if -e "$dirname/mixed.intra.fgr.idx2branch";
unlink "$dirname/mixed.intra.fgr.idx2site" if -e "$dirname/mixed.intra.fgr.idx2site";
unlink "$dirname/mixed.intra.fgr.mmtx" if -e "$dirname/mixed.intra.fgr.mmtx";

#test shuffle_incidence_matrix.R
$cmd_str="$ENV{EPISTAT_HOME}/shuffle_incidence_matrix.R $dirname/mixed.intra.fgr.mmtx.1 $dirname";
print STDERR "\n$cmd_str";
system($cmd_str);
print_child_termination_status($cmd_str);
unless(-s "$dirname/1.mlist"){
	fail("shuffle_incidence_matrix.R");
	$t_run_epistat=0;
}else{
	pass("shuffle_incidence_matrix.R");
	#test sample2xparr.pl
	$cmd_str="$ENV{EPISTAT_HOME}/sample2xparr.pl -t $dirname/1.mlist $tmp_sample2xparr_prm>$dirname/1.xparr";
	print STDERR "\n$cmd_str";
	system($cmd_str);
	print_child_termination_status($cmd_str);
	unless(-s "$dirname/1.xparr"){
		fail("sample2xparr.pl");
	}else{
		pass("sample2xparr.pl");
	}
	unlink "$dirname/1.xparr" if -e "$dirname/1.xparr";
}
unlink "$dirname/1.mlist" if -e "$dirname/1.mlist";
unlink $tmp_sample2xparr_prm;

#test run_epistat.pl
SKIP: {
	skip "An installation error of epistat.pl or shuffle_incidence_matrix.R", 1 if $t_run_epistat==0;
	$cmd_str="$ENV{EPISTAT_HOME}/run_epistat.pl -x $dirname/mixed.xparr -m 2 -p 10 $dirname/mixed.epistat.prm $dirname/run_epistat.prm";
	print STDERR "\n$cmd_str";
	system($cmd_str);
	print_child_termination_status($cmd_str);
	my @out_file_types;
	{
		push @out_file_types,".intragene.unord_pairs.lower.pvalue";
		push @out_file_types,".intragene.unord_pairs.upper.pvalue";
		push @out_file_types,".stat.exp.0.804";
		push @out_file_types,".site_pairs";
		push @out_file_types,".mean";
		push @out_file_types,".var";
		push @out_file_types,".moments12";
		push @out_file_types,".intra.fgr.idx2branch";
		push @out_file_types,".intra.fgr.idx2site";
		push @out_file_types,".intra.fgr.mmtx";
	}
	my $t=1;
	foreach my $fn_ext(@out_file_types){
		unless(-s "$dirname/mixed$fn_ext"){
			$t=0;
			last;
		}
	}
	unless($t){
		fail("run_epistat.pl");
	}else{
		pass("run_epistat.pl");
		#system("rm -r $dirname/mixed");
	}
	foreach my $fn_ext(@out_file_types){
		unlink "$dirname/mixed$fn_ext" if -e "$dirname/mixed$fn_ext";
	}
	unlink "$dirname/mixed.sample2xparr.prm" if -e "$dirname/mixed.sample2xparr.prm";
}
#test cor2pcor.R
$cmd_str="$ENV{EPISTAT_HOME}/minvert/cor2pcor.R -f $tmp_block_mtx.block.mtx -n 0 -l 0.9 -s .n0.l90.cor2pcor.R.out";
print STDERR "\n$cmd_str";
system($cmd_str);
print_child_termination_status($cmd_str);
unless(-s "$dirname/mixed.n0.l90.cor2pcor.R.out"){
	fail("minvert/cor2pcor.R");
}else{
	pass("minvert/cor2pcor.R");
}
unlink "$dirname/mixed.n0.l90.cor2pcor.R.out" if -e "$dirname/mixed.n0.l90.cor2pcor.R.out";

#test psicov.R
$cmd_str="$ENV{EPISTAT_HOME}/minvert/psicov.R -f $tmp_block_mtx.block.mtx -n 0 -d 0.05 -s .n0.d05.psicov.R.out";
print STDERR "\n$cmd_str";
system($cmd_str);
print_child_termination_status($cmd_str);
unless(-s "$dirname/mixed.n0.d05.psicov.R.out"){
	fail("minvert/psicov.R");
}else{
	pass("psicov.R");
}
unlink "$dirname/mixed.n0.d05.psicov.R.out" if -e "$dirname/mixed.n0.d05.psicov.R.out";
unlink "$tmp_block_mtx.block.mtx";