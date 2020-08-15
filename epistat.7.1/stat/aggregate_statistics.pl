#!/usr/bin/env perl
#This script aggregates data from raw files into a table
#Usage: [options] <input_folder>
#<input_folder> - A folder with input files
#[options]
#	-p <fn_prefix_list> - A file name which contains a prefixes of files with data. Each prefix corresponds to a row in the final table
#	-o <file_name> - Output file name
#	-n <file_name_parsing_rules> - A file with input a file's name parsing rules
#	-c <content_parsing_rules> - A file with input a file's content parsing rules
#<file_name_parsing_rules> Grammatics:
#	<SeparatorDecl>
#	<AssignmentTable>
#	<SeparatorDecl> <- Separator=\"STRING\" - Splits a file name suffix by <STR>
#	<AssignmentTable> <- <AssignmentRow>+
#	<AssignmentRow> <- <ElementIdx>\t<TableColumn>
#	<ElementIdx> <- UINT - An index of an element in array containing of the splitted file name suffix
#	<TableColumn> <- UINT - A column number in the table where element should be placed
#<content_parsing_rules> Grammatics:
#	<AssignmentTable>
#	<AssignmentTable> <- <AssignmentRow>+
#	<AssignmentRow> <- <FileSearchString>\t<ElementSearchString>\t<ElementType>\t<TableColumn>
#	<FileSearchString> <- \"STRING\" - String used for identification of input file name to apply a rule. Empty string applies a rule to any file.
#	<ElementSearchString> <- \"STRING\" - String used for element searching. Element starts just after the search string
#	<ElementType> <- (NUMBER|LOGIC|WORD|STRING)

use strict;
use Getopt::Std;
use File::Basename;
use Class::Struct;

my %args;
if(!getopts('o:n:c:p:',\%args)){
	die "\nError in option string!";
}

my $indir=$ARGV[0];
die "\nThe folder with input data is required!\n\tUse '.' to specify the current folder!" unless defined $indir;
$indir.="/" unless $indir=~/\/$/;
my $fnprefix_fn=$args{p};
my $out_fn=$args{o};
my $name_parsing_fn=$args{n};
my $content_parsing_fn=$args{c};
my @infiles;
my %table;

struct TableElement=>{
	column_idx => '$',
	value => '$'
};

#begin script
if(defined $fnprefix_fn){
	open INPF, "<$fnprefix_fn" or die "\nUnable to open input file $fnprefix_fn!";
	while(<INPF>){
		chomp;
		if(/\S/){
			$_=~s/^\s+//;
			$_=~s/\s+$//;
			push @infiles, [($_,1)];
		}
	}
	close INPF;
}
{
	my $file;
	opendir(DIR, $indir) or die "can't opendir $indir: $!";
	while (defined($file = readdir(DIR))) {
		my ($basename,$dir,$ext) = fileparse($file,'\.[^.]*$');
		my $t=0;
		push @infiles, [($basename.$ext,0,[])] unless ($basename eq "")||($basename=~/^\.+$/)||!defined($basename);
	}
	closedir(DIR);
}
@infiles=sort {$a->[0] cmp $b->[0]} @infiles;
my $key;
my $n;
for(my $i=0;$i<@infiles;$i++){
	my $rfn=$infiles[$i];
	if(defined $fnprefix_fn){
		if($rfn->[1]){
			$key=$rfn->[0];
			$n=length $key;
			$table{$key}=[] unless defined $table{$key};
		}else{
			my $str=substr $rfn->[0],0,$n;
			push @{$table{$key}},$i if($str eq $key);
		}
	}else{
		$key=$rfn->[0];
		$table{$key}=[($i)] unless defined $table{$key};
	}	
}

my $ncol=0;

if(defined $name_parsing_fn){
	my @fn_parse;
	my $fn_spliter;
	open INPF, "<$name_parsing_fn" or die "\nUnable to open input file $name_parsing_fn!";
	while(<INPF>){
		$_=$` if(/#/);
		chomp;
		s/\s+$//;
		if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
			my ($key,$value)=($1,$2);
			if($key eq "Separator"){
				$fn_spliter=$value;
				$fn_spliter=~s/\./\\\./;
			}
		}elsif(/^(\d+)\t(\d+)/){
			push @fn_parse,[($1,$2)];
			$ncol=$2 if $ncol<$2;
		}
	}
	close INPF;
	my @fn_parse=sort {$a->[0] <=> $b->[0]} @fn_parse;
	my %tmp_table;
	foreach my $key(keys %table){
		foreach my $i(@{$table{$key}}){
			my $fn=$infiles[$i]->[0];
			my $str=$key;
			my @line=split /$fn_spliter/,$fn;
			foreach my $r(@fn_parse){
				my $j=$r->[0];
				if($j<=@line){
					$str.="\\".$line[$j-1];
					my $ne=TableElement->new();
					$ne->column_idx($r->[1]-1);
					$ne->value($line[$j-1]);
					push @{$infiles[$i]->[2]},$ne;
				}else{
					die "\nError index $j of a file name element in $name_parsing_fn: out of range in the $fn input file!"
				}
			}
			$tmp_table{$str}=[] unless(defined $tmp_table{$str});
			push @{$tmp_table{$str}},$i;
		}
	}
	%table=%tmp_table;
}

if(defined $content_parsing_fn){
	my %fn_parse;
	open INPF, "<$content_parsing_fn" or die "\nUnable to open input file $content_parsing_fn!";
	while(<INPF>){
		$_=$` if(/#/);
		chomp;
		s/\s+$//;
		if(/^\"(.*?)\"\t\"(.+?)\"\t(NUMBER|LOGIC|WORD|STRING)\t(\d+)/){
			$fn_parse{$1}=[] unless defined $fn_parse{$1};
			push @{$fn_parse{$1}},[($2,$3,$4)];
			$ncol=$4 if $ncol<$4;
		}
	}
	close INPF;
	foreach my $key(keys %table){
		foreach my $i(@{$table{$key}}){
			my $fn=$infiles[$i]->[0];
			foreach my $ss(keys %fn_parse){
				if(($ss eq "")||(index($fn, $ss) != -1)){
					open INPF, "<$indir".$fn or die "\nUnable to open input file $fn!";
					while(<INPF>){
						chomp;
						if(/\S/){
							foreach my $r(@{$fn_parse{$ss}}){
								my $idx=index($_,$r->[0]);
								if($idx!=-1){
									my $val;
									my $str=substr $_,$idx+length($r->[0]);
									$str=~s/^\s+//;
									if($r->[1] eq "STRING"){
										$val=$1 if $str=~/\"(.*?)\"/;
									}else{
										$str=~s/\s+$//;
										$val=$1 if $str=~/(\S+)/;
									}
									if(defined $val){
										my $ne=TableElement->new();
										$ne->column_idx($r->[2]-1);
										$ne->value($val);
										push @{$infiles[$i]->[2]},$ne;
									}
								}
							}
						}
					}
					close INPF;
				}
			}
		}
	}
}

if(defined $out_fn){
	open OUTF, ">>$out_fn" or die "\nUnable to open output file $out_fn!";
}else{
	open(OUTF, ">&STDOUT") or die "Couldn't dup STDOUT: $!";
}

my @keys_srt=sort {$a cmp $b} keys %table;
foreach $key(@keys_srt){
	my $str=$key;
	$str=~s/\\.+$//;
	my @tline=("") x $ncol;
	foreach my $i(@{$table{$key}}){
		foreach my $ne(@{$infiles[$i]->[2]}){
			my $j=$ne->column_idx;
			my $val=$ne->value;
			if($tline[$j] ne ""){
				$tline[$j]=$tline[$j].";".$val if $tline[$j] ne $val;
			}else{
				$tline[$j]=$val;
			}
		}
	}
	print OUTF "\n$str";
	$str=join "\t",@tline;
	print OUTF "\t$str";
}
close OUTF if(defined $out_fn);