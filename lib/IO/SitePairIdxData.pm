package IO::SitePairIdxData;
#This package is a collection of methods for reading from and writing to a file the data indexed by a number of site pairs in a '*.site_pair' file.
#The file format of in/out files:
#[<column_names>]
#<index>\t<data>
#	<index>=uint>0
#	<data> - tab delimited
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
# Symbols to autoexport (:DEFAULT tag)
@EXPORT = qw(); 
@EXPORT_OK = qw(hread_file aread_file airead_file);

sub read_data_lines_number{
	my $fh=shift;
	my ($key,$value);
	while(readline($fh)){
		if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
			($key,$value)=($1,$2);
			last if($key eq "NMaxDataLines");
		}
	}
	return $value;
}

#data extracted from the input file are stored column-wise in the array of hash references @{$ra_data}.
#the columns to extract could be specified in the '@{$ra_sel_cols}' array starting from 1 (0 is an index column).
#returns the maximal index
sub hread_file{
	my ($in_fname,$ncol,$ra_data,$ra_sel_cols)=@_;
	open INPF, "<$in_fname" or die "\nUnable to open input file:$in_fname!";
	my $n=$ncol-1;
	my %scols;
	if(defined $ra_sel_cols){
		$n=0;
		for(my $i=0;$i<@{$ra_sel_cols};$i++){ 
			my $cn=$ra_sel_cols->[$i];
			if($cn>0){
				$scols{$cn}=$i;
				$n++;
			}
		}
	}
	die "\nhread_file(): The number of hash references not matches to the number of requested columns!" unless @{$ra_data}==$n;
	for(my $i=0;$i<$n;$i++){
		%{$ra_data->[$i]}=();
	}
	my $nlines=read_data_lines_number(*INPF);
	die "\nhread_file(): Unable to identify the number of data rows in $in_fname!" unless defined $nlines;
	while(<INPF>){
		chomp;
		my @line=split "\t",$_,-1;
		$line[0]=~s/^\s+//;
		$line[0]=~s/\s+$//;
		if($line[0]=~m/^\d+/){
			my $I=$line[0]-1;
			die "\nError in hread_file(): zero is not allowed value for site pair indexing!" unless $I>=0;
			die "\nhread_file(): Unexpected number of columns in the input file $in_fname!" unless @line==$ncol;
			die "\nhread_file(): The wrong number of lines declared in $in_fname!" unless $I<$nlines;
			for(my $j=1;$j<$ncol;$j++){
				my $J=$j-1;
				if(defined $ra_sel_cols){
					if(defined $scols{$j}){
						$J=$scols{$j};
					}else{$J=undef;}
				}
				$ra_data->[$J]->{$I}=$line[$j] if defined $J;
			}
		}
	}
	close INPF;
	return $nlines;
}

#data extracted from the input file are stored column-wise in the array of array references @{$ra_data}.
#the columns to extract could be specified in the '@{$ra_sel_cols}' array starting from 1 (0 is an index column).
#for indexes which are absent in the input file the stub values are used from an array specified by @{$ra_stub_values}.
#The number of elements in the '@{$ra_stub_values}' have to be equal to the number of elements in '@{$ra_sel_cols}' if defined, or '$ncol-1' otherwise.
#returns the  the maximal index
sub aread_file{
	my ($in_fname,$ncol,$ra_data,$ra_stub_values,$ra_sel_cols)=@_;
	open INPF, "<$in_fname" or die "\nUnable to open input file:$in_fname!";
	my $n=$ncol-1;
	my %scols;
	if(defined $ra_sel_cols){
		$n=0;
		for(my $i=0;$i<@{$ra_sel_cols};$i++){ 
			my $cn=$ra_sel_cols->[$i];
			if($cn>0){
				$scols{$cn}=$i;
				$n++;
			}
		}
	}
	die "\nError in aread_file(): The number of stub values have to be equal to $n!" unless @{$ra_stub_values}==$n;
	die "\naread_file(): The number of array references ".scalar(@{$ra_data})." not matches to the number of requested columns $n!" unless @{$ra_data}==$n;
	for(my $i=0;$i<$n;$i++){
		@{$ra_data->[$i]}=();
	}
	my $I=0;
	my $nlines=read_data_lines_number(*INPF);
	while(<INPF>){
		chomp;
		my @line=split "\t",$_,-1;
		$line[0]=~s/^\s+//;
		$line[0]=~s/\s+$//;
		if($line[0]=~m/^\d+/){
			die "\naread_file(): Unexpected number of columns in the input file $in_fname!" unless @line==$ncol;
			die "\nError in aread_file(): zero is not allowed value for site pair indexing!" unless $line[0]>0;
			die "\naread_file(): Wrong number of lines declared in $in_fname!" unless $I<$nlines;
			for(my $j=1;$j<$ncol;$j++){
				my $J=$j-1;
				if(defined $ra_sel_cols){
					if(defined $scols{$j}){
						$J=$scols{$j};
					}else{$J=undef;}
				}
				if(defined $J){
					unless($line[0]-1==$I){
						for(my $i=$I;$i<$line[0]-1;$i++){
							push @{$ra_data->[$J]},$ra_stub_values->[$J];
						}
					}
					push @{$ra_data->[$J]},$line[$j];
				}
			}
			$I=$line[0];
		}
	}
	if($I<$nlines){
		for(my $j=1;$j<$ncol;$j++){
			my $J=$j-1;
			if(defined $ra_sel_cols){
				if(defined $scols{$j}){
					$J=$scols{$j};
				}else{$J=undef;}
			}
			if(defined $J){
				for(my $i=$I;$i<$nlines;$i++){
					push @{$ra_data->[$J]},$ra_stub_values->[$J];
				}
			}
		}
	}
	close INPF;
	return $nlines;
}

#data extracted from the input file are stored column-wise in the array of array references @{$ra_data}.
#the corresponding indices are stored in the array @{$ra_indices}.
#the columns to extract could be specified in the '@{$ra_sel_cols}' array starting from 1 (0 is an index column).
#returns the maximal index
sub airead_file{
	my ($in_fname,$ncol,$ra_indices,$ra_data,$ra_sel_cols)=@_;
	open INPF, "<$in_fname" or die "\nUnable to open input file:$in_fname!";
	my $n=$ncol-1;
	my %scols;
	if(defined $ra_sel_cols){
		$n=0;
		for(my $i=0;$i<@{$ra_sel_cols};$i++){ 
			my $cn=$ra_sel_cols->[$i];
			if($cn>0){
				$scols{$cn}=$i;
				$n++;
			}
		}
	}
	die "\nairead_file(): The number of array references not matches to the number of requested columns!" unless @{$ra_data}==$n;
	for(my $i=0;$i<$n;$i++){
		@{$ra_data->[$i]}=();
	}
	@{$ra_indices}=();
	my $nlines=read_data_lines_number(*INPF);
	while(<INPF>){
		chomp;
		my @line=split "\t",$_,-1;
		$line[0]=~s/^\s+//;
		$line[0]=~s/\s+$//;
		if($line[0]=~m/^\d+/){
			my $I=$line[0]-1;
			push @{$ra_indices},$I;
			die "\nairead_file(): Unexpected number of columns in the input file $in_fname!" unless @line==$ncol;
			die "\nError in airead_file(): zero is not allowed value for site pair indexing!" if $I<0;
			for(my $j=1;$j<=$ncol;$j++){
				my $J=$j-1;
				if(defined $ra_sel_cols){
					if(defined $scols{$j}){
						$J=$scols{$j};
					}else{$J=undef;}
				}
				if(defined $J){
					push @{$ra_data->[$J]},$line[$j];
				}
			}
		}
	}
	close INPF;
	return $nlines;
}

1;
