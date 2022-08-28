#!/usr/bin/env perl
#This script prepares list of edges for a site graph
#Usage: options <Matrix>
#options:
#	[-u] - Unordered site pairs. 
#		If the option is accounted, it extracts from the matrix file only those pairs which first coordinate was less than the second 
#	[-f] <Filters> - Filters for matrix rows
#		<Filters>="<Filter>(,<Filter>)*"
#		<Filter>=<ColIdx>((<|>)=?|!=|==)(int|float)
#		<ColIdx>=uint
#	[-a] <FN> - Add columns to the matrix file. File in the argument contains tables with data to add.
#		Format:
#			(Path="ProjectPath")?
#			Pairs="FileName" - list of all pairs of markers (sites)
#			<FileName>\t<SkipNLines>\t<KeyColIdx>\t<ColumnList>\t<ColumnNameList>\t<StubValuesList>\n
#				<ColumnList>=<ColIdx><Filter>?(,<ColIdx><Filter>)*
#				<ColIdx>=uint
#				<SkipNLines>=uint|"" - number of lines to skip. Default=0
#				<KeyColIdx>=uint|"" - a column that contains indices of rows from the <Pairs> file
#				<StubValuesList>=number(,number)*||"" - list of values in columns for rows missed in the data file.
#					It is required if <KeyColIdx> is not equal to "" 
#				<Filter>= (h(ide)?)?((<|>)=?|!=|==)(int|float)
#				<ColumnNameList>=string(,string)*
#	[-s] <uint> - Skip specified number of lines from the top of the matrix file
#Matrix
#	Format:
#		(<InfoLine>\n)*
#		<SiteID>\t<SiteID>(\tColumnValue)*\n
#			<InfoLine>=string
#			<SiteId>=uint
#			ColumnValue=(int|float|string)
use strict;
use Getopt::Std;
use Class::Struct;
use FileHandle;
#constants
use constant {
	FILE_NAME =>0,
	SKIP_NLINES => 1,
	KEY_COL_IDX => 2,
	COL_LIST => 3,
	COL_NAME_LIST => 4,
	STUB_VAL_LIST => 5,
};

my %args;
if(!getopts('a:us:f:',\%args)){
	die "\nError in option string!";
}
my $f_ord_pairs=1;
$f_ord_pairs=0 if $args{u};
my $flist=$args{a};
my $skip_nlines=0;
if(defined $args{s}){
	$skip_nlines=$args{s};
	die "\nWrong argument specified for -s option: $skip_nlines!" unless $skip_nlines=~/^\d+/;
}
my @mtx_filters;
if(defined $args{f}){
	my @line=split ',',$args{f};
	foreach my $str(@line){
		if($str=~/(\d+)((<|>)=?|!=|==)(\S+)/){
			my $idx=$1-1;
			die "\nNumbering of columns starts from 1!" unless $idx>=0;
			push @mtx_filters,[($idx,$2,$4)];
		}else{
			die "\nUnable to parse selection condition: $str!";
		}
	}
}
my %pair2idx;
my %add_data;
my @matrix;
struct DataSourceInfo =>{
	name => '$',
	key => '$',
	skip_nlines => '$',
	col_indices => '@',
	col_names => '@',
	col_thresholds => '@',
	f_hide_columns => '@',
	stub_values => '@'
};

struct DataRow => {
	index => '$',
	data => '@'
};

sub extract_columns{
	my ($ra_line,$ra_indices,$ra_out)=@_;
	my $ncol=@{$ra_indices};
	for(my $i=0;$i<$ncol;$i++){
		my $idx=$ra_indices->[$i];
		return 0 unless $idx<@{$ra_line};
		push @{$ra_out},$ra_line->[$idx];
	}
	return 1;
}

sub get_line_from_sources{
	my ($ra_fh,$ra_data_sources,$rs_line_index,$ra_out)=@_;
	my $T=1;
	my $n=0;
	my $line_index=${$rs_line_index};
	for(my $i=0;$i<@{$ra_fh};$i++){
		my $fh=$ra_fh->[$i];
		my $data_source=$ra_data_sources->[$i];
		my $name=$data_source->name;
		my $pos;
		if(defined $data_source->key){
			$pos=tell $fh;
		}
		while(<$fh>){
			chomp;
			if(/\S+/){
				my @line=split '\t',$_,-1;
				my $key=$data_source->key;
				if(defined $key){
					$key=$line[$key]-1; #read the current line index from the key column if it is defined
				}else{
					$key=$line_index;
				}
				unless($key==$line_index){
					for(my $i=0;$i<@{$data_source->col_indices};$i++){
						my $ci=$data_source->col_indices($i);
						$line[$ci]=$data_source->stub_values($i);
					}
					#return back a reading position in the file buffer of $fh
					seek $fh, $pos, 0 or die "Couldn't seek to $pos: $!\n" if defined $pos;
				}
				if(@{$data_source->col_thresholds}){
					$T=$T&&check_row_conditions(\@line,$data_source->col_indices,$data_source->col_thresholds);
				}
				my @out;
				extract_columns(\@line,$data_source->col_indices,\@out)
					or die "\nIndex out of range in source $name!";
				for(my $j=0;$j<@out;$j++){
					push @{$ra_out},$out[$j] unless $data_source->f_hide_columns($j);
				}
				$n++;
				if(defined($data_source->key)&&$fh->eof()){
					seek $fh, $pos, 0 or die "Couldn't seek to $pos: $!\n" if defined $pos;
				}
				last;
			}else{
				if(defined $data_source->key){
					$pos=tell $fh;
				}
			}
		}
	}
	$T=undef if $n<@{$ra_fh};
	${$rs_line_index}++ if defined $T;
	return $T;
}

sub check_filter{
	my ($value,$operation,$threshold)=@_;
	my $t=0;
	if($operation =~/^>/){
		$t=($value>$threshold);
	}elsif($operation =~/^</){
		$t=($value<$threshold);
	}
	if($operation eq "!="){
		$t=($value!=$threshold);
	}elsif($operation =~/=/){
		$t=$t||($value==$threshold);
	}
	return $t;
}

sub check_row_conditions{
	my ($ra_line,$ra_indices,$ra_thresholds)=@_;
	return 1 unless defined $ra_thresholds;
	my @columns;
	extract_columns($ra_line,$ra_indices,\@columns) or die "\nIndex out of range!";
	my $T=1;
	for(my $i=0;$i<@columns;$i++){
		if(defined $ra_thresholds->[$i]){
			my $op;
			my $th;
			if($ra_thresholds->[$i]=~/^([<>]=?|!=|==)\s*(\S+)/){
				$op=$1;
				$th=$2;
			}else{
				die "\nUnable to parse selection condition: $ra_thresholds->[$i]!";
			}
			$T&&=check_filter($columns[$i],$op,$th);
		}
	}
	return $T;
}

if(defined $flist){
	my @data_sources;
	my @data;
	my %data_sources_names;
	my @fh_data_sources;
	my $prj_path;
	my $pairs_fn;
	my $N=6;
	open INPF, "<$flist" or die "\nUnable to open input file: $flist!";
	while(<INPF>){
		$_=$` if(/#/);
		if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
			my $name=$1;
			if($name eq "Path"){
				$prj_path=$2;
				$prj_path.="/" unless $prj_path=~/\/$/;
			}elsif($name eq "Pairs"){
				$pairs_fn=$2;
			}
		}else{
			my @line=split '\t',$_,-1;
			die "\nError in format of the project description file!" unless @line==$N;
			for(my $i=0;$i<$N;$i++){
				$line[$i]=~s/^\s+//;
				$line[$i]=~s/\s+$//;
			}
			my $name=$line[FILE_NAME];
			my $data_source;
			if(defined $data_sources_names{$name}){
				die "\nThe name of data source $name repeated in the $flist file!";
			}else{
				$data_source=DataSourceInfo->new();
				$data_source->name($name);
				$data_sources_names{$name}=1;
				my $n=0;
				$n=$line[SKIP_NLINES] if $line[SKIP_NLINES]=~/^\d+/;
				$data_source->skip_nlines($n);
				$data_source->key($line[KEY_COL_IDX]-1) if($line[KEY_COL_IDX] ne "");
			}
			my $ncol;
			for(my $i=COL_LIST;$i<$N;$i++){
				my @items=split(',',$line[$i]);
				if($i==COL_LIST){
					die "\nSelectors of columns weren't specified in $flist file!" unless @items;
					$ncol=@items;
					for(my $i=0;$i<$ncol;$i++){
						$items[$i]=~s/^\s+//;
						$items[$i]=~s/\s+$//;
						my $site;
						my $f_hide_col=0;
						my $filter_str;
						if($items[$i]=~/^(\d+)\s*(h(ide)*\s*)?(.*)/i){
							$site=$1;
							die "\nWrong column index $items[$i] in the source $name in file $flist!" if $site<1;
							$site--;
							$f_hide_col=1 if $2 ne "";
							$data_source->col_thresholds->[$i]=$4 if $4 ne "";
						}else{
							die "\nNo site id in selector $items[$i] in $flist!";
						}
						push @{$data_source->col_indices},$site;
						push @{$data_source->f_hide_columns},$f_hide_col;
					}
				}elsif($i==COL_NAME_LIST){
					if(@items){
						die "\nNumber of specified names of columns doesn't equal to the number of indices!" unless @items==$ncol;
						@{$data_source->col_names}=@items;
					}
				}elsif($i==STUB_VAL_LIST){
					if(@items){
						die "\nNumber of specified stubs doesn't equal to the number of indices!" unless @items==$ncol;
						@{$data_source->stub_values}=@items;
					}
				}
			}
			if($data_source->key){
				die "\nUndefined stub values for indexed data file!" unless @{$data_source->stub_values};
			}
			push @data_sources,$data_source;
		}
	}
	close INPF;
	die "\nThe file with pairs of associations wasn't specified in $flist!" unless defined $pairs_fn;
	my $header=DataRow->new();
	for(my $i=0;$i<@data_sources;$i++){
		my $fname=$prj_path.$data_sources[$i]->name;
		my $ncol=@{$data_sources[$i]->col_indices};
		my $fh = FileHandle->new();
		$fh->open("< $fname") or die "\nUnable to open input file $fname!";
		my $n=$data_sources[$i]->skip_nlines;
		#skip top lines if required
		while($n--){
			$_=<$fh>;
		}
		unless(@{$data_sources[$i]->col_names}){
			#extract file names from the first line
			while(<$fh>){
				chomp;
				if(/\S+/){
					my @line=split '\t';
					warn "\nIn the source $fname from $flist first row is interpreted as names of columns!";
					extract_columns(\@line,$data_sources[$i]->col_indices,$data_sources[$i]->col_names)
						or die "\nIndex out of range in source $fname!";
					last;
				}
			}
		}
		for(my $j=0;$j<@{$data_sources[$i]->col_names};$j++){
			push @{$header->data},$data_sources[$i]->col_names->[$j] unless $data_sources[$i]->f_hide_columns->[$j];
		}
		push @fh_data_sources,$fh;
	}
	$add_data{HEADER}=$header;
	my $I=0;
	my $T=1;
	while($T){
		#reading all data files
		my @line;
		my $t=get_line_from_sources(\@fh_data_sources,\@data_sources,\$I,\@line);
		if($t){
			my $row=DataRow->new();
			my $idx=$I-1;
			$row->index($idx);
			@{$row->data}=@line;
			$add_data{$idx}=$row;
		}
		last unless (defined $t);
		foreach my $fh(@fh_data_sources){
			if($fh->eof()){
				$T=0;
				last;
			}
		}
	}
	foreach my $fh(@fh_data_sources){
		$fh->close();
	}
	my $n=$I;
	$I=0;
	#reading a file with pairs
	open INPF, "<$pairs_fn" or die "\nUnable to open input file $pairs_fn!";
	while(<INPF>){
		chomp;
		s/^\s+//;
		s/\s+$//;
		if(/^(\d+\t\d+)\t/){
			$pair2idx{$1}=$I++;
		}
	}
	close INPF;
	die "\nNumber of rows in $pairs_fn equals to $I which is inconsistent with number of rows $n extracted from the record in additional data $flist!"
		unless ($n==$I);
}

open INPF,"<$ARGV[0]" or die "\nUnable to open input file: $ARGV[0]!";
my $ncol=0;
my $I=$skip_nlines;
while(<INPF>){
	chomp;
	if($I){
		$I--;
		next;
	}
	s/^\s+//;
	s/\s+$//;
	if(/\S+/){
		my @line=split '\t';
		if(@matrix==0){
			$ncol=@line;
		}else{
			die "\nWrong format of the matrix file: $ARGV[0]!" unless /^\d+\t\d+/;
			next if $line[0]==$line[1];
			next unless $f_ord_pairs||$line[0]<$line[1];
		}
		if(@line==$ncol){
			my $t=1;
			if(@matrix){
				foreach my $f(@mtx_filters){
					my $idx=$f->[0];
					die "\nColumn index $idx in $ARGV[0] is out of range!" unless $idx<$ncol;
					$t&&=check_filter($line[$idx],$f->[1],$f->[2]);
				}
			}
			if($t&&defined($flist)){
				my $idx;
				if(@matrix){
					$idx=$pair2idx{"$line[0]\t$line[1]"};
				}else{
					$idx="HEADER";
				}
				my $row;
				if(defined $idx){
					$row=$add_data{$idx};
				}else{
					warn "\nMatrix file $ARGV[0] contains a pair $line[0] $line[1] that is absent in the pairs file!";
				}
				if(defined $row){
					push @line,@{$row->data};
				}else{
					$t=0;
				}
			}
			push @matrix,[@line] if $t;
		}else{
			die "\nUnequal number of colums lines!";
		}
	}
}
close INPF;
my @header_line=@{shift @matrix};
$ncol=@header_line;
for(my $i=0;$i<$ncol;$i++){ 
	my $str=$header_line[$i];
	print $str;
	print "\t" if $i<$ncol-1;
}
for(my $j=0;$j<@matrix;$j++){
	print "\n";
	for(my $i=0;$i<$ncol;$i++){ 
		my $str=$matrix[$j]->[$i];
		print $str;
		print "\t" if $i<$ncol-1;
	}
}
