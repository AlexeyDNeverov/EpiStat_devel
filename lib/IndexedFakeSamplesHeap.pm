package IndexedFakeSamplesHeap;
#Emulates folder for storage of huge amount (billions) of samples from some process which generates files of the same format.
#The samples are required to be indexed starting from 1

use strict;
use Class::Struct;

struct StorageNodeInfo =>{
	id => '$',
	level => '$',
	parent => '$',
	children => '@'
};

sub new{
	my $class=shift;
	my $self={};
	bless($self,$class);
	$self->_init(@_);
	return $self;
}

sub sample_idx2path{
	my $self=shift;
	my ($file_idx)=@_;
	die "\nError IndexedFakeSamplesHeap::sample_idx2path: index must be an unsigned integer!" unless $file_idx=~/^[1-9]\d*$/;
	$file_idx--;
	my $storage_idx=int($file_idx/$self->{FOLDER_CAPACITY});
	my $storage_folder="";
	$storage_folder=$self->{PATHS}->[$storage_idx] if @{$self->{PATHS}}>1;
	return $storage_folder;
}

sub get_max_samples_indir{
	my $self=shift;
	return $self->{FOLDER_CAPACITY};
}

sub get_storage_path_list{
	my $self=shift;
	return @{$self->{PATHS}};
}

sub _init{
	my $self=shift;
	my ($nsamples,$max_items_inside_folder)=@_;
	die "\nError IndexedFakeSamplesHeap::_init: Number of samples must be an unsigned integer!" unless $nsamples=~/^[1-9]\d*$/;
	die "\nError IndexedFakeSamplesHeap::_init: Number of samples in a folder must be an unsigned integer!" unless $max_items_inside_folder=~/^[1-9]\d*$/;
	$self->{FOLDER_CAPACITY}=$max_items_inside_folder;
	$self->{NSAMPLES}=$nsamples;
	$self->{PATHS}=[];
	my @folders0=_make_storage_tree($nsamples,$max_items_inside_folder);
	if(@folders0>1){
		foreach my $fi(@folders0){
			my $path="L".$fi->level."I".$fi->id;
			my $pfi=$fi->parent;
			while(defined $pfi->parent){
				my $str="L".$pfi->level."I".$pfi->id."/";
				$path=$str.$path;
				$pfi=$fi->parent;
			}
			push @{$self->{PATHS}},$path;
		}
	}
}

sub _make_storage_tree{
	my ($nsamples,$max_items_inside)=@_;
	my @folders;
	my $n=int(($nsamples-1)/$max_items_inside);
	$n++;
	my $level=0;
	for(my $i=0;$i<$n;$i++){
		my $fi=StorageNodeInfo->new();
		$fi->id($i);
		$fi->level($level);
		push @folders,$fi;
	}
	my @folders0=@folders;
	while($n-1){
		my $n1=int(($n-1)/$max_items_inside);
		$n1++;
		$level++;
		my @folders1;
		for(my $i=0;$i<$n1;$i++){
			my $fi=StorageNodeInfo->new();
			$fi->id($i);
			$fi->level($level);
			push @folders1,$fi;
		}
		for(my $j=0;$j<$n;$j++){
			my $i=int($j/$max_items_inside);
			$folders[$j]->parent($folders1[$i]);
			push @{$folders1[$i]->children},$folders[$j];
		}
		$n=$n1;
		@folders=@folders1;
	}
	if(wantarray()){
		return @folders0;
	}
	return $folders[0];
}

1;