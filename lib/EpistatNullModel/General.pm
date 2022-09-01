package EpistatNullModel::General;
#This module provides functions for calculating statistics on set of trees with randomized mutation distributions
use strict;
use SitePairMatrix;
use IO::SitePairIdxData;
use Class::Struct;
use IndexedFakeSamplesHeap;
use Time::Progress;
#use Scalar::Util qw(blessed);
#interface declaration
#constructor
sub new{
	my $class=shift;
	my $self={};
	bless($self,$class);
	$self->_init(@_);
	return $self;
}

struct SitePairStat =>{
	summ => '$',
	summSQ => '$',
	lower_count => '$',
	upper_count => '$',
	unordered_lower_count => '$',
	unordered_upper_count => '$',
	summXY => '$'
};

struct SiteStat =>{
	summ => '$',
	summSQ => '$'
};

sub _init{
	my $self=shift;
	my %args = @_;
	my $nperm=0;
	my ($samples_dir,$stats_ext,$pairs_ext,$stats_fn,$pairs_fn,$f_intragene);
	$self->{SKIP_ID}=undef;
	my $max_samples_indir=10000000;
	my $samples_heap;
	my $f_square_mtx;
	my $bgr_site_num;
	my $fgr_site_num;
	my @line2site_indices;
	my @line2compl_line;
	my @summ;
	my @lower_count;
	my @upper_count;
	my @unordered_lower_count;
	my @unordered_upper_count;
	my @bgr_summ;
	my @fgr_summ;
	while(my ($k,$v)=each %args){
		$v=~s/^\s*//;
		$v=~s/^\s*$//;
		if($k eq "-size"){
			$nperm=$v if $v=~m/^\d+/;
		}elsif($k eq "-samples_dir"){
			$v=~s/\/$//;
			$samples_dir=$v;
		}elsif($k eq "-max_samples_indir"){
			$max_samples_indir=$v;
		}elsif($k eq "-stats_ext"){
			$stats_ext=$v;
		}elsif($k eq "-pairs_ext"){
			$pairs_ext=$v;
		}elsif($k eq "-obs_stats_fn"){
			$stats_fn=$v;
		}elsif($k eq "-obs_pairs_fn"){
			$pairs_fn=$v;
		}elsif($k eq "-is_intragene"){
			$f_intragene=$v;
		}elsif($k eq "-skip_obs_id"){
			$self->{SKIP_ID}=$v;
		}else{
			die "\nUnknown parameter: $k!";
		};
	}
	my @samples_storage;
	if($nperm>$max_samples_indir){
		$samples_heap=IndexedFakeSamplesHeap->new($nperm,$max_samples_indir);
		my @tmp=$samples_heap->get_storage_path_list;
		for(my $i=0;$i<@tmp-1;$i++){
			push @samples_storage,[($max_samples_indir*$i+1,$max_samples_indir*($i+1),"/".$tmp[$i])];
		}
		push @samples_storage,[($max_samples_indir*(@tmp-1)+1,$nperm,"/".$tmp[-1])];
	}else{
		push @samples_storage,[(1,$nperm,"")];
	}
	die "\nParameter -stats_ext has no default value!" unless defined $stats_ext;
	if(defined $self->{SKIP_ID}){
		my $indir=$samples_dir;
		if(defined $samples_heap){
			$indir.="/".$samples_heap->sample_idx2path($self->{SKIP_ID});
		}
		$stats_fn=$indir."/".$self->{SKIP_ID}.$stats_ext unless defined($stats_fn);
		$pairs_fn=$indir."/".$self->{SKIP_ID}.$pairs_ext if(defined($pairs_ext)&&(!defined($pairs_fn)));
	}
	$self->{SAMPLE_SIZE}=$nperm;
	$nperm-- if defined $self->{SKIP_ID};
	die "\nError package EpistatNullModel::_init(): No permutations sampled!" unless $nperm>0;
	$self->{OBS_STAT}={};
	$self->{SAMPLE_STAT}=[];
	$self->{SITE_PAIR_MTX}=undef;
	$self->{BGR_SITE_STAT}=undef;
	$self->{FGR_SITE_STAT}=undef;
	#########################
	my @tmp_data=($self->{OBS_STAT});
	my @col_idx=(1);
	my $nlines=IO::SitePairIdxData::hread_file($stats_fn,5,\@tmp_data,\@col_idx);
	#########################
	if(defined($pairs_fn)&&defined($f_intragene)){
		$self->{BGR_SITE_STAT}=[];
		$self->{FGR_SITE_STAT}=[];
		$self->{SITE_PAIR_MTX}=SitePairMatrix->new();
		$self->{SITE_PAIR_MTX}->init_sites($pairs_fn,$f_intragene);
		$bgr_site_num=$self->{SITE_PAIR_MTX}->get_bgr_site_num;
		$fgr_site_num=$self->{SITE_PAIR_MTX}->get_fgr_site_num;
		for(my $i=0;$i<$bgr_site_num;$i++){
			push @{$self->{BGR_SITE_STAT}},SiteStat->new();
		}
		for(my $i=0;$i<$fgr_site_num;$i++){
			push @{$self->{FGR_SITE_STAT}},SiteStat->new();
		}
		$f_square_mtx=$self->{SITE_PAIR_MTX}->is_square;
		for(my $j=0;$j<$nlines;$j++){
			my ($bgr_idx,$fgr_idx)=$self->{SITE_PAIR_MTX}->line2sites_idxs($j);
			push @line2site_indices,[($bgr_idx,$fgr_idx)];
			if($f_square_mtx){
				my $symj=$self->{SITE_PAIR_MTX}->site_idx_pair2line($fgr_idx,$bgr_idx);
				push @line2compl_line,$symj;
			}
		}
	}
	print STDERR "\nEpistatNullModel::init_(): Start reading data files:";
	my $p = Time::Progress->new(min => 1, max => $self->{SAMPLE_SIZE});
	print STDERR "\nFirst pass:\n";
	my @stub_values=(0);
	for(my $si=0;$si<@samples_storage;$si++){
		my $from=$samples_storage[$si]->[0];
		my $to=$samples_storage[$si]->[1];
		my $folder=$samples_storage[$si]->[2];
		for(my $i=$from;$i<=$to;$i++){
			print STDERR $p->report("\r%20b  ETA: %E", $i);
			next if $i==$self->{SKIP_ID};
			my $str=$samples_dir.$folder."/".$i.$stats_ext;
			my @stat_buff;
			@tmp_data=(\@stat_buff);
			IO::SitePairIdxData::aread_file($str,5,\@tmp_data,\@stub_values,\@col_idx);
			die "\nError in EpistatNullModel::_init(): Unexpected number of records ".scalar(@stat_buff)." in the file $str!" unless @stat_buff==$nlines;
			my @bgr_site_stat;
			my @fgr_site_stat;
			if(defined($self->{SITE_PAIR_MTX})){
				for(my $i=0;$i<$bgr_site_num;$i++){
					push @bgr_site_stat,0;
				}
				for(my $i=0;$i<$fgr_site_num;$i++){
					push @fgr_site_stat,0;
				}
			}
			for(my $j=0;$j<$nlines;$j++){
				my $obs_stat=0;
				if(exists $self->{OBS_STAT}->{$j}){
					$obs_stat=$self->{OBS_STAT}->{$j};
				}
				if($obs_stat<=$stat_buff[$j]){
					$upper_count[$j]++;
				}
				if($obs_stat>=$stat_buff[$j]){
					$lower_count[$j]++;
				}
				if($stat_buff[$j]>0){
					$summ[$j]+=$stat_buff[$j];
				}
				if(defined($self->{SITE_PAIR_MTX})){
					my ($bgr_idx,$fgr_idx)=@{$line2site_indices[$j]};
					$bgr_site_stat[$bgr_idx]+=$stat_buff[$j];
					$fgr_site_stat[$fgr_idx]+=$stat_buff[$j];
					if($f_square_mtx){
						my $symj=$line2compl_line[$j];
						$obs_stat+=$self->{OBS_STAT}->{$symj} if(exists $self->{OBS_STAT}->{$symj});
						if($obs_stat<=$stat_buff[$j]+$stat_buff[$symj]){
							$unordered_upper_count[$j]++;
						}
						if($obs_stat>=$stat_buff[$j]+$stat_buff[$symj]){
							$unordered_lower_count[$j]++;
						}
					}
				}
			}
			if(defined($self->{SITE_PAIR_MTX})){
				for(my $i=0;$i<$bgr_site_num;$i++){
					$bgr_summ[$i]+=$bgr_site_stat[$i];
				}
				for(my $i=0;$i<$fgr_site_num;$i++){
					$fgr_summ[$i]+=$fgr_site_stat[$i];
				}
			}
		}
	}
	for(my $j=0;$j<$nlines;$j++){
		my $sps=SitePairStat->new();
		$sps->summ($summ[$j]);
		$sps->upper_count($upper_count[$j]);
		$sps->lower_count($lower_count[$j]);
		if($f_square_mtx){
			$sps->unordered_upper_count($unordered_upper_count[$j]);
			$sps->unordered_lower_count($unordered_lower_count[$j]);
		}
		push @{$self->{SAMPLE_STAT}},$sps;
	}
	if(defined($self->{SITE_PAIR_MTX})){
		for(my $i=0;$i<$bgr_site_num;$i++){
			my $sps=$self->{BGR_SITE_STAT}->[$i];
			$sps->summ($bgr_summ[$i]);
		}
		for(my $i=0;$i<$fgr_site_num;$i++){
			my $sps=$self->{FGR_SITE_STAT}->[$i];
			$sps->summ($fgr_summ[$i]);
		}
	}
	unless(defined $self->{SKIP_ID}){
		my @summSQ;
		my @summXY;
		my @bgr_summSQ;
		my @fgr_summSQ;
		print STDERR "\nEpistatNullModel::init_(): Start reading data files:";
		$p = Time::Progress->new(min => 1, max => $self->{SAMPLE_SIZE});
		print STDERR "\nSecond pass:\n";
		for(my $si=0;$si<@samples_storage;$si++){
			my $from=$samples_storage[$si]->[0];
			my $to=$samples_storage[$si]->[1];
			my $folder=$samples_storage[$si]->[2];
			for(my $i=$from;$i<=$to;$i++){
				print STDERR $p->report("\r%20b  ETA: %E", $i);
				next if $i==$self->{SKIP_ID};
				my $str=$samples_dir.$folder."/".$i.$stats_ext;
				my @stat_buff;
				@tmp_data=(\@stat_buff);
				IO::SitePairIdxData::aread_file($str,5,\@tmp_data,\@stub_values,\@col_idx);
				die "\nError in EpistatNullModel::_init(): Unexpected number of records in the file $str!" unless @stat_buff==$nlines;
				my @bgr_site_stat;
				my @fgr_site_stat;
				if(defined($self->{SITE_PAIR_MTX})){
					for(my $i=0;$i<$bgr_site_num;$i++){
						push @bgr_site_stat,0;
					}
					for(my $i=0;$i<$fgr_site_num;$i++){
						push @fgr_site_stat,0;
					}
				}
				for(my $j=0;$j<$nlines;$j++){
					my $delta=($stat_buff[$j]-$summ[$j]/$nperm);
					$summSQ[$j]+=$delta**2;
					if(defined($self->{SITE_PAIR_MTX})){
						my ($bgr_idx,$fgr_idx)=@{$line2site_indices[$j]};
						$bgr_site_stat[$bgr_idx]+=$stat_buff[$j];
						$fgr_site_stat[$fgr_idx]+=$stat_buff[$j];
						if($f_square_mtx){
							my $symj=$line2compl_line[$j];
							$summXY[$j]+=$delta*($stat_buff[$symj]-$summ[$symj]/$nperm);
						}
					}
				}
				if(defined($self->{SITE_PAIR_MTX})){
					for(my $i=0;$i<$bgr_site_num;$i++){
						$bgr_summSQ[$i]+=($bgr_site_stat[$i]-$bgr_summ[$i]/$nperm)**2;
					}
					for(my $i=0;$i<$fgr_site_num;$i++){
						$fgr_summSQ[$i]+=($fgr_site_stat[$i]-$fgr_summ[$i]/$nperm)**2;
					}
				}
			}
		}
		for(my $j=0;$j<$nlines;$j++){
			my $sps=$self->{SAMPLE_STAT}->[$j];
			$sps->summSQ($summSQ[$j]);
			$sps->summXY($summXY[$j]) if($f_square_mtx);
		}
		if(defined($self->{SITE_PAIR_MTX})){
			for(my $i=0;$i<$bgr_site_num;$i++){
				my $sps=$self->{BGR_SITE_STAT}->[$i];
				$sps->summSQ($bgr_summSQ[$i]);
			}
			for(my $i=0;$i<$fgr_site_num;$i++){
				my $sps=$self->{FGR_SITE_STAT}->[$i];
				$sps->summSQ($fgr_summSQ[$i]);
			}
		}
	}
	print STDERR "\nEpistatNullModel::init_(): done!";
}

sub is_square{
	my $self=shift;
	return undef unless defined $self->{SITE_PAIR_MTX};
	return $self->{SITE_PAIR_MTX}->is_square;
}
	
sub calc_pvalues{
	my $self=shift;
	my ($ra_lower_pval,$ra_upper_pval)=@_;
	die "\nError: undefined value not allowed for lower pvalue ARRAYREF!" unless defined $ra_lower_pval;
	die "\nError: undefined value not allowed for upper pvalue ARRAYREF!" unless defined $ra_upper_pval;
	@{$ra_lower_pval}=();
	@{$ra_upper_pval}=();
	my $nperm=$self->{SAMPLE_SIZE};
	$nperm-- if defined $self->{SKIP_ID};
	for(my $idx=0;$idx<@{$self->{SAMPLE_STAT}};$idx++){
		$ra_lower_pval->[$idx]=$self->{SAMPLE_STAT}->[$idx]->lower_count/$nperm;
		$ra_upper_pval->[$idx]=$self->{SAMPLE_STAT}->[$idx]->upper_count/$nperm;
	}
}

sub calc_unordered_pairs_pvalues{
	my $self=shift;
	my ($ra_lower_pval,$ra_upper_pval)=@_;
	die "\nError: undefined value not allowed for lower pvalue ARRAYREF!" unless defined $ra_lower_pval;
	die "\nError: undefined value not allowed for upper pvalue ARRAYREF!" unless defined $ra_upper_pval;
	@{$ra_lower_pval}=();
	@{$ra_upper_pval}=();
	if($self->is_square){
		my $nperm=$self->{SAMPLE_SIZE};
		$nperm-- if defined $self->{SKIP_ID};
		for(my $idx=0;$idx<@{$self->{SAMPLE_STAT}};$idx++){
			$ra_lower_pval->[$idx]=$self->{SAMPLE_STAT}->[$idx]->unordered_lower_count/$nperm;
			$ra_upper_pval->[$idx]=$self->{SAMPLE_STAT}->[$idx]->unordered_upper_count/$nperm;
		}
		return 1;
	}
	return 0;
}

sub calc_moments12{
	my $self=shift;
	my ($ra_avg,$ra_var)=@_;
	die "\nNot applicable for fake data: ".$self->{SKIP_ID} if defined $self->{SKIP_ID};
	die "\nError: undefined value not allowed for average ARRAYREF!" unless defined $ra_avg;
	die "\nError: undefined value not allowed for variance ARRAYREF!" unless defined $ra_var;
	@{$ra_avg}=();
	@{$ra_var}=();
	my $nperm=$self->{SAMPLE_SIZE};
	$nperm-- if defined $self->{SKIP_ID};
	for(my $idx=0;$idx<@{$self->{SAMPLE_STAT}};$idx++){
		$ra_avg->[$idx]=$self->{SAMPLE_STAT}->[$idx]->summ/$nperm;
		$ra_var->[$idx]=$self->{SAMPLE_STAT}->[$idx]->summSQ;
		$ra_var->[$idx]/=($nperm-1) if($nperm>1);
	}
}

sub calc_ordered_site_pairs_covariances{
	my $self=shift;
	my ($ra_cov)=@_;
	die "\nNot applicable for fake data: ".$self->{SKIP_ID} if defined $self->{SKIP_ID};
	@{$ra_cov}=();
	if($self->is_square){
		my $nperm=$self->{SAMPLE_SIZE};
		$nperm-- if defined $self->{SKIP_ID};
		for(my $idx=0;$idx<@{$self->{SAMPLE_STAT}};$idx++){
			$ra_cov->[$idx]=$self->{SAMPLE_STAT}->[$idx]->summXY;
			$ra_cov->[$idx]/=($nperm-1) if($nperm>1);
		}
		return 1;
	}
	return 0;
}

sub calc_bgr_moments12{
	my $self=shift;
	my ($ra_avg,$ra_var)=@_;
	die "\nNot applicable for fake data: ".$self->{SKIP_ID} if defined $self->{SKIP_ID};
	die "\nError: undefined value not allowed for average ARRAYREF!" unless defined $ra_avg;
	die "\nError: undefined value not allowed for variance ARRAYREF!" unless defined $ra_var;
	@{$ra_avg}=();
	@{$ra_var}=();
	return unless defined $self->{BGR_SITE_STAT};
	my $nperm=$self->{SAMPLE_SIZE};
	$nperm-- if defined $self->{SKIP_ID};
	for(my $idx=0;$idx<@{$self->{BGR_SITE_STAT}};$idx++){
		$ra_avg->[$idx]=$self->{BGR_SITE_STAT}->[$idx]->summ/$nperm;
		$ra_var->[$idx]=$self->{BGR_SITE_STAT}->[$idx]->summSQ;
		$ra_var->[$idx]/=($nperm-1) if($nperm>1);
	}
}

sub calc_fgr_moments12{
	my $self=shift;
	my ($ra_avg,$ra_var)=@_;
	die "\nNot applicable for fake data: ".$self->{SKIP_ID} if defined $self->{SKIP_ID};
	die "\nError: undefined value not allowed for average ARRAYREF!" unless defined $ra_avg;
	die "\nError: undefined value not allowed for variance ARRAYREF!" unless defined $ra_var;
	@{$ra_avg}=();
	@{$ra_var}=();
	return unless defined $self->{FGR_SITE_STAT};
	my $nperm=$self->{SAMPLE_SIZE};
	$nperm-- if defined $self->{SKIP_ID};
	for(my $idx=0;$idx<@{$self->{FGR_SITE_STAT}};$idx++){
		$ra_avg->[$idx]=$self->{FGR_SITE_STAT}->[$idx]->summ/$nperm;
		$ra_var->[$idx]=$self->{FGR_SITE_STAT}->[$idx]->summSQ;
		$ra_var->[$idx]/=($nperm-1) if($nperm>1);
	}
}

sub get_sites{
	my $self=shift;
	my ($ra_bgr_sites,$ra_fgr_sites)=@_;
	if(defined $self->{SITE_PAIR_MTX}){
		$self->{SITE_PAIR_MTX}->get_sites($ra_bgr_sites,$ra_fgr_sites);
		return 1;
	}
	return 0;
}

sub DESTROY {
	my $self = shift;
	undef($self->{SAMPLE_STAT});
	undef($self->{OBS_STAT});
	if(defined $self->{BGR_SITE_STAT}){
		undef($self->{BGR_SITE_STAT});
	}
	if(defined $self->{FGR_SITE_STAT}){
		undef($self->{FGR_SITE_STAT});
	}
	undef($self->{SITE_PAIR_MTX});
	print STDERR "\nEpistatNullModel DESTROYED";
}
	
sub get_observation_stats{
	my $self=shift;
	my @obs_stat;
	for(my $idx=0;$idx<@{$self->{SAMPLE_STAT}};$idx++){
		my $obs_stat=0;
		if(exists $self->{OBS_STAT}->{$idx}){
			$obs_stat=$self->{OBS_STAT}->{$idx};
		}
		push @obs_stat,$obs_stat;
	}
	return @obs_stat;
};

1;