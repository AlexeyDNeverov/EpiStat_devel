package EpistatNullModel::Symmetric;
@ISA = ("EpistatNullModel::General");
#This module provides functions for calculating statistics on set of trees with randomized mutation distributions
use strict;
use SitePairMatrix;
use IO::SitePairIdxData;
use Class::Struct;
use IndexedFakeSamplesHeap;
use Time::Progress; 
#interface declaration
#constructor
sub new{
	my $class=shift;
	my $self={};
	bless($self,$class);
	$self->_init(@_);
	return $self;
}

#struct SitePairStat =>{
#	summ => '$',
#	summSQ => '$',
#	lower_count => '$',
#	upper_count => '$',
#};
#
#struct SiteStat =>{
#	summ => '$',
#	summSQ => '$'
#};

sub _init{
	my $self=shift;
	my %args = @_;
	my $nperm=0;
	my ($samples_dir,$stats_ext,$pairs_ext,$stats_fn,$pairs_fn);
	$self->{SKIP_ID}=undef;
	my $max_samples_indir=10000000;
	my $samples_heap;
	my $f_intragene=1; #!!!
	my $site_num;
	my @line2site_indices;
	my @line2compl_line;
	my @summ;
	my @lower_count;
	my @upper_count;
	my @site_summ;
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
	$self->{SITE_STAT}=[];
	$self->{IN2OUT_LINE}=[];
	#########################
	$self->{SITE_PAIR_MTX}=SitePairMatrix->new();
	$self->{SITE_PAIR_MTX}->init_sites($pairs_fn,$f_intragene);
	$site_num=$self->{SITE_PAIR_MTX}->get_bgr_site_num;
	die "\nThe Symmetric epistat null model is applicable only for symmetric statistics!" unless $site_num==$self->{SITE_PAIR_MTX}->get_fgr_site_num;
	for(my $i=0;$i<$site_num;$i++){
		push @{$self->{SITE_STAT}},SiteStat->new();
	}
	for(my $i=0;$i<$site_num;$i++){
		for(my $j=$i+1;$j<$site_num;$j++){
			push @line2site_indices,[($i,$j)];
			push @{$self->{IN2OUT_LINE}},$self->{SITE_PAIR_MTX}->site_idx_pair2line($i,$j);
		}
	}
	##########################
	my $nlines;
	my @tmp_data;
	my @tmp_inds;
	my @col_idx=(1);
	{
		my @stat_buff;
		my @stat_inds;
		@tmp_data=(\@stat_buff);
		$nlines=IO::SitePairIdxData::airead_file($stats_fn,5,\@stat_inds,\@tmp_data,\@col_idx);
		die "\nThe wrong number of index lines in $stats_fn!" unless $nlines==$site_num*($site_num-1)/2;
		my $r_ind=-1;
		my $r_diff=0;
		for(my $i=0;$i<@stat_inds;$i++){
			my $j=$stat_inds[$i];
			my $obs_stat=$stat_buff[$i];
			#convertion of coordinates in the full matrix array into coordinates of upper-triangle dioganal matrix array
			my $r=int($j/($site_num-1));
			unless($r==$r_ind){
				$r_ind=$r;
				$r_diff=$r_ind*($r_ind+1)/2;
			}
			$j-=$r_diff;
			$self->{OBS_STAT}->{$j}=$obs_stat;
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
			my @stat_inds;
			@tmp_data=(\@stat_buff);
			IO::SitePairIdxData::airead_file($str,5,\@stat_inds,\@tmp_data,\@col_idx);
			my @site_stat;
			for(my $i=0;$i<$site_num;$i++){
				push @site_stat,0;
			}
			my $r_ind=-1;
			my $r_diff=0;
			my $from=0;
			for(my $k=0;$k<@stat_inds;$k++){
				my $j=$stat_inds[$k];
				my $r=int($j/($site_num-1));
				unless($r==$r_ind){
					$r_ind=$r;
					$r_diff=$r_ind*($r_ind+1)/2;
				}
				$j-=$r_diff;
				my $obs_stat=0;
				for(my $i=$from;$i<$j;$i++){
					$obs_stat=0;
					if(exists $self->{OBS_STAT}->{$i}){
						$obs_stat=$self->{OBS_STAT}->{$i};
					}
					if($obs_stat<=0){
						$upper_count[$i]++;
					}
					if($obs_stat>=0){
						$lower_count[$i]++;
					}
				}
				$from=$j+1;
				$obs_stat=0;
				if(exists $self->{OBS_STAT}->{$j}){
					$obs_stat=$self->{OBS_STAT}->{$j};
				}
				if($obs_stat<=$stat_buff[$k]){
					$upper_count[$j]++;
				}
				if($obs_stat>=$stat_buff[$k]){
					$lower_count[$j]++;
				}
				$summ[$j]+=$stat_buff[$k];
				my ($bgr_idx,$fgr_idx)=@{$line2site_indices[$j]};
				$site_stat[$bgr_idx]+=$stat_buff[$k];
				$site_stat[$fgr_idx]+=$stat_buff[$k];
			}
			for(my $i=$from;$i<$nlines;$i++){
				my $obs_stat=0;
				if(exists $self->{OBS_STAT}->{$i}){
					$obs_stat=$self->{OBS_STAT}->{$i};
				}
				if($obs_stat<=0){
					$upper_count[$i]++;
				}
				if($obs_stat>=0){
					$lower_count[$i]++;
				}
			}
			for(my $i=0;$i<$site_num;$i++){
				$site_summ[$i]+=$site_stat[$i];
			}
		}
	}
	for(my $j=0;$j<$nlines;$j++){
		my $sps=SitePairStat->new();
		$sps->summ($summ[$j]);
		$sps->upper_count($upper_count[$j]);
		$sps->lower_count($lower_count[$j]);
		push @{$self->{SAMPLE_STAT}},$sps;
	}
	for(my $i=0;$i<$site_num;$i++){
		my $sps=$self->{SITE_STAT}->[$i];
		$sps->summ($site_summ[$i]);
	}

	unless(defined $self->{SKIP_ID}){
		my @summSQ;
		my @summXY;
		my @site_summSQ;
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
				my @stat_inds;
				@tmp_data=(\@stat_buff);
				IO::SitePairIdxData::airead_file($str,5,\@stat_inds,\@tmp_data,\@col_idx);
				my @site_stat;
				for(my $i=0;$i<$site_num;$i++){
					push @site_stat,0;
				}
				my $r_ind=-1;
				my $r_diff=0;
				my $from=0;
				for(my $k=0;$k<@stat_inds;$k++){
					my $j=$stat_inds[$k];
					my $r=int($j/($site_num-1));
					unless($r==$r_ind){
						$r_ind=$r;
						$r_diff=$r_ind*($r_ind+1)/2;
					}
					$j-=$r_diff;
					my $delta=0;
					for(my $i=$from;$i<$j;$i++){
						$delta=-$summ[$i]/$nperm;
						$summSQ[$i]+=$delta**2;
					}
					$from=$j+1;
					$delta=($stat_buff[$k]-$summ[$j]/$nperm);
					$summSQ[$j]+=$delta**2;
					my ($bgr_idx,$fgr_idx)=@{$line2site_indices[$j]};
					$site_stat[$bgr_idx]+=$stat_buff[$k];
					$site_stat[$fgr_idx]+=$stat_buff[$k];
				}
				for(my $i=$from;$i<$nlines;$i++){
					my $delta=-$summ[$i]/$nperm;
					$summSQ[$i]+=$delta**2;
				}
				for(my $i=0;$i<$site_num;$i++){
					$site_summSQ[$i]+=($site_stat[$i]-$site_summ[$i]/$nperm)**2;
				}
			}
		}
		for(my $j=0;$j<$nlines;$j++){
			my $sps=$self->{SAMPLE_STAT}->[$j];
			$sps->summSQ($summSQ[$j]);
		}
		for(my $i=0;$i<$site_num;$i++){
			my $sps=$self->{SITE_STAT}->[$i];
			$sps->summSQ($site_summSQ[$i]);
		}
	}
	print STDERR "\nEpistatNullModel::init_(): done!";
}

sub is_square{
	my $self=shift;
	return 1;
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
	for(my $i=0;$i<@{$self->{SAMPLE_STAT}};$i++){
		my $idx=$self->{IN2OUT_LINE}->[$i];
		$ra_lower_pval->[$idx]=$self->{SAMPLE_STAT}->[$i]->lower_count/$nperm;
		$ra_upper_pval->[$idx]=$self->{SAMPLE_STAT}->[$i]->upper_count/$nperm;
		my ($a1,$a2)=$self->{SITE_PAIR_MTX}->line2sites_idxs($idx);
		my $sym_idx=$self->{SITE_PAIR_MTX}->site_idx_pair2line($a2,$a1);
		$ra_lower_pval->[$sym_idx]=$ra_lower_pval->[$idx];
		$ra_upper_pval->[$sym_idx]=$ra_upper_pval->[$idx];
	}
}

sub calc_unordered_pairs_pvalues{
	my $self=shift;
	return $self->calc_pvalues(@_);
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
	for(my $i=0;$i<@{$self->{SAMPLE_STAT}};$i++){
		my $idx=$self->{IN2OUT_LINE}->[$i];
		$ra_avg->[$idx]=$self->{SAMPLE_STAT}->[$i]->summ/$nperm;
		$ra_var->[$idx]=$self->{SAMPLE_STAT}->[$i]->summSQ;
		$ra_var->[$idx]/=($nperm-1) if($nperm>1);
		my ($a1,$a2)=$self->{SITE_PAIR_MTX}->line2sites_idxs($idx);
		my $sym_idx=$self->{SITE_PAIR_MTX}->site_idx_pair2line($a2,$a1);
		$ra_avg->[$sym_idx]=$ra_avg->[$idx];
		$ra_var->[$sym_idx]=$ra_var->[$idx];
	}
}

sub calc_ordered_site_pairs_covariances{
	my $self=shift;
	die "\nNot applicable: for symmetric matrices covariances between symmetric ordered pairs equal to variances!";
	return 0;
}

sub calc_site_moments12{
	my $self=shift;
	my ($ra_avg,$ra_var)=@_;
	die "\nNot applicable for fake data: ".$self->{SKIP_ID} if defined $self->{SKIP_ID};
	die "\nError: undefined value not allowed for average ARRAYREF!" unless defined $ra_avg;
	die "\nError: undefined value not allowed for variance ARRAYREF!" unless defined $ra_var;
	@{$ra_avg}=();
	@{$ra_var}=();
	my $nperm=$self->{SAMPLE_SIZE};
	$nperm-- if defined $self->{SKIP_ID};
	for(my $idx=0;$idx<@{$self->{SITE_STAT}};$idx++){
		$ra_avg->[$idx]=$self->{SITE_STAT}->[$idx]->summ/$nperm;
		$ra_var->[$idx]=$self->{SITE_STAT}->[$idx]->summSQ;
		$ra_var->[$idx]/=($nperm-1) if($nperm>1);
	}
}

sub calc_bgr_moments12{
	my $self=shift;
	die "\nNot applicable: use calc_site_moments12() instead!";
}

sub calc_fgr_moments12{
	my $self=shift;
	my $self=shift;
	die "\nNot applicable: use calc_site_moments12() instead!";
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
	undef($self->{SITE_STAT});
	undef($self->{SITE_PAIR_MTX});
	undef($self->{IN2OUT_LINE});
	print STDERR "\nEpistatNullModel DESTROYED";
}
	
sub get_observation_stats{
	my $self=shift;
	my @obs_stat;
	for(my $i=0;$i<@{$self->{SAMPLE_STAT}};$i++){
		my $idx=$self->{IN2OUT_LINE}->[$i];
		my $obs_stat=0;
		if(exists $self->{OBS_STAT}->{$i}){
			$obs_stat=$self->{OBS_STAT}->{$i};
		}
		$obs_stat[$idx]=$obs_stat;
		my ($a1,$a2)=$self->{SITE_PAIR_MTX}->line2sites_idxs($idx);
		my $sym_idx=$self->{SITE_PAIR_MTX}->site_idx_pair2line($a2,$a1);
		$obs_stat[$sym_idx]=$obs_stat;
	}
	return @obs_stat;
}

1;