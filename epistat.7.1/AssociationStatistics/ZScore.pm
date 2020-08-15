package AssociationStatistics::ZScore;
#Module for calculating covariation/correlation-like association statistics for a site pair on the base of epistatic statistics
use AssociationStatistics::BaseAssociationMeasure;
@ISA = ("AssociationStatistics::BaseAssociationMeasure");

sub _init{
	my $self=shift;
	my ($dataset_desc,$rh_site_pair_filters)=@_;
	$self->SUPER::_init($dataset_desc,$rh_site_pair_filters);
	#check sufficiency of the dataset description
	die "\nError in ZScore::_init(): undefined file with pairs of sites!" unless defined $dataset_desc->variance_fn;
	###
	$self->{VARIANCE}=[];
	my $variance_fn=$dataset_desc->variance_fn;
	open INPF, "<$variance_fn" or die "\nUnable to open input file $variance_fn!";
	while(<INPF>){
		chomp;
		$_=~s/\s*$//;
		if($_ ne ""){
			push @{$self->{VARIANCE}},$_;
		}
	}
	close INPF;
	die "\nNumber of lines in the file $variance_fn isn't equal to number of lines in files with epistatic statistics!" unless @{$self->{VARIANCE}}==$self->{NLINES};
}

sub _init_copy{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in ZScore::_init_copy(): copy constructor is disabled!";
}

sub _init_transp{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in ZScore::_init_transp(): undefined matrix argument!" unless defined $sp_matrix;
	$self->SUPER::_init_transp($sp_matrix,@_);
	$self->{VARIANCE}=[];
	for(my $i=0;$i<$self->{NLINES};$i++){
		my ($bg_idx,$fg_idx)=$self->line2sites_idxs($i);
		my $ti=$sp_matrix->site_idx_pair2line($fg_idx,$bg_idx);
		$self->{VARIANCE}->[$i]=$sp_matrix->{VARIANCE}->[$ti];
	}
}

#interface declaration

sub mtx_diag_value{
	my $self=shift;
	return 1.0;
}

#$norm_const==undef - default normalization
#$norm_const==0 - no normalization
#$norm_const!=0 - the value is used for normalization
sub get_statistics{
	my $self=shift;
	my ($norm_const)=@_;
	my @stat=$self->SUPER::get_statistics;
	$M=0;
	for(my $i=0;$i<@stat;$i++){
		if(defined($self->{VARIANCE}->[$i])&&($self->{VARIANCE}->[$i]>0)){
			$stat[$i]/=sqrt($self->{VARIANCE}->[$i]);
			$M=abs($stat[$i]) if(abs($stat[$i])>$M);
		}else{
			$stat[$i]=undef;
		}
	}
	$M=$norm_const if defined $norm_const;
	if($M){
		for(my $i=0;$i<@stat;$i++){
			$stat[$i]/=$M if defined $stat[$i];
		}
	}
	return @stat;
}

1;
	