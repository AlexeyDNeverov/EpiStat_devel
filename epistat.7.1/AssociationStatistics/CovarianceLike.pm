package AssociationStatistics::CovarianceLike;
#Module for calculating covariance/correlation-like association statistics for a site pair on the base of epistatic statistics
use AssociationStatistics::BaseAssociationMeasure;
@ISA = ("AssociationStatistics::BaseAssociationMeasure");

sub _init{
	my $self=shift;
	my ($dataset_desc,$general_settings)=@_;
	$self->SUPER::_init($dataset_desc,$general_settings,2);
}

sub _init_copy{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in CovarianceLike::_init_copy(): copy constructor is disabled!";
}

sub _init_transp{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in CovarianceLike::_init_transp(): matrix transposition is disabled!";
}

#interface declaration

sub mtx_diag_value{
	my $self=shift;
	return 1.0;
}

sub get_statistics{
	my $self=shift;
	my @stat=$self->SUPER::get_statistics;
	for(my $i=0;$i<@stat;$i++){
		my $k=$self->_line2fgr_idx($i);
		my $m=$self->{FGR_MUT_NUMBERS}->[$k]-$self->{MEAN}->[$i];
		$m=$self->{MEAN}->[$i] if $m<$self->{MEAN}->[$i];
		$stat[$i]/=$m;
	}
	return @stat;
}

1;
	