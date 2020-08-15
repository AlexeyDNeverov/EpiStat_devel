package AssociationStatistics::BasicIndexFilter;
#Do not use this class directly!
use AssociationStatistics::basicAssociationFilter;
@ISA = ("AssociationStatistics::basicAssociationFilter");

sub _filter2flags{
	my $self=shift;
	my ($filter)=@_;
	my $n=$filter->{HOST_MATRIX}->{NLINES};
	my @tmp_flags=(1) x $n;
	@tmp_flags=$filter->apply(\@tmp_flags);
	return @tmp_flags;
}

sub _init{
	my $self=shift;
	my ($my_matrix,$stub_value)=@_;
	$self->SUPER::_init(@_);
	$self->{HOST_MATRIX}=$my_matrix;
	$self->{STUB_VALUE}=$stub_value;
	$self->{FLAGS};
}

sub _init_copy{
	my $self=shift;
	my ($filter,$my_matrix)=@_;
	$self->SUPER::_init_copy(@_);
	$self->{HOST_MATRIX}=$my_matrix;
	$self->{STUB_VALUE}=$filter->{STUB_VALUE};
	$self->{FLAGS};
	my @tmp_flags=$self->_filter2flags($filter);
	if(@tmp_flags){
		$self->{FLAGS}=[];
		@{$self->{FLAGS}}=@tmp_flags;
	}
}

sub _init_transp{
	my $self=shift;
	my ($filter,$my_matrix)=@_;
	$self->SUPER::_init_transp(@_);
	$self->{HOST_MATRIX}=$my_matrix;
	$self->{STUB_VALUE}=$filter->{STUB_VALUE};
	$self->{FLAGS};
	my @tmp_flags=$self->_filter2flags($filter);
	if(@tmp_flags){
		$self->{FLAGS}=[];
		for(my $i=0;$i<$my_matrix->{NLINES};$i++){
			($bg_idx,$fg_idx)=$my_matrix->line2sites_idxs($i);
			my $ti=$filter->{HOST_MATRIX}->site_idx_pair2line($fg_idx,$bg_idx);
			$self->{FLAGS}->[$i]=$tmp_flags[$ti];
		}
	}
}

#interface

#performs job
sub apply{
	my $self=shift;
	$self->SUPER::apply(@_);
	my ($ra_inarray)=@_;
	my @tmp=@{$ra_inarray};
	for(my $i=0;$i<@tmp;$i++){
		$tmp[$i]=$self->{STUB_VALUE} unless $self->{FLAGS}->[$i];
	}
	return @tmp;
}

sub get_my_matrix{ 
	my $self=shift;
	return $self->{HOST_MATRIX};
}

sub stub_value{
	return $self->{STUB_VALUE};
}

1;