package AssociationStatistics::SiteMutationNumberFilter;
#Do not use this class directly!
use AssociationStatistics::basicAssociationFilter;
@ISA = ("AssociationStatistics::basicAssociationFilter");

sub new{
	my $class=shift;
	my $self={};
	bless($self,$class);
	$self->_init(@_);
	return $self;
}

sub _init{
	my $self=shift;
	my ($my_matrix,$bgr_min_muts,$fgr_min_muts,$stub_value)=@_;
	$self->SUPER::_init(@_);
	$self->{HOST_MATRIX}=$my_matrix;
	$self->{STUB_VALUE}=0;
	$self->{STUB_VALUE}=$stub_value if(defined($stub_value));
	$self->{BGR_MIN_MUT}=0;
	$self->{FGR_MIN_MUT}=0;
	if(defined $bgr_min_muts){
		$self->{BGR_MIN_MUT}=$bgr_min_muts;
	}
	if(defined $fgr_min_muts){
		$self->{FGR_MIN_MUT}=$fgr_min_muts;
	}
	if($self->{BGR_MIN_MUT}>0){
		die "\nError in AssociationStatistics::SiteMutationNumberFilter::_init(): Undefined numbers of mutations in background sites!" unless defined $my_matrix->{BGR_MUT_NUMBERS};
	}
	if($self->{FGR_MIN_MUT}>0){
		die "\nError in AssociationStatistics::SiteMutationNumberFilter::_init(): Undefined numbers of mutations in foreground sites!" unless defined $my_matrix->{FGR_MUT_NUMBERS};
	}
}

sub _init_copy{
	my $self=shift;
	my ($filter,$my_matrix)=@_;
	$self->SUPER::_init_copy(@_);
	$self->{HOST_MATRIX}=$my_matrix;
	$self->{STUB_VALUE}=$filter->{STUB_VALUE};
	$self->{BGR_MIN_MUT}=$filter->{BGR_MIN_MUT};
	$self->{FGR_MIN_MUT}=$filter->{FGR_MIN_MUT};
	if($self->{BGR_MIN_MUT}>0){
		die "\nError in AssociationStatistics::SiteMutationNumberFilter::_init_copy(): Undefined numbers of mutations in background sites!" unless defined $my_matrix->{BGR_MUT_NUMBERS};
	}
	if($self->{FGR_MIN_MUT}>0){
		die "\nError in AssociationStatistics::SiteMutationNumberFilter::_init_copy(): Undefined numbers of mutations in foreground sites!" unless defined $my_matrix->{FGR_MUT_NUMBERS};
	}
}

sub _init_transp{
	my $self=shift;
	my ($filter,$my_matrix)=@_;
	$self->SUPER::_init_transp(@_);
	$self->{HOST_MATRIX}=$my_matrix;
	$self->{STUB_VALUE}=$filter->{STUB_VALUE};
	$self->{BGR_MIN_MUT}=$filter->{FGR_MIN_MUT};
	$self->{FGR_MIN_MUT}=$filter->{BGR_MIN_MUT};
	if($self->{BGR_MIN_MUT}>0){
		die "\nError in AssociationStatistics::SiteMutationNumberFilter::_init_transp(): Undefined numbers of mutations in background sites!" unless defined $my_matrix->{BGR_MUT_NUMBERS};
	}
	if($self->{FGR_MIN_MUT}>0){
		die "\nError in AssociationStatistics::SiteMutationNumberFilter::_init_transp(): Undefined numbers of mutations in foreground sites!" unless defined $my_matrix->{FGR_MUT_NUMBERS};
	}
}

#returns a reference to corresponding matrix object
sub get_my_matrix{ 
	my $self=shift;
	return $self->{HOST_MATRIX};
}
#performs job
sub apply{
	my $self=shift;
	my ($ra_elements)=@_;
	$self->SUPER::apply(@_);
	my $my_matrix=$self->{HOST_MATRIX};
	my @elements;
	for(my $i=0;$i<@{$ra_elements};$i++){
		my $val=$ra_elements->[$i];
		my ($bgr_idx,$fgr_idx)=$my_matrix->line2sites_idxs($i);
		if($self->{BGR_MIN_MUT}>0){
			my $n=$my_matrix->{BGR_MUT_NUMBERS}->[$bgr_idx];
			die "\nError in AssociationStatistics::SiteMutationNumberFilter::apply(): Undefined mutations' numbers for line $i!" unless defined($n);
			$val=$self->{STUB_VALUE} if $n<$self->{BGR_MIN_MUT};
		}
		if($self->{FGR_MIN_MUT}>0){
			my $n=$my_matrix->{FGR_MUT_NUMBERS}->[$fgr_idx];
			die "\nError in AssociationStatistics::SiteMutationNumberFilter::apply(): Undefined mutations' numbers for line $i!" unless defined($n);
			$val=$self->{STUB_VALUE} if $n<$self->{FGR_MIN_MUT};
		}
		push @elements,$val;
	}
	return @elements;
}
#
sub stub_value{
	return $self->{STUB_VALUE};
}
1;