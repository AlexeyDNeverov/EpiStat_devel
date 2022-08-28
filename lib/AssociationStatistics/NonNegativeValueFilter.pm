package AssociationStatistics::NonNegativeValueFilter;
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
	my ($my_matrix,$stub_value)=@_;
	$self->SUPER::_init(@_);
	$self->{HOST_MATRIX}=$my_matrix;
	$self->{STUB_VALUE}=0;
	$self->{STUB_VALUE}=$stub_value if(defined($stub_value));
}

sub _init_copy{
	my $self=shift;
	my ($filter,$my_matrix)=@_;
	$self->SUPER::_init_copy(@_);
	$self->{HOST_MATRIX}=$my_matrix;
	$self->{STUB_VALUE}=$filter->{STUB_VALUE};
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
	my $stub_value=$self->{STUB_VALUE};
	my @elements;
	foreach my $val(@{$ra_elements}){
		if($val>=0){
			push @elements,$val;
		}else{
			push @elements,$stub_value;
		}
	}
	return @elements;
}
#
sub stub_value{
	return $self->{STUB_VALUE};
}
1;