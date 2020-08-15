package AssociationStatistics::SymNonNegativeValueFilter;
#Do not use this class directly!
use AssociationStatistics::NonNegativeValueFilter;
@ISA = ("AssociationStatistics::NonNegativeValueFilter");

sub _init{
	my $self=shift;
	#args: BaseAssociationMeasure,NonNegativeValueFilter
	my ($my_matrix,$filter1,$filter2)=@_;
	die "\nError: The host matrix is undefined!" unless defined $my_matrix;
	my $stub_value;
	if(defined $filter1){
		$stub_value=$filter1->{STUB_VALUE};
	}
	if(defined $filter2){
		if(!defined $stub_value){
			$stub_value=$filter2->{STUB_VALUE};
		}elsif($stub_value!=$filter2->{STUB_VALUE}){
			die "\nError: Different stub values are not allowed!";
		}
	}
	die "\nError: At least one filter must be defined!" unless defined($filter1)||defined($filter2);
	$self->SUPER::_init($my_matrix,$stub_value);
}

1;