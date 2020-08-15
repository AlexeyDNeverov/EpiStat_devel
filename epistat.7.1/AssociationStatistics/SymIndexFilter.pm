package AssociationStatistics::SymIndexFilter;
#Do not use this class directly!
use AssociationStatistics::BasicIndexFilter;
@ISA = ("AssociationStatistics::BasicIndexFilter");

sub _init{
	my $self=shift;
	my ($my_matrix,$filter1,$filter2)=@_;
	die "\nError: The host matrix is undefined!" unless defined $my_matrix;
	die "\nError AssociationStatistics::SymIndexFilter::_init(): at least one filter parameter is required!" unless defined($filter1)||defined($filter2);
	my $stub_value;
	my $f_intragene;
	my @tmp_flags1;
	my @tmp_flags2;
	if(defined $filter1){
		$stub_value=$filter1->{STUB_VALUE};
		$f_intragene=$filter1->{HOST_MATRIX}->{F_INTRAGENE};
		@tmp_flags1=$self->_filter2flags($filter1);
	}
	if(defined $filter2){
		if(!defined $stub_value){
			$stub_value=$filter2->{STUB_VALUE};
		}elsif($stub_value!=$filter2->{STUB_VALUE}){
			die "\nError: Different stub values are not allowed!";
		}
		die "\nError AssociationStatistics::SymIndexFilter::_init(): 2-nd filter argument is allowed only for intergenic case!" if $filter2->{HOST_MATRIX}->{F_INTRAGENE};
		$f_intragene=$filter2->{HOST_MATRIX}->{F_INTRAGENE} unless defined $f_intragene;
		@tmp_flags2=$self->_filter2flags($filter2);
	}
	$self->SUPER::_init($my_matrix,$stub_value);
	$self->{FLAGS}=[];
	for(my $i=0;$i<$my_matrix->{NLINES};$i++){
		($bgs,$fgs)=$my_matrix->line2site_pair($i);
		my ($i1,$i2);
		$i1=$filter1->{HOST_MATRIX}->site_pair2line($bgs,$fgs) if defined $filter1;
		my $flag=0;
		$flag=$tmp_flags1[$i1] if defined $i1;
		if(!$f_intragene){
			$i2=$filter2->{HOST_MATRIX}->site_pair2line($fgs,$bgs) if defined $filter2;
			$flag=($flag||$tmp_flags2[$i2]) if defined $i2;
		}else{
			$i2=$filter1->{HOST_MATRIX}->site_pair2line($fgs,$bgs) if defined $filter1;
			$flag=($flag||$tmp_flags1[$i2]) if defined $i2;
		}
		$self->{FLAGS}->[$i]=$flag;
	}
}

1;