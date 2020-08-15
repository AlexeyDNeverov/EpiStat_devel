package AssociationStatistics::PvalueFilter;
#Do not use this class directly!
use AssociationStatistics::basicAssociationFilter;
@ISA = ("AssociationStatistics::basicAssociationFilter");

sub _read_pvalues{
	my ($pvalues_fn)=@_;
	my @pvals;
	open INPF, "<$pvalues_fn" or die "\nUnable to open input file $pvalues_fn!";
	while(<INPF>){
		chomp;
		$_=~s/\s*$//;
		if($_ ne ""){
			push @pvals,$_;
		};
	}
	close INPF;
	return @pvals;
}

sub new{
	my $class=shift;
	my $self={};
	bless($self,$class);
	$self->_init(@_);
	return $self;
}

sub _init{
	my $self=shift;
	my ($my_matrix,$lower_pvalue_fn,$lower_pvalue_threshold,$upper_pvalue_fn,$upper_pvalue_threshold,$stub_value)=@_;
	$self->SUPER::_init(@_);
	$self->{HOST_MATRIX}=$my_matrix;
	$self->{LOWER_PVALUE_THRESHOLD};
	$self->{LOWER_PVALUES};
	$self->{UPPER_PVALUE_THRESHOLD};
	$self->{UPPER_PVALUES};
	$self->{STUB_VALUE}=0;
	die "\nError in AssociationStatistics::PvalueFilter::init_(): Undefined both pvalue thresholds!" unless defined($lower_pvalue_threshold)||defined($upper_pvalue_threshold);
	if(defined $lower_pvalue_threshold){
		die "\nError in AssociationStatistics::PvalueFilter::init_(): Unpropper pvalue threshold!" unless 0<=$lower_pvalue_threshold&&$lower_pvalue_threshold<=1;
		$self->{LOWER_PVALUE_THRESHOLD}=$lower_pvalue_threshold;
		$self->{LOWER_PVALUES}=[];
		@{$self->{LOWER_PVALUES}}=_read_pvalues($lower_pvalue_fn);
	}
	if(defined $upper_pvalue_threshold){
		die "\nError in AssociationStatistics::PvalueFilter::init_(): Unpropper pvalue threshold!" unless 0<=$upper_pvalue_threshold&&$upper_pvalue_threshold<=1;
		$self->{UPPER_PVALUE_THRESHOLD}=$upper_pvalue_threshold;
		$self->{UPPER_PVALUES}=[];
		@{$self->{UPPER_PVALUES}}=_read_pvalues($upper_pvalue_fn);
	}
	$self->{STUB_VALUE}=$stub_value if(defined($stub_value));
}

sub _init_copy{
	my $self=shift;
	my ($filter,$my_matrix)=@_;
	$self->SUPER::_init_copy(@_);
	$self->{HOST_MATRIX}=$my_matrix;
	$self->{LOWER_PVALUE_THRESHOLD}=$filter->{LOWER_PVALUE_THRESHOLD};
	$self->{LOWER_PVALUES};
	if(defined $filter->{LOWER_PVALUES}){
		$self->{LOWER_PVALUES}=[];
		@{$self->{LOWER_PVALUES}}=@{$filter->{LOWER_PVALUES}};
	}
	$self->{UPPER_PVALUE_THRESHOLD}=$filter->{UPPER_PVALUE_THRESHOLD};
	$self->{UPPER_PVALUES};
	if(defined $filter->{UPPER_PVALUES}){
		$self->{UPPER_PVALUES}=[];
		@{$self->{UPPER_PVALUES}}=@{$filter->{UPPER_PVALUES}};
	}
	$self->{STUB_VALUE}=$filter->{STUB_VALUE};
}

sub _init_transp{
	my $self=shift;
	my ($filter,$my_matrix)=@_;
	$self->{HOST_MATRIX}=$my_matrix;
	$self->{LOWER_PVALUE_THRESHOLD}=$filter->{LOWER_PVALUE_THRESHOLD};
	$self->{LOWER_PVALUES};
	if(defined $filter->{LOWER_PVALUES}){
		$self->{LOWER_PVALUES}=[];
	}
	$self->{UPPER_PVALUE_THRESHOLD}=$filter->{UPPER_PVALUE_THRESHOLD};
	$self->{UPPER_PVALUES};
	if(defined $filter->{UPPER_PVALUES}){
		$self->{UPPER_PVALUES}=[];
	}
	$self->{STUB_VALUE}=$filter->{STUB_VALUE};
	for(my $i=0;$i<$my_matrix->{NLINES};$i++){
		my ($bgs,$fgs)=$my_matrix->line2site_pair($i);
		my $ti=$filter->{HOST_MATRIX}->site_pair2line($fgs,$bgs);
		$self->{LOWER_PVALUES}->[$i]=$filter->{LOWER_PVALUES}->[$ti] if defined $self->{LOWER_PVALUES};
		$self->{UPPER_PVALUES}->[$i]=$filter->{UPPER_PVALUES}->[$ti] if defined $self->{UPPER_PVALUES};
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
	my $stub_value=$self->{STUB_VALUE};
	if(defined $self->{LOWER_PVALUE_THRESHOLD}){
		die "\nError pvalue_filter(): The size of element array is not equal to the size of pvalue array!" unless @{$ra_elements}==@{$self->{LOWER_PVALUES}};
	}
	if(defined $self->{UPPER_PVALUE_THRESHOLD}){
		die "\nError pvalue_filter(): The size of element array is not equal to the size of pvalue array!" unless @{$ra_elements}==@{$self->{UPPER_PVALUES}};
	}
	my @elements;
	if(defined($self->{LOWER_PVALUE_THRESHOLD}) || defined($self->{UPPER_PVALUE_THRESHOLD})){
		for(my $i=0;$i<@{$ra_elements};$i++){
			my $t;
			if(defined $self->{LOWER_PVALUE_THRESHOLD}){
				$t=1 if $self->{LOWER_PVALUES}->[$i]<=$self->{LOWER_PVALUE_THRESHOLD};
			}
			if(defined $self->{UPPER_PVALUE_THRESHOLD}){
				$t=1 if $self->{UPPER_PVALUES}->[$i]<=$self->{UPPER_PVALUE_THRESHOLD};
			}
			if($t){
				push @elements, $ra_elements->[$i];
			}else{
				push @elements, $stub_value;
			}
		}
	}else{
		@elements=@{$ra_elements};
	}
	return @elements;
}
#
sub stub_value{
	my $self=shift;
	return $self->{STUB_VALUE};
}
1;