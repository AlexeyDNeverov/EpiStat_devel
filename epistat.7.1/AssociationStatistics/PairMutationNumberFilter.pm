package AssociationStatistics::PairMutationNumberFilter;
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
	my ($my_matrix,$epistat_fn,$bgr_min_muts,$fgr_min_muts,$stub_value)=@_;
	$self->SUPER::_init(@_);
	$self->{HOST_MATRIX}=$my_matrix;
	$self->{STUB_VALUE}=0;
	$self->{STUB_VALUE}=$stub_value if(defined($stub_value));
	$self->{BGR_MIN_MUT}=0;
	$self->{FGR_MIN_MUT}=0;
	$self->{BGR_MUT_NUMBERS};
	$self->{FGR_MUT_NUMBERS};
	if(defined $bgr_min_muts){
		$self->{BGR_MIN_MUT}=$bgr_min_muts;
		$self->{BGR_MUT_NUMBERS}=[];
	}
	if(defined $fgr_min_muts){
		$self->{FGR_MIN_MUT}=$fgr_min_muts;
		$self->{FGR_MUT_NUMBERS}=[];
	}
	open INPF, "<$epistat_fn" or die "\nUnable to open input file $epistat_fn!";
	while(<INPF>){
		chomp;
		$_=~s/\s*$//;
		if($_ ne ""){
			my @lines=split "\t";
			die "\nError in AssociationStatistics::MutationNumberFilter::_init(): expacts four columns in the file $epistat_fn!" unless @lines==4;
			push @{$self->{BGR_MUT_NUMBERS}},$lines[2] if defined $self->{BGR_MUT_NUMBERS};
			push @{$self->{FGR_MUT_NUMBERS}},$lines[3] if defined $self->{FGR_MUT_NUMBERS};
		}
	}
	close INPF;
	if($self->{BGR_MIN_MUT}>0){
		die "\nError in AssociationStatistics::MutationNumberFilter::_init(): Undefined numbers of mutations in background sites!" unless @{$self->{BGR_MUT_NUMBERS}}==$my_matrix->{NLINES};
	}
	if($self->{FGR_MIN_MUT}>0){
		die "\nError in AssociationStatistics::MutationNumberFilter::_init(): Undefined numbers of mutations in foreground sites!" unless @{$self->{FGR_MUT_NUMBERS}}==$my_matrix->{NLINES};
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
	$self->{BGR_MUT_NUMBERS};
	$self->{FGR_MUT_NUMBERS};
	if(defined $filter->{BGR_MUT_NUMBERS}){
		$self->{BGR_MUT_NUMBERS}=[];
		@{$self->{BGR_MUT_NUMBERS}}=@{$filter->{BGR_MUT_NUMBERS}};
	}
	if(defined $filter->{FGR_MUT_NUMBERS}){
		$self->{FGR_MUT_NUMBERS}=[];
		@{$self->{FGR_MUT_NUMBERS}}=@{$filter->{FGR_MUT_NUMBERS}};
	}
	if($self->{BGR_MIN_MUT}>0){
		die "\nError in AssociationStatistics::MutationNumberFilter::_init_copy(): Undefined numbers of mutations in background sites!" unless @{$self->{BGR_MUT_NUMBERS}}==$my_matrix->{NLINES};
	}
	if($self->{FGR_MIN_MUT}>0){
		die "\nError in AssociationStatistics::MutationNumberFilter::_init_copy(): Undefined numbers of mutations in foreground sites!" unless @{$self->{FGR_MUT_NUMBERS}}==$my_matrix->{NLINES};
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
	$self->{BGR_MUT_NUMBERS};
	$self->{FGR_MUT_NUMBERS};
	if(defined $filter->{BGR_MUT_NUMBERS}){
		$self->{FGR_MUT_NUMBERS}=[];
	}
	if(defined $filter->{FGR_MUT_NUMBERS}){
		$self->{BGR_MUT_NUMBERS}=[];
	}
	for(my $i=0;$i<$my_matrix->{NLINES};$i++){
		my ($bgs,$fgs)=$my_matrix->line2site_pair($i);
		my $ti=$filter->{HOST_MATRIX}->site_pair2line($fgs,$bgs);
		$self->{FGR_MUT_NUMBERS}->[$i]=$filter->{BGR_MUT_NUMBERS}->[$ti] if defined $self->{FGR_MUT_NUMBERS};
		$self->{BGR_MUT_NUMBERS}->[$i]=$filter->{FGR_MUT_NUMBERS}->[$ti] if defined $self->{BGR_MUT_NUMBERS};
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
		if($self->{BGR_MIN_MUT}>0){
			my $n=$self->{BGR_MUT_NUMBERS}->[$i];
			die "\nError in AssociationStatistics::MutationNumberFilter::apply(): Undefined mutations' numbers for line $i!" unless defined($n);
			$val=$self->{STUB_VALUE} if $n<$self->{BGR_MIN_MUT};
		}
		if($self->{FGR_MIN_MUT}>0){
			my $n=$self->{FGR_MUT_NUMBERS}->[$i];
			die "\nError in AssociationStatistics::MutationNumberFilter::apply(): Undefined mutations' numbers for line $i!" unless defined($n);
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