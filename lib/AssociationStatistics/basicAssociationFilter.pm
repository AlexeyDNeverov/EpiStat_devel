package AssociationStatistics::basicAssociationFilter;
#Do not use this class directly!
#This module provides interface for a filter object of association matrix elements
#Filter checks some logical condition on each element of association matrix and if it is wrong, drop association off - sets it's value to a stab (usually 0)
#Filter object is instantiated only by corresponding matrix object. 
sub _init{
	my $self=shift;
	my ($my_matrix,$stub_value)=@_;
	die "\nError basicAssociationFilter::_init(): The host matrix is undefined!" unless defined $my_matrix;
}

sub _init_copy{
	my $self=shift;
	my ($filter,$my_matrix)=@_;
	die "\nError basicAssociationFilter::_init_copy(): The host matrix is undefined!" unless defined $my_matrix;
}

sub _init_transp{
	my $self=shift;
	my ($filter,$my_matrix)=@_;
	$self->_init_copy(@_);
}

#interface
sub new{
	my $class=shift;
	my $self={};
	bless($self,$class);
	$self->_init(@_);
	return $self;
}

sub copy{
	my $proto=shift;
	my $class = ref($proto) || $proto;
	my $parent = ref($proto) && $proto;
	my $child={};
	bless($child,$class);
	$child->_init_copy($parent,@_);
	return $child;
}

sub transpose{
	my $proto=shift;
	my $class = ref($proto) || $proto;
	my $parent = ref($proto) && $proto;
	my $child={};
	bless($child,$class);
	$child->_init_transp($parent,@_);
	return $child;
}

#returns a reference to corresponding matrix object
sub get_my_matrix{ 
	my $self=shift;
	return undef;
}
#performs job
sub apply{
	my $self=shift;
	my ($ra_inarray)=@_;
	my $my_mtx=$self->get_my_matrix;
	if(defined $ra_inarray){
		die "\nError in AssociationStatistics::basicAssociationFilter::apply(): Input array have to be the same size as the array of corresponding matrix!" unless $my_mtx->{NLINES}==@{$ra_inarray};
	}
	return ();
}
#
sub stub_value{
	return undef;
}

1;