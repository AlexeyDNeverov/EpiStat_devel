package AssociationStatistics::SelectedSitesFilter;
#Do not use this class directly!
use AssociationStatistics::BasicIndexFilter;
@ISA = ("AssociationStatistics::BasicIndexFilter");

#Undefined values of '$bg_site_list_fn' and '$fg_site_list_fn' are interpretated as no restrictions on sites has been applied
sub _init{
	my $self=shift;
	my ($my_matrix,$stub_value,$bg_site_list_fn,$fg_site_list_fn,$f_pairs_ordering,$f_and)=@_;
	$bg_site_list_fn=~s/^\s+//;
	$bg_site_list_fn=~s/\s+$//;
	$fg_site_list_fn=~s/^\s+//;
	$fg_site_list_fn=~s/\s+$//;
	$f_pairs_ordering=1 unless defined $f_pairs_ordering;
	die "\nAt lest one list of sites for the background or for the foreground is required!" unless defined($bg_site_list_fn)||defined($fg_site_list_fn);
	if(!$f_pairs_ordering){
		die "\nAssociationStatistics::SelectedSitesFilter: the criterion for selection of unordered pairs is not determined!" unless defined $f_and;
	}
	$self->SUPER::_init(@_);
	$self->{HOST_MATRIX}=$my_matrix;
	$self->{STUB_VALUE}=$stub_value;
	$self->{FLAGS}=[];
	my $n=$my_matrix->{NLINES};
	@{$self->{FLAGS}}=(0) x $n;
	my %bg_sites;
	my %fg_sites;
	$my_matrix->get_sites(\%bg_sites,\%fg_sites,-nofilters=>1);
	my %bg_selected_sites;
	my %fg_selected_sites;
	if(defined $bg_site_list_fn){
		open INFILE, "<$bg_site_list_fn" or die "\nUnable to open input file: $bg_site_list_fn!";
		while(<INFILE>){
			while(m/(\d+)/g){
				my $site=$1;
				if(defined $bg_sites{$site}){
					$bg_selected_sites{$site}=1;
				}else{
					warn "\nThe site cordinate $site is absent in the background!";
				}
			}
		}
		close INFILE;
	}elsif($f_pairs_ordering){
		%bg_selected_sites=%bg_sites;
	}
	if(defined $fg_site_list_fn){
		if(!defined($bg_site_list_fn)||($bg_site_list_fn ne $fg_site_list_fn)){
			open INFILE, "<$fg_site_list_fn" or die "\nUnable to open input file: $fg_site_list_fn!";
			while(<INFILE>){
				while(m/(\d+)/g){
					my $site=$1;
					if(defined $fg_sites{$site}){
						$fg_selected_sites{$site}=1;
					}else{
						die "\nThe site cordinate $site is absent in the foreground!";
					}
				}
			}
			close INFILE;
		}elsif($f_pairs_ordering){
			%fg_selected_sites=%bg_selected_sites if defined($bg_site_list_fn);
		}
	}elsif($f_pairs_ordering){
		%fg_selected_sites=%fg_sites;
	}
	if($f_pairs_ordering){
		for(my $i=0;$i<$n;$i++){
			my ($bgs,$fgs)=$my_matrix->line2site_pair($i);
			$self->{FLAGS}->[$i]=1 if defined($bg_selected_sites{$bgs})&&defined($fg_selected_sites{$fgs});
		}
	}else{
		if($f_and){
			for(my $i=0;$i<$n;$i++){
				my ($bgs,$fgs)=$my_matrix->line2site_pair($i);
				$self->{FLAGS}->[$i]=1 if (defined($bg_selected_sites{$bgs})||defined($fg_selected_sites{$bgs}))&&
					(defined($bg_selected_sites{$fgs})||defined($fg_selected_sites{$fgs}));
			}
		}else{
			for(my $i=0;$i<$n;$i++){
				my ($bgs,$fgs)=$my_matrix->line2site_pair($i);
				$self->{FLAGS}->[$i]=1 if defined($bg_selected_sites{$bgs})||defined($fg_selected_sites{$bgs})||
					defined($bg_selected_sites{$fgs})||defined($fg_selected_sites{$fgs});
			}
		}
	}
}

#interface
sub new{
	my $class=shift;
	my $self={};
	bless($self,$class);
	$self->_init(@_);
	return $self;
}

1;
