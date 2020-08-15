package AssociationStatistics::SymBaseAssociationMeasure;
#Module for calculating covariation/correlation-like association statistics for a site pair on the base of epistatic statistics
use AssociationStatistics::BaseAssociationMeasure;
use Scalar::Util;
use AssociationStatistics::SymNonNegativeValueFilter;
use AssociationStatistics::SymIndexFilter;
@ISA = ("AssociationStatistics::BaseAssociationMeasure");

sub _init{
	my $self=shift;
	#args: BaseAssociationMeasure,BaseAssociationMeasure,hashref
	my ($spm_dataset1,$spm_dataset2)=@_;
	die "\nError SymBaseAssociationMeasure::_init(): the first parameter is required!" unless defined $spm_dataset1;
	my %bgr_sites;
	my %fgr_sites;
	$spm_dataset1->get_sites(\%bgr_sites,\%fgr_sites,'-nofilters' => 1);
	if(defined $spm_dataset2){
		die "\nError SymBaseAssociationMeasure::_init(): if the second parameter is defined than both passed matrices have to be intergenic!" if $spm_dataset1->{F_INTRAGENE}||$spm_dataset2->{F_INTRAGENE};
		my %bgr_sites2;
		my %fgr_sites2;
		$spm_dataset2->get_sites(\%bgr_sites2,\%fgr_sites2,'-nofilters' => 1);
		foreach my $site(keys %bgr_sites2){
			$fgr_sites{$site}++;
		}
		foreach my $site(keys %fgr_sites2){
			$bgr_sites{$site}++;
		}
	}else{
		foreach my $site(keys %fgr_sites){
			$bgr_sites{$site}++;
		}
		%fgr_sites=();
	}
	foreach my $site(keys %bgr_sites){
		delete $bgr_sites->{$site} unless $bgr_sites->{$site} == 2;
	}
	foreach my $site(keys %fgr_sites){
		delete $fgr_sites->{$site} unless $fgr_sites->{$site} == 2;
	}
	if(defined $spm_dataset2){
		$self->SUPER::_init_copy($spm_dataset1,\%bgr_sites,\%fgr_sites,'-nofilters' => 1);
	}else{
		%fgr_sites=%bgr_sites;
		$self->SUPER::_init_copy($spm_dataset1,\%bgr_sites,undef,'-nofilters' => 1);
	}
	#Summation statistics for elements and their transpositions
	my @tmp_obs;
	my $j=0;
	my $sp_matrix=$spm_dataset2;
	$sp_matrix=$spm_dataset1 unless defined $sp_matrix;
	for(my $i=0;$i<$sp_matrix->{NLINES};$i++){
		($bgs,$fgs)=$sp_matrix->line2site_pair($i);
		if(exists($bgr_sites{$fgs})&&exists($fgr_sites{$bgs})){
			my $k=$self->site_pair2line($fgs,$bgs);
			if(defined $k){
				$self->{MEAN}->[$k]+=$sp_matrix->{MEAN}->[$i];
				if($sp_matrix->{OBS_LINES}->[$j]==$i){
					push @tmp_obs,[($k,$sp_matrix->{OBS}->[$j])];
				}
			}
		}
		$j++ if($sp_matrix->{OBS_LINES}->[$j]==$i);
	}
	for(my $i=0;$i<@{$self->{OBS_LINES}};$i++){
		push @tmp_obs,[($self->{OBS_LINES}->[$i],$self->{OBS}->[$i])];
	};
	@tmp_obs=sort {$a->[0]<=>$b->[0]} @tmp_obs;
	@{$self->{OBS_LINES}}=();
	@{$self->{OBS}}=();
	$j=0;
	while($j<@tmp_obs){
		my $obs=$tmp_obs[$j]->[1];
		my $k=$tmp_obs[$j]->[0];
		if($j<$#tmp_obs&&$tmp_obs[$j]->[0]==$tmp_obs[$j+1]->[0]){
			$obs+=$tmp_obs[$j+1]->[1];
			$j++;
		}
		push @{$self->{OBS}},$obs;
		push @{$self->{OBS_LINES}},$k;
		$j++;
	}
	#Init filters
	@{$self->{FILTERS}}=();
	my %tmp_filters;
	foreach my $f(@{$spm_dataset1->{FILTERS}}){
		my $type_name=Scalar::Util::blessed($f);
		if(defined $type_name){
			$tmp_filters{$type_name}=[] unless exists $tmp_filters{$type_name};
			$tmp_filters{$type_name}->[0]=$f;
		}
	}
	if(defined $spm_dataset2){
		foreach my $f(@{$spm_dataset2->{FILTERS}}){
			my $type_name=Scalar::Util::blessed($f);
			if(defined $type_name){
				$tmp_filters{$type_name}=[] unless exists $tmp_filters{$type_name};
				$tmp_filters{$type_name}->[1]=$f;
			}
		}
	}
	foreach my $type_name(keys %tmp_filters){
		my $f;
		if($type_name =~/::NonNegativeValueFilter/){
			$f=AssociationStatistics::SymNonNegativeValueFilter->new($self,@{$tmp_filters{$type_name}});
		}elsif($type_name =~/::PvalueFilter/ ||
			$type_name =~/::PairMutationNumberFilter/ ||
			$type_name =~/::SiteMutationNumberFilter/ ||
			$type_name =~/::SelectedSitesFilter/){
			$f=AssociationStatistics::SymIndexFilter->new($self,@{$tmp_filters{$type_name}});
		}else{
			die "\nError SymBaseAssociationMeasure::_init(): unknown filter class: $type_name!";
		}
		push @{$self->{FILTERS}},$f;
	}
}

1;