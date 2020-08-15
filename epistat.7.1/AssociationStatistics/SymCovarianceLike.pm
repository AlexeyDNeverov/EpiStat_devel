package AssociationStatistics::SymCovarianceLike;
#Module for calculating covariance/correlation-like association statistics for a site pair on the base of epistatic statistics
use AssociationStatistics::SymBaseAssociationMeasure;
use AssociationStatistics::CovarianceLike;
@ISA = ("AssociationStatistics::SymBaseAssociationMeasure");

sub _init{
	my $self=shift;
	my ($dataset1_desc,$dataset2_desc,$general_settings)=@_;
	die "\nError AssociationStatistics::SymCovarianceLike::_init(): the description of the first dataset is required!" unless defined $dataset1_desc;
	my $f_intragene=$dataset1_desc->f_intragene_pairs;
	if(defined $dataset2_desc){
		if($dataset1_desc==$dataset2_desc){
			$dataset2_desc=undef;
		}elsif(!$f_intragene){
			die "\nError AssociationStatistics::SymCovarianceLike::_init(): unmatched datasets accounted!" unless (!$dataset2_desc->f_intragene_pairs)&&
				($dataset1_desc->bgr_gene_name eq $dataset2_desc->fgr_gene_name)&&
				($dataset1_desc->fgr_gene_name eq $dataset2_desc->bgr_gene_name);
		}else{
			die "\nError AssociationStatistics::SymCovarianceLike::_init(): unknown 2-nd dataset accounted for intragene matrix!";
		}
	}
	my ($matrix1,$matrix2);
	$matrix1=AssociationStatistics::CovarianceLike->new($dataset1_desc,$general_settings);
	$matrix2=AssociationStatistics::CovarianceLike->new($dataset2_desc,$general_settings) if defined $dataset2_desc;
	$self->SUPER::_init($matrix1,$matrix2);
	$self->{BGR_SITE2MUT_NUMBERS}={};
	$self->{FGR_SITE2MUT_NUMBERS}={};
	for(my $i=0;$i<$matrix1->get_fgr_site_num;$i++){
		my $site=$matrix1->{FGR_SITES}->[$i];
		$self->{FGR_SITE2MUT_NUMBERS}->{$site}=$matrix1->{FGR_MUT_NUMBERS}->[$i];
	}
	if(defined $matrix2){
		for(my $i=0;$i<$matrix2->get_fgr_site_num;$i++){
			my $site=$matrix2->{FGR_SITES}->[$i];
			$self->{BGR_SITE2MUT_NUMBERS}->{$site}=$matrix2->{FGR_MUT_NUMBERS}->[$i];
		}
	}else{
		%{$self->{BGR_SITE2MUT_NUMBERS}}=%{$self->{FGR_SITE2MUT_NUMBERS}};
	}
}

sub _init_copy{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in SymCovarianceLike::_init_copy(): undefined matrix argument!" unless defined $sp_matrix;
	$self->SUPER::_init_copy($sp_matrix,@_);
	%{$self->{BGR_SITE2MUT_NUMBERS}}=%{$sp_matrix->{BGR_SITE2MUT_NUMBERS}};
	%{$self->{FGR_SITE2MUT_NUMBERS}}=%{$sp_matrix->{FGR_SITE2MUT_NUMBERS}};
}

sub _init_transp{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in SymCovarianceLike::_init_copy(): undefined matrix argument!" unless defined $sp_matrix;
	$self->SUPER::_init_transp($sp_matrix,@_);
	%{$self->{BGR_SITE2MUT_NUMBERS}}=%{$sp_matrix->{FGR_SITE2MUT_NUMBERS}};
	%{$self->{FGR_SITE2MUT_NUMBERS}}=%{$sp_matrix->{BGR_SITE2MUT_NUMBERS}};
}

#interface declaration

sub mtx_diag_value{
	my $self=shift;
	return 1.0;
}

sub get_statistics{
	my $self=shift;
	my @stat=$self->SUPER::get_statistics;
	for(my $i=0;$i<@stat;$i++){
		my ($bgs,$fgs)=$self->line2site_pair($i);
		my $m=$self->{FGR_SITE2MUT_NUMBERS}->{$fgs}+$self->{BGR_SITE2MUT_NUMBERS}->{$bgs}-1-$self->{MEAN}->[$i];
		$m=$self->{MEAN}->[$i] if $m<$self->{MEAN}->[$i];
		$stat[$i]/=$m;
	}
	return @stat;
}

1;
	