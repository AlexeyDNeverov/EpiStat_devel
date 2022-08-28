package AssociationStatistics::SymCorrelationLike;
#Module for calculating covariance/correlation-like association statistics for a site pair on the base of epistatic statistics
use AssociationStatistics::SymBaseAssociationMeasure;
use AssociationStatistics::CorrelationLike;
@ISA = ("AssociationStatistics::SymBaseAssociationMeasure");

sub _init{
	my $self=shift;
	my ($dataset1_desc,$dataset2_desc,$general_settings)=@_;
	die "\nError AssociationStatistics::SymCorrelationLike::_init(): the description of the first dataset is required!" unless defined $dataset1_desc;
	my $f_intragene=$dataset1_desc->f_intragene_pairs;
	if(defined $dataset2_desc){
		if($dataset1_desc==$dataset2_desc){
			$dataset2_desc=undef;
		}elsif(!$f_intragene){
			die "\nError AssociationStatistics::SymCorrelationLike::_init(): unmatched datasets accounted!" unless (!$dataset2_desc->f_intragene_pairs)&&
				($dataset1_desc->bgr_gene_name eq $dataset2_desc->fgr_gene_name)&&
				($dataset1_desc->fgr_gene_name eq $dataset2_desc->bgr_gene_name);
		}else{
			die "\nError AssociationStatistics::SymCorrelationLike::_init(): unknown 2-nd dataset accounted for intragene matrix!";
		}
	}
	my ($matrix1,$matrix2);
	$matrix1=AssociationStatistics::CorrelationLike->new($dataset1_desc,$general_settings);
	$matrix2=AssociationStatistics::CorrelationLike->new($dataset2_desc,$general_settings) if defined $dataset2_desc;
	$self->SUPER::_init($matrix1,$matrix2);
	$self->{BGR_SITE_VAR}={};
	$self->{FGR_SITE_VAR}={};
	foreach my $site(keys %{$matrix1->{BGR_SITE_VAR}}){
		$self->{BGR_SITE_VAR}->{$site}=$matrix1->{BGR_SITE_VAR}->{$site};
		if(defined $matrix2){
			$self->{BGR_SITE_VAR}->{$site}+=$matrix2->{FGR_SITE_VAR}->{$site};
		}else{
			$self->{BGR_SITE_VAR}->{$site}+=$matrix1->{FGR_SITE_VAR}->{$site};
		}
	}
	foreach my $site(keys %{$matrix1->{FGR_SITE_VAR}}){
		$self->{FGR_SITE_VAR}->{$site}=$matrix1->{FGR_SITE_VAR}->{$site};
		if(defined $matrix2){
			$self->{FGR_SITE_VAR}->{$site}+=$matrix2->{BGR_SITE_VAR}->{$site};
		}else{
			$self->{FGR_SITE_VAR}->{$site}+=$matrix1->{BGR_SITE_VAR}->{$site};
		}
	}
}

sub _init_copy{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in SymCorrelationLike::_init_copy(): undefined matrix argument!" unless defined $sp_matrix;
	$self->SUPER::_init_copy($sp_matrix,@_);
	$self->{BGR_SITE_VAR}={};
	$self->{FGR_SITE_VAR}={};
	%{$self->{BGR_SITE_VAR}}=%{$sp_matrix->{BGR_SITE_VAR}};
	%{$self->{FGR_SITE_VAR}}=%{$sp_matrix->{FGR_SITE_VAR}};
}

sub _init_transp{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in SymCorrelationLike::_init_transp(): undefined matrix argument!" unless defined $sp_matrix;
	$self->SUPER::_init_transp($sp_matrix,@_);
	$self->{BGR_SITE_VAR}={};
	$self->{FGR_SITE_VAR}={};
	%{$self->{BGR_SITE_VAR}}=%{$sp_matrix->{FGR_SITE_VAR}};
	%{$self->{FGR_SITE_VAR}}=%{$sp_matrix->{BGR_SITE_VAR}};
}

#interface declaration

sub mtx_diag_value{
	my $self=shift;
	return 1.0;
}

#$norm_const==undef - default normalization
#$norm_const==0 - no normalization
#$norm_const!=0 - the value is used for normalization
sub get_statistics{
	my $self=shift;
	my ($norm_const)=@_;
	my @stat=$self->SUPER::get_statistics;
	my $M=0;
	for(my $i=0;$i<@stat;$i++){
		$N+=$stat[$i];
		my ($j,$k)=$self->line2site_pair($i);
		my $m=sqrt($self->{BGR_SITE_VAR}->{$j})*sqrt($self->{FGR_SITE_VAR}->{$k});
		if($m!=0){
			$stat[$i]/=$m;
			$M=abs($stat[$i]) if abs($stat[$i])>$M;
		}else{
			$stat[$i]="NA";
		}
	}
	$M=$norm_const if defined $norm_const;
	if($M){
		for(my $i=0;$i<@stat;$i++){
			$stat[$i]/=$M if($stat[$i] ne "NA");
		}
	}
	return @stat;
}

1;
	