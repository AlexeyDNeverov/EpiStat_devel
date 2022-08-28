package AssociationStatistics::SymNormZScore;
#Module for calculating covariation/correlation-like association statistics for a site pair on the base of epistatic statistics
use AssociationStatistics::SymBaseAssociationMeasure;
use AssociationStatistics::NormZScore;
@ISA = ("AssociationStatistics::SymBaseAssociationMeasure");

sub _init{
	my $self=shift;
	my ($dataset1_desc,$dataset2_desc,$general_settings)=@_;
	die "\nError AssociationStatistics::SymNormZScore::_init(): the description of the first dataset is required!" unless defined $dataset1_desc;
	my $f_intragene=$dataset1_desc->f_intragene_pairs;
	if(defined $dataset2_desc){
		if($dataset1_desc==$dataset2_desc){
			$dataset2_desc=undef;
		}elsif(!$f_intragene){
			die "\nError AssociationStatistics::SymNormZScore::_init(): unmatched datasets accounted!" unless (!$dataset2_desc->f_intragene_pairs)&&
				($dataset1_desc->bgr_gene_name eq $dataset2_desc->fgr_gene_name)&&
				($dataset1_desc->fgr_gene_name eq $dataset2_desc->bgr_gene_name);
		}else{
			die "\nError AssociationStatistics::SymNormZScore::_init(): unknown 2-nd dataset accounted for intragene matrix!";
		}
	}
	my ($matrix1,$matrix2);
	$matrix1=AssociationStatistics::NormZScore->new($dataset1_desc,$general_settings);
	$matrix2=AssociationStatistics::NormZScore->new($dataset2_desc,$general_settings) if defined $dataset2_desc;
	$self->SUPER::_init($matrix1,$matrix2);
	$self->{VARIANCE}=[];
	for(my $i=0;$i<$self->{NLINES};$i++){
		my ($bgs,$fgs)=$self->line2site_pair($i);
		my $i1=$matrix1->site_pair2line($bgs,$fgs);
		$self->{VARIANCE}->[$i]+=$matrix1->{VARIANCE}->[$i1] if defined $i1;
		my $i2;
		if(defined $matrix2){
			$i2=$matrix2->site_pair2line($fgs,$bgs);
			$self->{VARIANCE}->[$i]+=$matrix2->{VARIANCE}->[$i2] if defined $i2;
		}else{
			$i2=$matrix1->site_pair2line($fgs,$bgs);
			$self->{VARIANCE}->[$i]+=$matrix1->{VARIANCE}->[$i2] if defined $i2;
		}
	}
	$self->{FGR_SITE2MUT_NUMBERS}={};
	$self->{FGR_SITE_MEAN}={};
	$self->{FGR_SITE_VAR}={};
	$self->{BGR_SITE2MUT_NUMBERS}={};
	$self->{BGR_SITE_MEAN}={};
	$self->{BGR_SITE_VAR}={};
	for(my $i=0;$i<$matrix1->get_fgr_site_num;$i++){
		my $site=$matrix1->{FGR_SITES}->[$i];
		$self->{FGR_SITE2MUT_NUMBERS}->{$site}=$matrix1->{FGR_MUT_NUMBERS}->[$i];
		$self->{FGR_SITE_MEAN}->{$site}=$matrix1->{FGR_SITE_MEAN}->{$site};
		$self->{FGR_SITE_VAR}->{$site}=$matrix1->{FGR_SITE_VAR}->{$site};
	}
	if(defined $matrix2){
		for(my $i=0;$i<$matrix2->get_fgr_site_num;$i++){
			my $site=$matrix2->{FGR_SITES}->[$i];
			$self->{BGR_SITE2MUT_NUMBERS}->{$site}=$matrix2->{FGR_MUT_NUMBERS}->[$i];
			$self->{BGR_SITE_MEAN}->{$site}=$matrix2->{FGR_SITE_MEAN}->{$site};
			$self->{BGR_SITE_VAR}->{$site}=$matrix2->{FGR_SITE_VAR}->{$site};
		}
	}else{
		%{$self->{BGR_SITE2MUT_NUMBERS}}=%{$self->{FGR_SITE2MUT_NUMBERS}};
		%{$self->{BGR_SITE_MEAN}}=%{$self->{FGR_SITE_MEAN}};
		%{$self->{BGR_SITE_VAR}}=%{$self->{FGR_SITE_VAR}};
	}
}

sub _init_copy{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in SymNormZScore::_init_copy(): undefined matrix argument!" unless defined $sp_matrix;
	$self->SUPER::_init_copy($sp_matrix,@_);
	$self->{VARIANCE}=[];
	@{$self->{VARIANCE}}=@{$sp_matrix->{VARIANCE}};
	%{$self->{FGR_SITE2MUT_NUMBERS}}=%{$sp_matrix->{FGR_SITE2MUT_NUMBERS}};
	%{$self->{FGR_SITE_MEAN}}=%{$sp_matrix->{FGR_SITE_MEAN}};
	%{$self->{FGR_SITE_VAR}}=%{$sp_matrix->{FGR_SITE_VAR}};
	%{$self->{BGR_SITE2MUT_NUMBERS}}=%{$sp_matrix->{BGR_SITE2MUT_NUMBERS}};
	%{$self->{BGR_SITE_MEAN}}=%{$sp_matrix->{BGR_SITE_MEAN}};
	%{$self->{BGR_SITE_VAR}}=%{$sp_matrix->{BGR_SITE_VAR}};
}

sub _init_transp{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in SymNormZScore::_init_transp(): undefined matrix argument!" unless defined $sp_matrix;
	$self->SUPER::_init_transp($sp_matrix,@_);
	$self->{VARIANCE}=[];
	for(my $i=0;$i<$self->{NLINES};$i++){
		my ($bg_idx,$fg_idx)=$self->line2sites_idxs($i);
		my $ti=$sp_matrix->site_idx_pair2line($fg_idx,$bg_idx);
		$self->{VARIANCE}->[$i]=$sp_matrix->{VARIANCE}->[$ti];
	}
	%{$self->{FGR_SITE2MUT_NUMBERS}}=%{$sp_matrix->{BGR_SITE2MUT_NUMBERS}};
	%{$self->{FGR_SITE_MEAN}}=%{$sp_matrix->{BGR_SITE_MEAN}};
	%{$self->{FGR_SITE_VAR}}=%{$sp_matrix->{BGR_SITE_VAR}};
	%{$self->{BGR_SITE2MUT_NUMBERS}}=%{$sp_matrix->{FGR_SITE2MUT_NUMBERS}};
	%{$self->{BGR_SITE_MEAN}}=%{$sp_matrix->{FGR_SITE_MEAN}};
	%{$self->{BGR_SITE_VAR}}=%{$sp_matrix->{FGR_SITE_VAR}};
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
	#my $M=0;
	for(my $i=0;$i<$self->{NLINES};$i++){
		my ($bgs,$fgs)=$self->line2site_pair($i);
		my $m=($self->{FGR_SITE2MUT_NUMBERS}->{$fgs}+$self->{BGR_SITE2MUT_NUMBERS}->{$bgs}-1)-$self->{BGR_SITE_MEAN}->{$bgs}-$self->{FGR_SITE_MEAN}->{$fgs};
		#$mbgr=($self->{BGR_SITE2MUT_NUMBERS}->{$bgs}-$self->{BGR_SITE_MEAN}->{$bgs})/sqrt($self->{BGR_SITE_VAR}->{$bgs});
		#$mfgr=($self->{FGR_SITE2MUT_NUMBERS}->{$fgs}-$self->{FGR_SITE_MEAN}->{$fgs})/sqrt($self->{FGR_SITE_VAR}->{$fgs});
		#my $m=sqrt($mbgr*$mfgr);
		if($self->{VARIANCE}->[$i]){
			$stat[$i]/=$m*sqrt($self->{VARIANCE}->[$i]/($self->{FGR_SITE_VAR}->{$fgs}+$self->{BGR_SITE_VAR}->{$bgs}));
			$stat[$i]=1 if  $stat[$i]>1;
			#$stat[$i]/=$m*sqrt($self->{VARIANCE}->[$i]);
			#$M=abs($stat[$i]) if(abs($stat[$i])>$M);
		}else{
			$stat[$i]=undef;
		}
	}
	#$M=$norm_const if defined $norm_const;
	#if($M){
	#	for(my $i=0;$i<@stat;$i++){
	#		$stat[$i]/=$M if defined $stat[$i];
	#	}
	#}
	return @stat;
}

1;
	