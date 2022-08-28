package AssociationStatistics::SymZScore;
#Module for calculating covariation/correlation-like association statistics for a site pair on the base of epistatic statistics
use AssociationStatistics::SymBaseAssociationMeasure;
use AssociationStatistics::ZScore;
@ISA = ("AssociationStatistics::SymBaseAssociationMeasure");

sub _init{
	my $self=shift;
	my ($dataset1_desc,$dataset2_desc,$general_settings)=@_;
	die "\nError AssociationStatistics::SymZScore::_init(): the description of the first dataset is required!" unless defined $dataset1_desc;
	my $f_intragene=$dataset1_desc->f_intragene_pairs;
	if(defined $dataset2_desc){
		if($dataset1_desc==$dataset2_desc){
			$dataset2_desc=undef;
		}elsif(!$f_intragene){
			die "\nError AssociationStatistics::SymZScore::_init(): unmatched datasets accounted!" unless (!$dataset2_desc->f_intragene_pairs)&&
				($dataset1_desc->bgr_gene_name eq $dataset2_desc->fgr_gene_name)&&
				($dataset1_desc->fgr_gene_name eq $dataset2_desc->bgr_gene_name);
		}else{
			die "\nError AssociationStatistics::SymZScore::_init(): unknown 2-nd dataset accounted for intragene matrix!";
		}
	}
	my ($matrix1,$matrix2);
	$matrix1=AssociationStatistics::ZScore->new($dataset1_desc,$general_settings);
	$matrix2=AssociationStatistics::ZScore->new($dataset2_desc,$general_settings) if defined $dataset2_desc;
	$self->SUPER::_init($matrix1,$matrix2);
	$self->{VARIANCE}=[];
	my @osp_cov;
	if($f_intragene&& defined($dataset1_desc->ord_pairs_covariance_fn)){
		my $str=$dataset1_desc->ord_pairs_covariance_fn;
		open INPF, "<$str" or die  "\nUnable to open input file $str!";
		while(<INPF>){
			chomp;
			$_=~s/\s*$//;
			if($_ ne ""){
				push @osp_cov,$_;
			}
		}
		close INPF;
		die "\nError AssociationStatistics::SymZScore::_init(): Number of lines in covariance file $str inconsistent with the dataset ".$dataset1_desc->name."!" unless @osp_cov==$matrix1->{NLINES};
	}
	for(my $i=0;$i<$self->{NLINES};$i++){
		my ($bgs,$fgs)=$self->line2site_pair($i);
		my $i1=$matrix1->site_pair2line($bgs,$fgs);
		if(defined $i1){
			$self->{VARIANCE}->[$i]+=$matrix1->{VARIANCE}->[$i1];
			if(@osp_cov){
				$self->{VARIANCE}->[$i]+=$osp_cov[$i1];
			}
		}
		my $i2;
		if(defined $matrix2){
			$i2=$matrix2->site_pair2line($fgs,$bgs);
			$self->{VARIANCE}->[$i]+=$matrix2->{VARIANCE}->[$i2] if defined $i2;
		}else{
			$i2=$matrix1->site_pair2line($fgs,$bgs);
			if(defined $i2){
				$self->{VARIANCE}->[$i]+=$matrix1->{VARIANCE}->[$i2];
				if(@osp_cov){
					$self->{VARIANCE}->[$i]+=$osp_cov[$i2];
				}
			}
		}
	}
}

sub _init_copy{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in SymZScore::_init_copy(): undefined matrix argument!" unless defined $sp_matrix;
	$self->SUPER::_init_copy($sp_matrix,@_);
	$self->{VARIANCE}=[];
	@{$self->{VARIANCE}}=@{$sp_matrix->{VARIANCE}};
}

sub _init_transp{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in SymZScore::_init_transp(): undefined matrix argument!" unless defined $sp_matrix;
	$self->SUPER::_init_transp($sp_matrix,@_);
	$self->{VARIANCE}=[];
	for(my $i=0;$i<$self->{NLINES};$i++){
		my ($bg_idx,$fg_idx)=$self->line2sites_idxs($i);
		my $ti=$sp_matrix->site_idx_pair2line($fg_idx,$bg_idx);
		$self->{VARIANCE}->[$i]=$sp_matrix->{VARIANCE}->[$ti];
	}
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
		if($self->{VARIANCE}->[$i]){
			$stat[$i]/=sqrt($self->{VARIANCE}->[$i]);
			$M=abs($stat[$i]) if(abs($stat[$i])>$M);
		}else{
			$stat[$i]=undef;
		}
	}
	$M=$norm_const if defined $norm_const;
	if($M){
		for(my $i=0;$i<@stat;$i++){
			$stat[$i]/=$M if defined $stat[$i];
		}
	}
	return @stat;
}

1;
	