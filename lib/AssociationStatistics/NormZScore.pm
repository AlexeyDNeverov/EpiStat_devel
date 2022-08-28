package AssociationStatistics::NormZScore;
#Module for calculating covariation/correlation-like association statistics for a site pair on the base of epistatic statistics
use AssociationStatistics::BaseAssociationMeasure;
@ISA = ("AssociationStatistics::BaseAssociationMeasure");

sub _init{
	my $self=shift;
	my ($dataset_desc,$rh_site_pair_filters)=@_;
	$self->SUPER::_init($dataset_desc,$rh_site_pair_filters,2);
	#check sufficiency of the dataset description
	die "\nError in NormZScore::_init(): undefined file with pairs of sites!" unless defined $dataset_desc->variance_fn;
	die "\nError in NormZScore::_init(): undefined file with pairs of sites!" unless defined $dataset_desc->fgr_site_moments12_fn;
	###
	$self->{VARIANCE}=[];
	my $input_fn=$dataset_desc->variance_fn;
	open INPF, "<$input_fn" or die "\nUnable to open input file $input_fn!";
	while(<INPF>){
		chomp;
		$_=~s/\s*$//;
		if($_ ne ""){
			push @{$self->{VARIANCE}},$_;
		}
	}
	close INPF;
	die "\nNumber of lines in the file $input_fn isn't equal to number of lines in files with epistatic statistics!" unless @{$self->{VARIANCE}}==$self->{NLINES};
	$self->{FGR_SITE_MEAN}={};
	$self->{FGR_SITE_VAR}={};
	$input_fn=$dataset_desc->fgr_site_moments12_fn;
	open INPF, "<$input_fn" or die "\nUnable to open input file $input_fn!";
	<INPF>; #skip header
	my $norm=scalar(@{$self->{BGR_SITES}});
	if($self->{F_INTRAGENE}){
		$norm--;
	}
	my $n=0;
	while(<INPF>){
		chomp;
		$_=~s/\s*$//;
		if($_ ne ""){
			my @line=split '\t';
			$self->{FGR_SITE_MEAN}->{$line[0]}=$line[1]/$norm;
			$self->{FGR_SITE_VAR}->{$line[0]}=$line[2]/$norm;
			$n++;
		}
	}
	close INPF;
	$_=$self->get_fgr_site_num;
	die "\nNumber of sites $n in the file $input_fn doesn't match the number of foreground sites $_!" unless $n==$_;
}

sub _init_copy{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in NormZScore::_init_copy(): copy constructor is disabled!";
}

sub _init_transp{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in NormZScore::_init_transp(): matrix transposition is disabled!";
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
		if($self->{VARIANCE}->[$i]){
			my ($bgs,$fgs)=$self->line2site_pair($i);
			my $fgi=$self->_line2fgr_idx($i);
			$stat[$i]/=($self->{FGR_MUT_NUMBERS}->[$fgi]-$self->{FGR_SITE_MEAN}->{$fgs})*sqrt($self->{VARIANCE}->[$i]/$self->{FGR_SITE_VAR}->{$fgs});
			$stat[$i]=1 if  $stat[$i]>1;
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
	