package AssociationStatistics::CorrelationLike;
#Module for calculating covariance/correlation-like association statistics for a site pair on the base of epistatic statistics
use AssociationStatistics::BaseAssociationMeasure;
@ISA = ("AssociationStatistics::BaseAssociationMeasure");

sub _init{
	my $self=shift;
	my ($dataset_desc)=@_;
	$self->SUPER::_init(@_);
	#check sufficiency of the dataset description
	die "\nError in CorrelationLike::_init(): undefined file with pairs of sites!" unless defined $dataset_desc->pairs_fn;
	die "\nError in CorrelationLike::_init(): undefined file with marginal background sites' moments!" unless defined $dataset_desc->bgr_site_moments12_fn;
	die "\nError in CorrelationLike::_init(): undefined file with marginal foreground sites' moments!" unless defined $dataset_desc->fgr_site_moments12_fn;
	###
	$self->{FGR_MUT_NUMBER}=0;
	my $input_fn=$self->{SITE_PAIRS_FN};
	open INPF, "<$input_fn" or die "\nUnable to open input file $input_fn!";
	<INPF>; #skip header
	my $str="";
	my $i=-1;
	my $j=0;
	my $f_intragene=$self->{F_INTRAGENE};
	while(<INPF>){
		chomp;
		$_=~s/\s*$//;
		if($_ ne ""){
			my @line=split '\t';
			if($str ne $line[0]){
				$str=$line[0];
				$i++;
			}
			if($i>0){
				$self->{FGR_MUT_NUMBER}+=$line[3] if($line[3]>0&&$f_intragene);
				last;				
			}elsif($i==0){
				if($line[3]>0){
					$self->{FGR_MUT_NUMBER}+=$line[3];
				}
			}
			$j++;
		}
	}
	close INPF;
	
	my $norm=scalar(@{$self->{BGR_SITES}});
	if($self->{F_INTRAGENE}){
		$norm--;
	}
	$norm=$self->{FGR_MUT_NUMBER}*$norm;
	
	$self->{BGR_SITE_VAR}={};
	$self->{FGR_SITE_VAR}={};
	$input_fn=$dataset_desc->bgr_site_moments12_fn;
	open INPF, "<$input_fn" or die "\nUnable to open input file $input_fn!";
	<INPF>; #skip header
	my $n=0;
	while(<INPF>){
		chomp;
		$_=~s/\s*$//;
		if($_ ne ""){
			my @line=split '\t';
			$self->{BGR_SITE_VAR}->{$line[0]}=$line[2]/($norm**2);
			$n++;
		}
	}
	close INPF;
	$_=$self->get_bgr_site_num;
	die "\nNumber of sites $n in the file $input_fn doesn't match the number of background sites $_!" unless $n==$_;
	
	$input_fn=$dataset_desc->fgr_site_moments12_fn;
	open INPF, "<$input_fn" or die "\nUnable to open input file $input_fn!";
	<INPF>; #skip header
	$n=0;
	while(<INPF>){
		chomp;
		$_=~s/\s*$//;
		if($_ ne ""){
			my @line=split '\t';
			$self->{FGR_SITE_VAR}->{$line[0]}=$line[2]/($norm**2);
			$n++;
		}
	}
	close INPF;
	$_=$self->get_fgr_site_num;
	die "\nNumber of sites $n in the file $input_fn doesn't match the number of foreground sites $_!" unless $n==$_;
	
	#normalization
	$j=0;
	for(my $i=0;$i<$self->{NLINES};$i++){
		$self->{MEAN}->[$i]/=$norm;
		if($self->{OBS_LINES}->[$j]==$i){
			$self->{OBS}->[$j]/=$norm;
			$j++
		}
	}
}

sub _init_copy{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in CorrelationLike::_init_copy(): copy constructor is disabled!";
}

sub _init_transp{
	my $self=shift;
	my $sp_matrix=shift;
	die "\nError in CorrelationLike::_init_transp(): matrix transposition is disabled!";
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
	