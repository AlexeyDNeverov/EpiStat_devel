package AssociationStatistics::BaseAssociationMeasure;
#Module for calculating covariation/correlation-like association statistics for a site pair on the base of epistatic statistics
use SitePairMatrix;
use AssociationStatistics::NonNegativeValueFilter;
use AssociationStatistics::PvalueFilter;
use AssociationStatistics::PairMutationNumberFilter;
use AssociationStatistics::SiteMutationNumberFilter;
use AssociationStatistics::SelectedSitesFilter;
@ISA = ("SitePairMatrix");

sub _init{
	my $self=shift;
	my ($dataset_desc,$general_settings,$infile_reading_mode)=@_;
	$infile_reading_mode=0 unless defined $infile_reading_mode;
	if($general_settings->f_mutation_numbers eq "sites"){
		if($general_settings->bgr_min_mutations>0){
			$infile_reading_mode++ unless $infile_reading_mode==3;
		}
		if($general_settings->fgr_min_mutations>0){
			$infile_reading_mode+=2 unless $infile_reading_mode>=2;
		}
	}
	$self->SUPER::_init($dataset_desc->pairs_fn,$dataset_desc->f_intragene_pairs,$infile_reading_mode);
	$self->SUPER::_init_sites;
	#check sufficiency of the dataset description
	die "\nError in Associations::BaseAssociationMeasure::_init(): undefined matrix name!" unless defined $dataset_desc->name;
	die "\nError in Associations::BaseAssociationMeasure::_init(): undefined file with epistatic statistics!" unless defined $dataset_desc->epistat_fn;
	die "\nError in Associations::BaseAssociationMeasure::_init(): undefined file with expected epistatic statistics!" unless defined $dataset_desc->mean_fn;
	###
	$self->{NAME}=$dataset_desc->name;
	$self->{OBS}=[];
	$self->{OBS_LINES}=[];
	$self->{MEAN}=[];
	$self->{FILTERS}=[];
	if($general_settings->no_negative_elements){
		my $f=AssociationStatistics::NonNegativeValueFilter->new($self);
		push @{$self->{FILTERS}},$f;
	}
	my $str=$dataset_desc->epistat_fn;
	open INPF, "<$str" or die "\nUnable to open input file $str!";
	my $i=0;
	while(<INPF>){
		chomp;
		$_=~s/\s*$//;
		if($_ ne ""){
			my @lines=split "\t";
			if($lines[0]>0){
				push @{$self->{OBS}},$lines[0];
				push @{$self->{OBS_LINES}},$i;
			}
			$i++;
		}
	}
	close INPF;
	die "\nError in AssociationStatistics::BaseAssociationMeasure::_init(): Unexpected number of elements in $str!" unless $self->{NLINES}==$i;
	my $str=$dataset_desc->mean_fn;
	open INPF, "<$str" or die "\nUnable to open input file $str!";
	while(<INPF>){
		chomp;
		$_=~s/\s*$//;
		if($_ ne ""){
			push @{$self->{MEAN}},$_;
		}
	}
	close INPF;
	if(defined($dataset_desc->lower_pvalue_fn)&&defined($dataset_desc->lower_pvalue_threshold) ||
		defined($dataset_desc->upper_pvalue_fn)&&defined($dataset_desc->upper_pvalue_threshold)){
		my $f=AssociationStatistics::PvalueFilter->new($self,$dataset_desc->lower_pvalue_fn,$dataset_desc->lower_pvalue_threshold,
			$dataset_desc->upper_pvalue_fn,$dataset_desc->upper_pvalue_threshold);
		push @{$self->{FILTERS}},$f;
	}
	if($general_settings->bgr_min_mutations>0||$general_settings->fgr_min_mutations>0){
		my $f;
		if($general_settings->f_mutation_numbers eq "sites"){
			$f=AssociationStatistics::SiteMutationNumberFilter->new($self,$general_settings->bgr_min_mutations,$general_settings->fgr_min_mutations);
		}elsif($general_settings->f_mutation_numbers eq "pairs"){
			$f=AssociationStatistics::PairMutationNumberFilter->new($self,$dataset_desc->epistat_fn,$general_settings->bgr_min_mutations,$general_settings->fgr_min_mutations);
		}
		push @{$self->{FILTERS}},$f;
	}
	if(defined($dataset_desc->bgr_selected_sites_fn)||defined($dataset_desc->bgr_selected_sites_fn)){
		my $f=AssociationStatistics::SelectedSitesFilter->new($self,0,$dataset_desc->bgr_selected_sites_fn,$dataset_desc->fgr_selected_sites_fn);
		push @{$self->{FILTERS}},$f;
	}
}

sub _init_copy{
	my $self=shift;
	#args: BaseAssociationMeasure,hashref,hashref
	#flags: -nofilters - do not copy filters
	my $sp_matrix=shift;
	my $rh_bgr_sites=shift;
	my $rh_fgr_sites=shift;
	my %flags=@_;
	die "\nError in BaseAssociationMeasure::_init_copy(): undefined matrix argument!" unless defined $sp_matrix;
	die "\nError in BaseAssociationMeasure::_init_copy(): a '(rh_bgr|rh_fgr)_sites' parameter couldn't be omitted here!" unless $sp_matrix->{F_INTRAGENE}||
		(defined($rh_bgr_sites)&&defined($rh_fgr_sites))||(!(defined($rh_bgr_sites)||defined($rh_fgr_sites)));
	if(defined($rh_bgr_sites)||defined($rh_fgr_sites)){
		$rh_fgr_sites=$rh_bgr_sites unless defined $rh_fgr_sites;
		$rh_bgr_sites=$rh_fgr_sites unless defined $rh_bgr_sites;
	}
	$self->SUPER::_init_copy($sp_matrix,$rh_bgr_sites,$rh_fgr_sites);
	$self->{NAME}=$sp_matrix->{NAME};
	$self->{OBS}=[];
	$self->{OBS_LINES}=[];
	$self->{MEAN}=[];
	$self->{FILTERS}=[];
	if(defined $rh_fgr_sites){
		my $obs_idx=0;
		my $nline=0;
		for(my $i=0;$i<$sp_matrix->{NLINES};$i++){
			($bgs,$fgs)=$sp_matrix->line2site_pair($i);
			if(exists($rh_bgr_sites->{$bgs})&&exists($rh_fgr_sites->{$fgs})){
				if($sp_matrix->{OBS_LINES}->[$obs_idx]==$i){
					push @{$self->{OBS}},$sp_matrix->{OBS}->[$obs_idx];
					push @{$self->{OBS_LINES}},$nline;
				}
				push @{$self->{MEAN}},$sp_matrix->{MEAN}->[$i];
				$nline++;
			}
			$obs_idx++ if($sp_matrix->{OBS_LINES}->[$obs_idx]==$i);
		}
		die "\nError in copy constructor!" unless $self->{NLINES}==$nline;
	}else{
		@{$self->{OBS}}=@{$sp_matrix->{OBS}};
		@{$self->{OBS_LINES}}=@{$sp_matrix->{OBS_LINES}};
		$self->{MEAN}=$sp_matrix->{MEAN};
	}
	if(!exists($flags{-nofilters})){
		foreach my $f(@{$sp_matrix->{FILTERS}}){
			push @{$self->{FILTERS}}, $f->copy($self);
		}
	}
}

sub _init_transp{
	my $self=shift;
	$self->SUPER::_init_transp(@_);
	#args: BaseAssociationMeasure
	my ($sp_matrix)=@_;
	my $str=$sp_matrix->{NAME};
	$self->{NAME}="tr(".$str.")";
	$self->{OBS}=[];
	$self->{OBS_LINES}=[];
	$self->{MEAN}=[];
	$self->{FILTERS}=[];
	my @tmp_obs;
	my $j=0;
	for(my $i=0;$i<$sp_matrix->{NLINES};$i++){
		($bg_idx,$fg_idx)=$sp_matrix->line2sites_idxs($i);
		my $k=$self->site_idx_pair2line($fg_idx,$bg_idx);
		if(defined $k){
			$self->{MEAN}->[$k]=$sp_matrix->{MEAN}->[$i];
			if($sp_matrix->{OBS_LINES}->[$j]==$i){
				push @tmp_obs,[($k,$sp_matrix->{OBS}->[$j])];
			}
		}
		$j++ if($sp_matrix->{OBS_LINES}->[$j]==$i);
	}
	@tmp_obs=sort {$a->[0]<=>$b->[0]} @tmp_obs;
	foreach my $i(@tmp_obs){
		push @{$self->{OBS_LINES}},$i->[0];
		push @{$self->{OBS}},$i->[1];
	}
	foreach my $f(@{$sp_matrix->{FILTERS}}){
		push @{$self->{FILTERS}}, $f->transpose($self);
	}
}

#interface declaration

sub mtx_diag_value{
	my $self=shift;
	return undef;
}

sub get_statistics{
	my $self=shift;
	my $j=0;
	my @stat;
	for(my $i=0;$i<$self->{NLINES};$i++){
		my $obs=0;
		if($i==$self->{OBS_LINES}->[$j]){
			$obs=$self->{OBS}->[$j++];
		}
		my $val=$obs-$self->{MEAN}->[$i];
		push @stat,$val;
	}
	foreach my $f(@{$self->{FILTERS}}){
		@stat=$f->apply(\@stat);
	}
	return @stat;
}

sub get_sites{
	my $self=shift;
	my $ref_bgr_sites=shift;
	my $ref_fgr_sites=shift;
	#flags: -nofilters - do not apply filters
	my %flags=@_;
	my %bgr_zeros;
	my %fgr_zeros;
	if((!exists($flags{-nofilters}))&&(@{$self->{FILTERS}}>0)){
		#remove sites with zero margins
		my $f_intragene=$self->{F_INTRAGENE};
		my @bgr_margs;
		my @fgr_margs;
		$n=@{$self->{BGR_SITES}};
		$m=@{$self->{FGR_SITES}};
		my @stat=$self->get_statistics;
		for(my $i=0;$i<$n;$i++){
			my $bgr_site=$self->{BGR_SITES}->[$i];
			for(my $j=0;$j<$m;$j++){
				my $fgr_site=$self->{FGR_SITES}->[$j];
				my $k=$j;
				if($bgr_site!=$fgr_site){
					if($f_intragene){
						$k-- if $j>$i;
						$k+=$i*($m-1);
					}else{
						$k+=$i*$m;
					}
					$bgr_margs[$i]+=$stat[$k];
					$fgr_margs[$j]+=$stat[$k];
				}
			}
		}
		for(my $i=0;$i<$n;$i++){
			my $bgr_site=$self->{BGR_SITES}->[$i];
			$bgr_zeros{$bgr_site}=1 if $bgr_margs[$i]==0;
		}
		for(my $i=0;$i<$m;$i++){
			my $fgr_site=$self->{FGR_SITES}->[$i];
			$fgr_zeros{$fgr_site}=1 if $fgr_margs[$i]==0;
		}
	}
	$self->SUPER::get_sites($ref_bgr_sites,$ref_fgr_sites,\%bgr_zeros,\%fgr_zeros);
}

sub print_matrix{
	my $self=shift;
	my ($file_handle,$rh_bgr_site_filter,$rh_fgr_site_filter)=@_;
	my %args=@_[3..$#_];
	my $norm_const=$args{-norm};
	$matrix_name=$self->{NAME};
	my $n;
	if(defined $rh_bgr_site_filter){
		$n=keys %{$rh_bgr_site_filter};
	}else{
		$n=@{$self->{BGR_SITES}};
	}
	my $m;
	if(defined $rh_fgr_site_filter){
		$m=keys %{$rh_fgr_site_filter};
	}else{
		$m=@{$self->{FGR_SITES}};
	}
	print $file_handle "Matrix: $matrix_name";
	print $file_handle "[$n,$m]\nbgr_site\tfgr_site\telement";
	my $f_intragene=$self->{F_INTRAGENE};
	my @stat=$self->get_statistics($norm_const);
	$n=@{$self->{BGR_SITES}} if(defined $rh_bgr_site_filter);
	$m=@{$self->{FGR_SITES}} if(defined $rh_fgr_site_filter);
	for(my $i=0;$i<$n;$i++){
		my $bgr_site=$self->{BGR_SITES}->[$i];
		if(defined $rh_bgr_site_filter){
			next unless defined $rh_bgr_site_filter->{$bgr_site};
		}
		for(my $j=0;$j<$m;$j++){
			my $fgr_site=$self->{FGR_SITES}->[$j];
			if(defined $rh_fgr_site_filter){
				next unless defined $rh_fgr_site_filter->{$fgr_site};
			}
			print $file_handle "\n$bgr_site\t$fgr_site";
			my $k=$self->site_idx_pair2line($i,$j);
			my $str=$self->mtx_diag_value;
			if(defined $k){
				$str=$stat[$k];
			}
			$str=sprintf("%.3f",$str);
			print $file_handle "\t$str";
		}
	}
}

sub get_name{
	my $self=shift;
	return $self->{NAME};
}

1;