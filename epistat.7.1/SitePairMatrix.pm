package SitePairMatrix;
#Interface for site pair matrix
#use List::BinarySearch qw(binsearch);

sub _line2fgr_idx{
	my $self=shift;
	my ($line)=@_;
	my $k;
	my $nfgr=@{$self->{FGR_SITES}};
	if($self->{F_INTRAGENE}){
		$k=$line%($nfgr-1);
		$k++ if $k>=int($line/($nfgr-1));
	}else{
		$k=$line%$nfgr;
	}
	return $k;
}

sub _init{
	my $self=shift;
	$self->{SITE_PAIRS_FN}="";
	$self->{F_INTRAGENE}=undef;
	$self->{BGR_SITES};
	$self->{FGR_SITES}=[];
	$self->{BGR_SITES2INDICES};
	$self->{FGR_SITES2INDICES}={};
	$self->{FGR_MUT_NUMBERS};
	$self->{BGR_MUT_NUMBERS};
	$self->{INFILE_READ_MODE};
	$self->{NLINES};
	my ($pairs_fn,$f_intragene_pairs,$infile_reading_mode)=@_;
	if(defined $f_intragene_pairs){
		$self->{F_INTRAGENE}=$f_intragene_pairs;
	}
	if(defined $pairs_fn){
		$self->{SITE_PAIRS_FN}=$pairs_fn;
	}
	if(!defined $infile_reading_mode){
		$infile_reading_mode=0; #do not read mutations' numbers
	}elsif(!($infile_reading_mode>=0 && $infile_reading_mode<=3)){
		die "\nError in SitePairMatrix::_init(): Unknown value $infile_reading_mode for parameter for the mutations' numbers reading mode !";
	}
	$self->{INFILE_READ_MODE}=$infile_reading_mode;
}

sub _init_sites{
	my $self=shift;
	my ($pairs_fn,$f_intragene_pairs,$infile_reading_mode)=@_;
	if(defined($pairs_fn)&&defined($f_intragene_pairs)&&($pairs_fn ne $self->{SITE_PAIRS_FN})){
		$self->{F_INTRAGENE}=$f_intragene_pairs;
		$self->{SITE_PAIRS_FN}=$pairs_fn;
		$self->{BGR_SITES}=undef;
		@{$self->{FGR_SITES}}=();
		$self->{BGR_MUT_NUMBERS}= undef;
		$self->{FGR_MUT_NUMBERS}= undef;
		$self->{NLINES}=undef;
		if(!defined $infile_reading_mode){
			$infile_reading_mode=0; #do not read mutations' numbers
		}elsif(!($infile_reading_mode>=0 && $infile_reading_mode<=3)){
			die "\nError in SitePairMatrix::_init_sites(): Unknown value $infile_reading_mode for parameter for the mutations' numbers reading mode !";
		}
		$self->{INFILE_READ_MODE}=$infile_reading_mode;
	}
	if(!$self->{NLINES}){
		my $f_intragene=$self->{F_INTRAGENE};
		my $pairs_fn=$self->{SITE_PAIRS_FN};
		open INPF, "<$pairs_fn" or die "\nUnable to open input file $pairs_fn!";
		<INPF>; #skip header
		$self->{BGR_SITES}=[] unless $f_intragene;
		if($self->{INFILE_READ_MODE}==1){
			$self->{BGR_MUT_NUMBERS}=[];
		}elsif($self->{INFILE_READ_MODE}==2){
			$self->{FGR_MUT_NUMBERS}=[];
		}elsif($self->{INFILE_READ_MODE}==3){
			$self->{FGR_MUT_NUMBERS}=[];
			$self->{BGR_MUT_NUMBERS}=[];
		}
		my $str;
		my $i=-1;
		my $j=0;
		my $nfgr;
		if($f_intragene){
			push @{$self->{FGR_SITES}}, undef;
			push @{$self->{FGR_MUT_NUMBERS}},undef if defined $self->{FGR_MUT_NUMBERS};
		}
		while(<INPF>){
			chomp;
			$_=~s/\s+$//;
			if(/\S+/){
				my @line=split "\t";
				if($str ne $line[0]){
					$str=$line[0];
					push @{$self->{BGR_SITES}},$str if defined $self->{BGR_SITES};
					push @{$self->{BGR_MUT_NUMBERS}},$line[2] if defined $self->{BGR_MUT_NUMBERS};
					$i++;
					if($f_intragene && $i>0){
						die "\nRows in the $pairs_fn are incorrectly sorted! Expected rows order is first sorted by background site, second by foreground site" unless $self->{FGR_SITES}->[$i]==$str;
					}
				}
				if($i>0){
					if($f_intragene && (!defined $self->{FGR_SITES}->[0])){
						$self->{FGR_SITES}->[0]=$line[1];
						$self->{FGR_MUT_NUMBERS}->[0]=$line[3] if defined $self->{FGR_MUT_NUMBERS};
					}
					my $k=$self->_line2fgr_idx($j);
					die "\nRows in the $pairs_fn are incorrectly sorted! Expected rows order is first sorted by background site, second by foreground site" unless $self->{FGR_SITES}->[$k]==$line[1];
				}elsif($i==0){
					if($j>0){
						die "\nRows in the $pairs_fn are incorrectly sorted! Expected rows order is first sorted by background site, second by foreground site" unless $self->{FGR_SITES}->[-1]<$line[1];
					}
					push @{$self->{FGR_SITES}},$line[1];
					push @{$self->{FGR_MUT_NUMBERS}},$line[3] if defined $self->{FGR_MUT_NUMBERS};
				}
				$j++;
			}
		}
		$self->{NLINES}=$j;
		close INPF;
		for(my $i=0;$i<@{$self->{FGR_SITES}};$i++){
			my $site=$self->{FGR_SITES}->[$i];
			$self->{FGR_SITES2INDICES}->{$site}=$i;
		}
		if($f_intragene){
			$self->{BGR_SITES}=$self->{FGR_SITES};
			$self->{BGR_SITES2INDICES}=$self->{FGR_SITES2INDICES};
		}else{
			$self->{BGR_SITES2INDICES}={};
			for(my $i=0;$i<@{$self->{BGR_SITES}};$i++){
				my $site=$self->{BGR_SITES}->[$i];
				$self->{BGR_SITES2INDICES}->{$site}=$i;
			}
		}
	}
}

sub _init_copy{
	my $self=shift;
	my ($sp_matrix,$rh_bgr_sites,$rh_fgr_sites)=@_;
	die "\nError in SitePairMatrix::_init_copy(): undefined matrix argument!" unless defined $sp_matrix;
	die "\nError in SitePairMatrix::_init_copy(): a '(rh_bgr|rh_fgr)_sites' parameter couldn't be omitted here!" unless $sp_matrix->{F_INTRAGENE}||
		(defined($rh_bgr_sites)&&defined($rh_fgr_sites))||(!(defined($rh_bgr_sites)||defined($rh_fgr_sites)));
	if(defined($rh_bgr_sites)||defined($rh_fgr_sites)){
		$rh_fgr_sites=$rh_bgr_sites unless defined $rh_fgr_sites;
		$rh_bgr_sites=$rh_fgr_sites unless defined $rh_bgr_sites;
	}
	$self->{F_INTRAGENE}=$sp_matrix->{F_INTRAGENE};
	$self->{SITE_PAIRS_FN}=undef;
	$self->{BGR_SITES}=undef;
	#Different filters for background and foreground sites for an INTRAGENE matrix are allowed
	$self->{BGR_SITES}=[] if (!$self->{F_INTRAGENE})||
				(defined($rh_bgr_sites)&&defined($rh_fgr_sites))&&($rh_bgr_sites!=$rh_fgr_sites);
	@{$self->{FGR_SITES}}=();
	$self->{BGR_MUT_NUMBERS}=[] if defined $sp_matrix->{BGR_MUT_NUMBERS};
	$self->{FGR_MUT_NUMBERS}=[] if defined $sp_matrix->{FGR_MUT_NUMBERS};
	$self->{NLINES}=undef;
	$self->{INFILE_READ_MODE}=$sp_matrix->{INFILE_READ_MODE};
	$self->{BGR_SITES2INDICES};
	$self->{FGR_SITES2INDICES}={};
	if(defined $rh_fgr_sites){
		@{$self->{FGR_SITES}}=keys %{$rh_fgr_sites};
		@{$self->{FGR_SITES}}=sort {$a<=>$b} @{$self->{FGR_SITES}};
	}elsif(defined $sp_matrix->{FGR_SITES}){
		@{$self->{FGR_SITES}}=@{$sp_matrix->{FGR_SITES}};
		@{$self->{FGR_MUT_NUMBERS}}=@{$sp_matrix->{FGR_MUT_NUMBERS}} if defined $sp_matrix->{FGR_MUT_NUMBERS};
	}
	if(defined $self->{BGR_SITES}){
		if(defined $rh_bgr_sites){
			@{$self->{BGR_SITES}}=keys %{$rh_bgr_sites};
			@{$self->{BGR_SITES}}=sort {$a<=>$b} @{$self->{BGR_SITES}};
		}else{
			@{$self->{BGR_SITES}}=@{$sp_matrix->{BGR_SITES}};
			@{$self->{BGR_MUT_NUMBERS}}=@{$sp_matrix->{BGR_MUT_NUMBERS}} if defined $sp_matrix->{BGR_MUT_NUMBERS};
		}
	}else{
		$self->{BGR_SITES}=$self->{FGR_SITES};
	}
	if(defined $rh_fgr_sites){
		if(defined $self->{BGR_MUT_NUMBERS}){
			for(my $k=0;$k<$sp_matrix->get_bgr_site_num;$k++){
				my $site=$sp_matrix->{BGR_SITES}->[$k];
				next unless exists($rh_bgr_sites->{$site});
				push @{$self->{BGR_MUT_NUMBERS}},$sp_matrix->{BGR_MUT_NUMBERS}->[$k];
			}
		}
		if(defined $self->{FGR_MUT_NUMBERS}){
			for(my $k=0;$k<$sp_matrix->get_fgr_site_num;$k++){
				my $site=$sp_matrix->{FGR_SITES}->[$k];
				next unless exists($rh_fgr_sites->{$site});
				push @{$self->{FGR_MUT_NUMBERS}},$sp_matrix->{FGR_MUT_NUMBERS}->[$k] ;
			}
		}
		my $nfg=$self->get_fgr_site_num;
		$nfg-- if $self->{F_INTRAGENE};
		$self->{NLINES}=scalar(@{$self->{BGR_SITES}})*$nfg;
	}else{
		$self->{NLINES}=$sp_matrix->{NLINES};
	}
	foreach my $site(keys %{$sp_matrix->{FGR_SITES2INDICES}}){
		$self->{FGR_SITES2INDICES}->{$site}=$sp_matrix->{FGR_SITES2INDICES}->{$site};
	}
	if($sp_matrix->{FGR_SITES2INDICES} == $sp_matrix->{BGR_SITES2INDICES}){
		$self->{BGR_SITES2INDICES}=$self->{FGR_SITES2INDICES};
	}else{
		$self->{BGR_SITES2INDICES}={};
		foreach my $site(keys %{$sp_matrix->{BGR_SITES2INDICES}}){
			$self->{BGR_SITES2INDICES}->{$site}=$sp_matrix->{BGR_SITES2INDICES}->{$site};
		}
	}
}

sub _init_transp{
	my $self=shift;
	my ($sp_matrix)=@_;
	my $read_mode=$sp_matrix->{INFILE_READ_MODE};
	if($sp_matrix->{INFILE_READ_MODE}==1){
		$read_mode=2;
	}elsif($sp_matrix->{INFILE_READ_MODE}==2){
		$read_mode=1;
	}
	$self->SitePairMatrix::_init(undef,$sp_matrix->{F_INTRAGENE},$read_mode);
	$self->{NLINES}=$sp_matrix->{NLINES};
	@{$self->{FGR_SITES}}=@{$sp_matrix->{BGR_SITES}};
	if($sp_matrix->{BGR_SITES}==$sp_matrix->{FGR_SITES}){
		$self->{BGR_SITES}=$self->{FGR_SITES};
	}else{
		$self->{BGR_SITES}=[];
		@{$self->{BGR_SITES}}=@{$sp_matrix->{FGR_SITES}};
	}
	if(defined $sp_matrix->{FGR_MUT_NUMBERS}){
		$self->{BGR_MUT_NUMBERS}=[];
		@{$self->{BGR_MUT_NUMBERS}}=@{$sp_matrix->{FGR_MUT_NUMBERS}};
	}
	if(defined $sp_matrix->{BGR_MUT_NUMBERS}){
		$self->{FGR_MUT_NUMBERS}=[];
		@{$self->{FGR_MUT_NUMBERS}}=@{$sp_matrix->{BGR_MUT_NUMBERS}};
	}
	$self->{BGR_SITES2INDICES};
	$self->{FGR_SITES2INDICES}={};
	foreach my $site(keys %{$sp_matrix->{BGR_SITES2INDICES}}){
		$self->{FGR_SITES2INDICES}->{$site}=$sp_matrix->{BGR_SITES2INDICES}->{$site};
	}
	if($sp_matrix->{FGR_SITES2INDICES} == $sp_matrix->{BGR_SITES2INDICES}){
		$self->{BGR_SITES2INDICES}=$self->{FGR_SITES2INDICES};
	}else{
		$self->{BGR_SITES2INDICES}={};
		foreach my $site(keys %{$sp_matrix->{FGR_SITES2INDICES}}){
			$self->{BGR_SITES2INDICES}->{$site}=$sp_matrix->{FGR_SITES2INDICES}->{$site};
		}
	}
}

#interface declaration
#constructor
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

sub init_sites{
	my $self=shift;
	$self->_init_sites(@_);
}

sub get_bgr_site_num{
	my $self=shift;
	return scalar(@{$self->{BGR_SITES}});
}

sub get_fgr_site_num{
	my $self=shift;
	return scalar(@{$self->{FGR_SITES}});
}

sub line2site_pair{
	my $self=shift;
	my ($line)=@_;
	my $nfgr=$self->get_fgr_site_num;
	$nfgr-- if($self->{F_INTRAGENE});
	my $bgr_idx=int($line/$nfgr);
	my $fgr_idx=$self->_line2fgr_idx($line);
	die "\nSitePairMatrix::line2site_pair(): Incorrect line index for a site pair matrix!" unless $bgr_idx<@{$self->{BGR_SITES}} && $fgr_idx<@{$self->{FGR_SITES}};;
	return ($self->{BGR_SITES}->[$bgr_idx],$self->{FGR_SITES}->[$fgr_idx]);
}

sub line2sites_idxs{
	my $self=shift;
	my ($line)=@_;
	my $nfgr=$self->get_fgr_site_num;
	$nfgr-- if($self->{F_INTRAGENE});
	my $bgr_idx=int($line/$nfgr);
	my $fgr_idx=$self->_line2fgr_idx($line);
	die "\nSitePairMatrix::line2sites_idxs(): Incorrect line index for a site pair matrix!" unless $bgr_idx<@{$self->{BGR_SITES}} && $fgr_idx<@{$self->{FGR_SITES}};;
	return ($bgr_idx,$fgr_idx);
}

sub site_idx_pair2line{
	my $self=shift;
	my ($i,$j)=@_;
	$n=$self->get_bgr_site_num;
	$m=$self->get_fgr_site_num;
	my $bgr_site=$self->{BGR_SITES}->[$i];
	my $fgr_site=$self->{FGR_SITES}->[$j];
	my $k=$j;
	if($self->{F_INTRAGENE}){
		if($bgr_site==$fgr_site){
			return undef;
		}else{
			$k-- if $j>$i;
			$k+=$i*($m-1);
		}
	}else{
		$k+=$i*$m;
	}
	return $k;
}

sub site_pair2line{
	my $self=shift;
	my ($bgr_site,$fgr_site)=@_;
	#my $bgr_idx=binsearch {$a <=> $b} $bgr_site, @{$self->{BGR_SITES}};
	my $bgr_idx=$self->{BGR_SITES2INDICES}->{$bgr_site};
	die "\nThe site $bgr_site was not found!" unless defined $bgr_idx;
	#my $fgr_idx=binsearch {$a <=> $b} $fgr_site, @{$self->{FGR_SITES}};
	my $fgr_idx=$self->{FGR_SITES2INDICES}->{$fgr_site};
	die "\nThe site $fgr_site was not found!" unless defined $fgr_idx;
	return $self->site_idx_pair2line($bgr_idx,$fgr_idx);
}

sub get_sites{
	my $self=shift;
	my ($ref_bgr_sites,$ref_fgr_sites,$rh_bgr_exclude,$rh_fgr_exclude)=@_;
	my $reftype_bgr=ref $ref_bgr_sites;
	die "\nError SitePairMatrix::get_sites(): Wrong argument. Expects array or hash reference !" unless $reftype_bgr eq "ARRAY" || $reftype_bgr eq "HASH";
	my $reftype_fgr=ref $ref_fgr_sites;
	die "\nError SitePairMatrix::get_sites(): Wrong argument. Expects array or hash reference !" unless $reftype_fgr eq "ARRAY" || $reftype_fgr eq "HASH";
	$self->_init_sites();
	if($reftype_bgr eq "HASH"){
		%{$ref_bgr_sites}=();
	}else{
		@{$ref_bgr_sites}=();
	}
	foreach my $site(@{$self->{BGR_SITES}}){
		if(!(defined($rh_bgr_exclude)&&defined($rh_bgr_exclude->{$site}))){
			if($reftype_bgr eq "HASH"){
				$ref_bgr_sites->{$site}=1;
			}else{
				push @{$ref_bgr_sites},$site;
			}
		}
	}
	if($reftype_fgr eq "HASH"){
		%{$ref_fgr_sites}=();
	}else{
		@{$ref_fgr_sites}=();
	}
	foreach my $site(@{$self->{FGR_SITES}}){
		if(!(defined($rh_fgr_exclude)&&defined($rh_fgr_exclude->{$site}))){
			if($reftype_bgr eq "HASH"){
				$ref_fgr_sites->{$site}=1;
			}else{
				push @{$ref_fgr_sites},$site;
			}
		}
	}
}

sub is_square{
	my $self=shift;
	return $self->{BGR_SITES}==$self->{FGR_SITES};
}

1;