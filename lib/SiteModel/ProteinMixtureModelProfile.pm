package SiteModel::ProteinMixtureModelProfile;
#This module computes mutation probabilities from the profile of character states in the site 

#$profile_fn - a file in IQTree *.sitefreq format
#$ra_alphabet - a reference on array with allowed site stattes
sub _init{
	my $self=shift;
	my ($profile_fn,$ra_alphabet)=@_;
	$self->{ALPHA_SIZE}=@{$ra_alphabet};
	$self->{ALPHABET}=[@{$ra_alphabet}];
	$self->{NSITES}=0;
	$self->{PROFILE}=[];
	open INPF,"<$profile_fn" or die "\nUnable to open input file: $profile_fn!";
	while(<INPF>){
		chomp;
		s/\s+$//;
		my @line=split '\s+';
		if(@line){
			my $i=shift @line;
			die "\nThe number of states in the line \n$_\n\t is not equal to the number of allowed states ".$self->{ALPHA_SIZE} unless $self->{ALPHA_SIZE}==@line;
			unless(defined $self->{PROFILE}->[$i-1]){
				my $norm=0;
				for(my $j=0;$j<$self->{ALPHA_SIZE};$j++){
					$norm+=$line[$j];
				}
				$self->{PROFILE}->[$i-1]=[@line];
				for(my $j=0;$j<$self->{ALPHA_SIZE};$j++){
					$self->{PROFILE}->[$i-1]->[$j]/=$norm;
				}
			}else{
				die "\nThe repetitive definition of a profile for the site $i!";
			}
			$self->{NSITES}++;
		}
	}
	close INPF;
	die "\nThe site numbering in the $profile_fn in not consequent!" unless $self->{NSITES}==@{$self->{PROFILE}};
}

#interface
sub new{
	my $class=shift;
	my $self={};
	bless($self,$class);
	$self->_init(@_);
	return $self;
}

#returns a stationary state probability
sub state_prob{
	my $self=shift;
	my ($state_idx,$site_idx)=@_;
	return $self->{PROFILE}->[$site_idx]->[$state_idx];
}

sub get_site_profile{
	my $self=shift;
	my ($site_idx,$anc_state_idx)=@_;
	my @tmp=@{$self->{PROFILE}->[$site_idx]};
	if(defined $anc_state_idx){
		my $u=$tmp[$anc_state_idx];
		for(my $i=0;$i<@tmp;$i++){
			$tmp[$i]/=(1.0-$u);
		}
		$tmp[$anc_state_idx]=0;
	}
	return @tmp;
}

#returns a probability to obtain a particilar state from any other conditioning on that a mutation has occured
sub derived_state_prob{
	my $self=shift;
	my ($state_idx,$site_idx)=@_;
	my $prob=0;
	for(my $i=0;$i<$self->{ALPHA_SIZE};$i++){
		next if $i==$state_idx;
		$prob+=$self->{PROFILE}->[$site_idx]->[$i]/(1.0-$self->{PROFILE}->[$site_idx]->[$i]);
	}
	$prob*=$self->{PROFILE}->[$site_idx]->[$state_idx];
	return $prob;
}

#probability of a mutation from the state $state1_idx to the $state2_idx conditioning on that a mutation from the state $state1_idx has occured
sub state_change_prob{
	my $self=shift;
	my ($state1_idx,$state2_idx,$site_idx)=@_;
	die "\nError state_change_prob(): different states expected!" if $state1_idx==$state2_idx;
	return $self->{PROFILE}->[$site_idx]->[$state2_idx]/(1.0-$self->{PROFILE}->[$site_idx]->[$state1_idx]);
}

sub get_states{
	my $self=shift;
	return @{$self->{ALPHABET}};
}

sub get_states_number{
	my $self=shift;
	return $self->{ALPHA_SIZE};
}

sub get_alphabet_size{
	my $self=shift;
	return $self->{ALPHA_SIZE};
}

sub get_site_num{
	my $self=shift;
	return $self->{NSITES};
}

sub treat_alleles{
	my $self=shift;
	my ($subst_info,$rh_site2idx)=@_;
	my $ra_bases=$subst_info->bases;
	return 1 if defined($ra_bases->[0])&&defined($ra_bases->[1]);
	return 0 unless defined($ra_bases->[0])||defined($ra_bases->[1]);
	my $idx=$rh_site2idx->{$subst_info->site};
	my @prof=$self->get_site_profile($idx);
	my $n=$self->get_alphabet_size;
	my $i;
	my $p=0;
	if(defined $ra_bases->[0]){
		for(my $k=0;$k<$n;$k++){
			next if $ra_bases->[0]==$k;
			if($prof[$k]>$p){
				$p=$prof[$k];
				$i=$k;
			}
		}
		$ra_bases->[1]=$i if $i<$n;
	}else{
		for(my $k=0;$k<$n;$k++){
			next if $ra_bases->[1]==$k;
			$_=$prof[$k]/(1-$prof[$k]);
			if($_>$p){
				$p=$_;
				$i=$k;
			}
		}
		$ra_bases->[0]=$i if $i<$n;
	}
	return 1 if $i<$n;
	return 0;
}

1;