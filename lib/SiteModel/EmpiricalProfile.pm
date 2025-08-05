package SiteModel::EmpiricalProfile;
#This module computes mutation probabilities from the profile of character states in the site

sub _init{
	my $self=shift;
	my ($tree,$rh_subst_map,$alpha_size,$rh_site2idx,$is_markovian)=@_;
	$self->{ALPHA_SIZE}=$alpha_size;
	$self->{MUT_NUMBER}=[];
	$self->{ALM_NSEQ};
	$self->{ALM_PROFILE}=[];
	$self->{MUT_PROFILE}=[];
	$self->{MUT_COUNTS}=[];
	#$self->{ALL_STATES}=[];
	$self->{DERIVED_STATES}=[];
	$self->{DERIVED_STATE_COUNTS}=[];
	$self->{ANCEST_STATE_COUNTS}=[];
	$self->{ROOT_SEQ}=[];
	$self->{ROOT_RECENT_MUTATIONS}={};
	my $max_site_idx=0;
	my @states;
##DEBUG
my @idx2site;
foreach my $site(keys %{$rh_site2idx}){
	$idx2site[$rh_site2idx->{$site}]=$site;
}
#####
	foreach my $nd_name(keys %{$rh_subst_map}){
		foreach my $site(keys %{$rh_subst_map->{$nd_name}}){
			my $si=$rh_subst_map->{$nd_name}->{$site};
			my $site_idx=$rh_site2idx->{$site};
			$max_site_idx=$site_idx if $max_site_idx<$site_idx;
			my @base_idxs=@{$si->bases};
			if(defined($base_idxs[0])&&defined($base_idxs[1])){
				$self->{MUT_PROFILE}->[$site_idx]={} unless defined $self->{MUT_PROFILE}->[$site_idx];
				$states[$site_idx]={} unless defined $states[$site_idx];
				$self->{DERIVED_STATE_COUNTS}->[$site_idx]={} unless defined $self->{DERIVED_STATE_COUNTS}->[$site_idx];
				$self->{ANCEST_STATE_COUNTS}->[$site_idx]={} unless defined $self->{ANCEST_STATE_COUNTS}->[$site_idx];
				$states[$site_idx]->{$base_idxs[0]}++;
				$states[$site_idx]->{$base_idxs[1]}++;
				$self->{DERIVED_STATE_COUNTS}->[$site_idx]->{$base_idxs[1]}++;
				$self->{ANCEST_STATE_COUNTS}->[$site_idx]->{$base_idxs[0]}++;
				$self->{MUT_PROFILE}->[$site_idx]->{$base_idxs[0]}={} unless defined $self->{MUT_PROFILE}->[$site_idx]->{$base_idxs[0]};
				$self->{MUT_PROFILE}->[$site_idx]->{$base_idxs[0]}->{$base_idxs[1]}++;
				$self->{MUT_NUMBER}->[$site_idx]++;
			}else{
				die "\nSome alleles were undefined: Node=$nd_name\tsite=$site \n\tEdit your input XPARR file!";
			}
		}
	}
	for(my $i=0;$i<=$#states;$i++){
		#$self->{ALL_STATES}->[$i]=[sort {$a<=>$b} keys %{$states[$i]}] if defined $states[$i];
		$self->{DERIVED_STATES}->[$i]=[sort {$a<=>$b} keys %{$self->{DERIVED_STATE_COUNTS}->[$i]}] if defined $self->{DERIVED_STATE_COUNTS}->[$i];
	}
	my %tips_counts;
	my %recent_mutations;
	$tree->visit_depth_first(
		-in   => sub {
			my $node=shift;
			my $name=$node->get_name;
			if($node->is_terminal()){
				$tips_counts{$node->get_name}=1;
			}else{
				$recent_mutations{$name}={};
				foreach my $chnode(@{$node->get_children}){
					my $chname=$chnode->get_name;
					$tips_counts{$name}+=$tips_counts{$chname};
					unless($chnode->is_terminal){
						foreach my $site (keys %{$recent_mutations{$chname}}){
							unless(defined($rh_subst_map->{$chname})&&defined($rh_subst_map->{$chname}->{$site})){
								$recent_mutations{$name}->{$site}=[] unless defined $recent_mutations{$name}->{$site};
								push @{$recent_mutations{$name}->{$site}},@{$recent_mutations{$chname}->{$site}};
							}
						}
					}
					if(defined $rh_subst_map->{$chname}){
						foreach my $site (keys %{$rh_subst_map->{$chname}}){
							$recent_mutations{$name}->{$site}=[] unless defined $recent_mutations{$name}->{$site};
							push @{$recent_mutations{$name}->{$site}},$chname;
						}
					}
				}
			}
		},
		-post   => sub {
			my $node=shift;
			#free allocated memory
			unless($node->is_terminal){
				foreach my $chnode(@{$node->get_children}){
					my $chname=$chnode->get_name;
					my @sites=keys %{$recent_mutations{$chname}};
					foreach my $site (@sites){
						$recent_mutations{$chname}->{$site}=undef;
						delete $recent_mutations{$chname}->{$site};
					}
					delete $recent_mutations{$chname};
				}
			}
		}
	);
	my $rname=$tree->get_root->get_name;
	#get ROOT sequence
	foreach my $site(keys %{$recent_mutations{$rname}}){
		my %alleles;
		my $anc_idx;
		my $site_idx=$rh_site2idx->{$site};
		$self->{ROOT_RECENT_MUTATIONS}->{$site_idx}=[];
		foreach my $nd_name(@{$recent_mutations{$rname}->{$site}}){
			push @{$self->{ROOT_RECENT_MUTATIONS}->{$site_idx}},$nd_name;
			my $si=$rh_subst_map->{$nd_name}->{$site};
			$anc_idx=$si->bases(0);
			die "\nIn the site $site the mutation occurred on the branch $nd_name has undefined uncestral state!" unless defined $anc_idx;
			$alleles{$anc_idx}++;
		}
		unless(scalar(keys %alleles)==1){
			my $n=$alleles{$anc_idx};
			print STDERR "\nWARNING: In the site $site the ambiguous ancestral alleles accounted";
			#print STDERR ":\n\tallele_idx\tcount";
			foreach my $aidx(sort {$a<=>$b} keys %alleles){
				if($n<$alleles{$aidx}){
					$n=$alleles{$aidx};
					$anc_idx=$aidx;
				}
				#print STDERR "\n\t".($aidx+1)."\t".$alleles{$aidx};
			}
		}
		$self->{ROOT_SEQ}->[$site_idx]=$anc_idx;
	}
	#calculate alignment profile
	$tree->visit_breadth_first(
		-in   => sub {
			my $node=shift;
			my $name=$node->get_name;
			if(defined $rh_subst_map->{$name}){
				foreach my $site (keys %{$rh_subst_map->{$name}}){
					my $site_idx=$rh_site2idx->{$site};
					my @base_idxs=@{$rh_subst_map->{$name}->{$site}->bases};
					unless(defined $self->{ALM_PROFILE}->[$site_idx]){
						$self->{ALM_PROFILE}->[$site_idx]={};
						$self->{ALM_PROFILE}->[$site_idx]->{$base_idxs[0]}=$tips_counts{$rname};
					}
					$self->{ALM_PROFILE}->[$site_idx]->{$base_idxs[0]}-=$tips_counts{$name};
					$self->{ALM_PROFILE}->[$site_idx]->{$base_idxs[1]}+=$tips_counts{$name};
				}
			}
		}
	);
	$self->{ALM_NSEQ}=$tips_counts{$rname};
	#makes MUT_PROFILE matrix a Markov matrix
	for(my $i=0;$i<=$max_site_idx;$i++){
		next unless defined $self->{ALM_PROFILE}->[$i];
		foreach my $s(keys %{$self->{ALM_PROFILE}->[$i]}){
			$self->{ALM_PROFILE}->[$i]->{$s}/=$self->{ALM_NSEQ};
		}
		$self->{MUT_COUNTS}->[$i]={};
		foreach my $as(keys %{$self->{MUT_PROFILE}->[$i]}){
			$self->{MUT_COUNTS}->[$i]->{$as}={};
			foreach my $ds(keys %{$self->{MUT_PROFILE}->[$i]->{$as}}){
				$self->{MUT_COUNTS}->[$i]->{$as}->{$ds}=$self->{MUT_PROFILE}->[$i]->{$as}->{$ds};
				$self->{MUT_PROFILE}->[$i]->{$as}->{$ds}/=$self->{ANCEST_STATE_COUNTS}->[$i]->{$as};
			}
		}
		if($is_markovian){
			foreach my $ds(keys %{$self->{DERIVED_STATE_COUNTS}->[$i]}){
				unless(defined $self->{ANCEST_STATE_COUNTS}->[$i]->{$ds}){
					#$ds is a dead-end state
					$self->{MUT_PROFILE}->[$i]->{$ds}={};
					%{$self->{MUT_PROFILE}->[$i]->{$ds}}=%{$self->{ALM_PROFILE}->[$i]};
					my $p=$self->{MUT_PROFILE}->[$i]->{$ds}->{$ds};
					delete $self->{MUT_PROFILE}->[$i]->{$ds}->{$ds};
					if(1.0-$p==0){
						warn "\n!!!WARNING Possible error in XPAR file: Constant site - ".$idx2site[$i]."\tderived allele id=".$ds;
						my @anc_states=keys %{$self->{ANCEST_STATE_COUNTS}->[$i]};
						die "\nError in XPAR file: in the constant site $idx2site[$i] the single ancestral state is expected. Found ".scalar(@anc_states)." states!" unless @anc_states==1;
						$self->{MUT_PROFILE}->[$i]->{$ds}->{$anc_states[0]}=1;
					}else{
						foreach my $state(keys %{$self->{MUT_PROFILE}->[$i]->{$ds}}){
							$self->{MUT_PROFILE}->[$i]->{$ds}->{$state}/=(1.0-$p);
						}
					}
				}
			}
		}
	}
}

#interface
sub new{
	my $class=shift;
	my $self={};
	bless($self,$class);
	$self->_init(@_);
	return $self;
}

sub get_root_seq{
	my $self=shift;
	return @{$self->{ROOT_SEQ}};
}

sub get_root_seq_ref{
	my $self=shift;
	return $self->{ROOT_SEQ};
}

sub get_root_recent_mutations{
	my $self=shift;
	my $site_idx=shift;
	return $self->{ROOT_RECENT_MUTATIONS}->{$site_idx}
}

#returns empirical probability of a particular state in the alignment
#if 'del_state_idx' is specified it is excluded from calculations
sub state_prob{
	my $self=shift;
	my ($state_idx,$site_idx)=@_;
	my $ndel=0;
	unless(defined $self->{ALM_PROFILE}->[$site_idx]){
		warn "\nSileModel::EmpiricalProfile::state_prob(): the site $site_idx is conservative!";
		return undef;
	}
	my $p=$self->{ALM_PROFILE}->[$site_idx]->{$state_idx};
	return 0 unless defined $p;
	return $p;
}

#returns a probability to obtain a particular state from any other state, conditioning on that a mutation has occurred
sub derived_state_prob{
	my $self=shift;
	my ($state_idx,$site_idx)=@_;
	unless(defined $self->{DERIVED_STATE_COUNTS}->[$site_idx]){
		warn "\nSileModel::EmpiricalProfile::derived_state_prob(): the site $site_idx is conservative!";
		return undef;
	}
	my $prob=$self->{DERIVED_STATE_COUNTS}->[$site_idx]->{$state_idx}/$self->{MUT_NUMBER}->[$site_idx];
	return $prob;
}

sub get_derived_state_count{
	my $self=shift;
	my ($state_idx,$site_idx)=@_;
	unless(defined $self->{DERIVED_STATE_COUNTS}->[$site_idx]){
		warn "\nSileModel::EmpiricalProfile::get_derived_state_count(): the site $site_idx is conservative!";
		return undef;
	}
	$cnt=$self->{DERIVED_STATE_COUNTS}->[$site_idx]->{$state_idx};
	unless(defined $cnt){
		die "\nSileModel::EmpiricalProfile::get_derived_state_count(): the derived state with the index $state_idx hasn't been observed in the site $site_idx!";
		return undef;
	}
	return $cnt;
}

sub get_mutation_count{
	my $self=shift;
	my ($site_idx)=@_;
	return $self->{MUT_NUMBER}->[$site_idx];
}

sub get_mutation_count{
	my $self=shift;
	my ($state1_idx,$state2_idx,$site_idx)=@_;
	unless(defined $self->{MUT_COUNTS}->[$site_idx]){
		die "\nSileModel::EmpiricalProfile::get_mutation_count(): the site $site_idx is conservative!";
		return undef;
	}
	unless(defined $self->{MUT_COUNTS}->[$site_idx]->{$state1_idx}){
		die "\nSileModel::EmpiricalProfile::get_mutation_count(): the state with the index $state1_idx hasn't been observed in the site $site_idx!";
		return undef;
	}
	
	my $cnt=$self->{MUT_COUNTS}->[$site_idx]->{$state1_idx}->{$state2_idx};
	return 0 unless defined $cnt;
	return $cnt;
}

sub get_ancestral_state_count{
	my $self=shift;
	my ($state_idx,$site_idx)=@_;
	unless(defined $self->{ANCEST_STATE_COUNTS}->[$site_idx]){
		die "\nSileModel::EmpiricalProfile::get_ancestral_state_count(): the site $site_idx is conservative!";
		return undef;
	}
	unless(defined $self->{ANCEST_STATE_COUNTS}->[$site_idx]->{$state_idx}){
		die "\nSileModel::EmpiricalProfile::get_ancestral_state_count(): the ancestral state with the index $state_idx hasn't been observed in the site $site_idx!";
		return undef;
	}
	return $self->{ANCEST_STATE_COUNTS}->[$site_idx]->{$state_idx};
}

#probability of a mutation from the state $state1_idx to $state2_idx, conditioning on that a mutation from the state $state1_idx has occurred
sub state_change_prob{
	my $self=shift;
	my ($state1_idx,$state2_idx,$site_idx)=@_;
	unless(defined $self->{MUT_PROFILE}->[$site_idx]){
		die "\nSileModel::EmpiricalProfile::state_change_prob(): the site $site_idx is conservative!";
		return undef;
	}
	unless(defined $self->{MUT_PROFILE}->[$site_idx]->{$state1_idx}){
		die "\nSileModel::EmpiricalProfile::state_change_prob(): the state with the index $state1_idx hasn't been observed in the site $site_idx!";
		return undef;
	}
	
	my $p=$self->{MUT_PROFILE}->[$site_idx]->{$state1_idx}->{$state2_idx};
	return 0 unless defined $p;
	return $p;
}

sub get_site_profile{
	my $self=shift;
	my ($site_idx,$anc_state_idx)=@_;
	my @tmp=(0) x $self->{ALPHA_SIZE};
	unless(defined $anc_state_idx){
		if(defined $self->{ALM_PROFILE}->[$site_idx]){
			foreach my $state(keys %{$self->{ALM_PROFILE}->[$site_idx]}){
				$tmp[$state]=$self->{ALM_PROFILE}->[$site_idx]->{$state};
			}
		}else{
			die "\nSileModel::EmpiricalProfile::get_site_profile(): wrong index or the site $site_idx is conservative!";
		}
	}else{
		if(defined $self->{MUT_PROFILE}->[$site_idx]){
			if(defined $self->{MUT_PROFILE}->[$site_idx]->{$anc_state_idx}){
				foreach my $state(keys %{$self->{MUT_PROFILE}->[$site_idx]->{$anc_state_idx}}){
					$tmp[$state]=$self->{MUT_PROFILE}->[$site_idx]->{$anc_state_idx}->{$state};
				}
			}else{
				die "\nSileModel::EmpiricalProfile::get_site_profile(): the state with the index $anc_state_idx hasn't been observed in the site $site_idx!";
			}
		}else{
			die "\nSileModel::EmpiricalProfile::get_site_profile(): the site $site_idx is conservative!";
		}
	}
	return @tmp;
}

sub get_states{
	my $self=shift;
	my ($site_idx)=@_;
	unless(defined $self->{ALM_PROFILE}->[$site_idx]){
		#warn "\nSileModel::EmpiricalProfile::get_states(): the site $site_idx is conservative!";
		return ();
	}
	return sort {$a <=> $b} keys %{$self->{ALM_PROFILE}->[$site_idx]};
}

sub get_ancestral_states{
	my $self=shift;
	my ($site_idx)=@_;
	unless(defined $self->{MUT_PROFILE}->[$site_idx]){
		#warn "\nSileModel::EmpiricalProfile::get_ancestral_states(): the site $site_idx is conservative!";
		return ();
	}
	return sort {$a <=> $b} keys %{$self->{MUT_PROFILE}->[$site_idx]};
}

sub get_derived_states{
	my $self=shift;
	my ($site_idx)=@_;
	unless(defined $self->{DERIVED_STATES}->[$site_idx]){
		#warn "\nSileModel::EmpiricalProfile::get_derived_states(): the site $site_idx is conservative!";
		return ();
	}
	return @{$self->{DERIVED_STATES}->[$site_idx]};
}

sub get_alphabet_size{
	my $self=shift;
	return $self->{ALPHA_SIZE};
}

#site indices are numbered from 1.
sub print_align_prof{
	my $self=shift;
	my ($fh,$rh_site2idx)=@_;
	$fh=*STDOUT unless defined $fh;
	my @idx2site;
	print $fh "site_idx";
	if(defined $rh_site2idx){
		my @sites=keys %{$rh_site2idx};
		print $fh "\tsite";
		foreach my $site(@sites){
			my $idx=$rh_site2idx->{$site};
			$idx2site[$idx]=$site;
		}
	}
	print $fh "\tallele_freqs";
	my $nd=int(log($self->{ALM_NSEQ})/log(10));
	$nd+=2;
	my $n=@{$self->{ALM_PROFILE}};
	for(my $idx=0;$idx<$n;$idx++){
		next unless defined $self->{ALM_PROFILE}->[$idx];
		print $fh "\n".($idx+1);
		print $fh "\t".$idx2site[$idx] if defined $rh_site2idx;
		my @states=$self->get_states($idx);
		my $p=0;
		for(my $i=0;$i<$self->{ALPHA_SIZE};$i++){
			if(($i==$states[$p])&&($p<@states)){
				print $fh "\t".sprintf("%.$nd"."f",$self->{ALM_PROFILE}->[$idx]->{$states[$p++]});
			}else{
				print $fh "\t";
			}
		}
	}
}

#site and allele indices are numbered from 1.
sub print_mut_prof{
	my $self=shift;
	my ($fh,$rh_site2idx)=@_;
	$fh=*STDOUT unless defined $fh;
	my @idx2site;
	print $fh "site_idx";
	if(defined $rh_site2idx){
		my @sites=keys %{$rh_site2idx};
		print $fh "\tsite";
		foreach my $site(@sites){
			my $idx=$rh_site2idx->{$site};
			$idx2site[$idx]=$site;
		}
	}
	print $fh "\tanc_allele_idx\tmutation_freqs";
	my $n=@{$self->{MUT_PROFILE}};
	my $nd=int(log($self->{ALM_NSEQ})/log(10));
	$nd+=2;
	for(my $idx=0;$idx<$n;$idx++){
		next unless defined $self->{MUT_PROFILE}->[$idx];
		my @states=$self->get_states($idx);
		foreach my $p(@states){
			print $fh "\n".($idx+1);
			print $fh "\t".$idx2site[$idx] if defined $rh_site2idx;
			print $fh "\t".($p+1);
			my @der_states=sort {$a<=>$b} keys %{$self->{MUT_PROFILE}->[$idx]->{$p}};
			my $q=0;
			for(my $j=0;$j<$self->{ALPHA_SIZE};$j++){
				if(($j==$der_states[$q])&&($q<@der_states)){
					print $fh "\t".sprintf("%.$nd"."f",$self->{MUT_PROFILE}->[$idx]->{$p}->{$der_states[$q++]});
				}else{
					print $fh "\t";
				}
			}
		}
	}
}

1;