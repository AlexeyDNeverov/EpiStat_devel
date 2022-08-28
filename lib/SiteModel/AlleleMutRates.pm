package SiteModel::AlleleMutRates;
#This module computes mutation probabilities from the profile of character states in the site

sub _init{
	my $self=shift;
	my ($tree,$rh_subst_map,$alpha_size,$rh_site2idx,$is_markovian)=@_;
	$self->{ALPHA_SIZE}=$alpha_size;
	$self->{MUT_RATES}=[];
	$self->{MAX_SITE_IDX}=0;
	$self->{IDX2SITE}=[];
	$self->{ROOT_SEQ}=[];
	my $max_site_idx=0;
	my %subtree_length;
	my %recent_mutations;
	my @mut_counts;
	my @mut_subtree_length;
	my $tree_length=0;
	foreach my $site(keys %{$rh_site2idx}){
		my $i=$rh_site2idx->{$site};
		$self->{IDX2SITE}->[$i]=$site;
		$max_site_idx=$i if $max_site_idx<$i;
	}
	$tree->visit_depth_first(
		-in   => sub {
			my $node=shift;
			$tree_length+=$node->get_branch_length unless $node->is_root;
			my $name=$node->get_name;
			if($node->is_terminal){
				$subtree_length{$name}=0;
			}else{
				$recent_mutations{$name}={};
				foreach my $chnode(@{$node->get_children}){
					my $chname=$chnode->get_name;
					$subtree_length{$name}+=$subtree_length{$chname}+$chnode->get_branch_length;
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
				unless($node->is_root){
					if(defined $rh_subst_map->{$name}){
						foreach my $site (keys %{$rh_subst_map->{$name}}){
							my $si=$rh_subst_map->{$name}->{$site};
							my $site_idx=$rh_site2idx->{$site};
							die "\nUnable to convert the site label $site into the site index!" unless defined $site_idx;
							my $anc_idx=$si->bases(1); #ancestral allele for downstream mutations
							$mut_counts[$site_idx]={} unless defined $mut_counts[$site_idx];
							$mut_subtree_length[$site_idx]={} unless defined $mut_subtree_length[$site_idx];
							my $l=$subtree_length{$name};
							my $n=0;
							foreach my $nd_name(@{$recent_mutations{$name}->{$site}}){
								$l-=$subtree_length{$nd_name};
								$n++;
							}
							$mut_subtree_length[$site_idx]->{$anc_idx}+=$l;
							$mut_counts[$site_idx]->{$anc_idx}+=$n;
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
		foreach my $nd_name(@{$recent_mutations{$rname}->{$site}}){
			my $si=$rh_subst_map->{$nd_name}->{$site};
			$anc_idx=$si->bases(0);
			die "\nIn the site $site the mutation occurred on the branch $nd_name has undefined uncestral state!" unless defined $anc_idx;
			$alleles{$anc_idx}++;
		}
		unless(scalar(keys %alleles)==1){
			my $n=$alleles{$anc_idx};
			print STDERR "\nWARNING: In the site $site the ambiguous ancestral alleles acconted";
			#print STDERR ":\n\tallele_idx\tcount";
			foreach my $aidx(sort {$a<=>$b} keys %alleles){
				if($n<$alleles{$aidx}){
					$n=$alleles{$aidx};
					$anc_idx=$aidx;
				}
				#print STDERR "\n\t".($aidx+1)."\t".$alleles{$aidx};
			}
		}
		my $site_idx=$rh_site2idx->{$site};
		$self->{ROOT_SEQ}->[$site_idx]=$anc_idx;
	}
	foreach my $site (keys %{$rh_site2idx}){
		my $site_idx=$rh_site2idx->{$site};
		die "\nUnable to convert the site label $site into the site index!" unless defined $site_idx;
		my $anc_idx=$self->{ROOT_SEQ}->[$site_idx];
		$mut_counts[$site_idx]={} unless defined $mut_counts[$site_idx];
		$mut_subtree_length[$site_idx]={} unless defined $mut_subtree_length[$site_idx];
		my $l=$subtree_length{$rname};
		my $n=0;
		foreach my $nd_name(@{$recent_mutations{$rname}->{$site}}){
			$l-=$subtree_length{$nd_name};
			$n++;
		}
		$mut_subtree_length[$site_idx]->{$anc_idx}+=$l;
		$mut_counts[$site_idx]->{$anc_idx}+=$n;
	}
	my @site_mut_counts;
	foreach my $name (keys %{$rh_subst_map}){
		if(defined $rh_subst_map->{$name}){
			foreach my $site (keys %{$rh_subst_map->{$name}}){
				my $si=$rh_subst_map->{$name}->{$site};
				my $site_idx=$rh_site2idx->{$site};
				my ($anc_idx,$der_idx)=@{$si->bases};
				$self->{MUT_RATES}->[$site_idx]={} unless defined $self->{MUT_RATES}->[$site_idx];
				$self->{MUT_RATES}->[$site_idx]->{$anc_idx}=0;
				$self->{MUT_RATES}->[$site_idx]->{$der_idx}=0;
				$site_mut_counts[$site_idx]++;
			}
		}
	}
	for(my $i=0;$i<=$max_site_idx;$i++){
		if(defined($self->{MUT_RATES}->[$i])){
			foreach my $a_idx(keys %{$self->{MUT_RATES}->[$i]}){
				#neutral site-specific mutation rate
				$self->{MUT_RATES}->[$i]->{$a_idx}=$site_mut_counts[$i]/$tree_length if $is_markovian;
				#$self->{MUT_RATES}->[$i]->{$a_idx}=$site_mut_counts[$i] if $is_markovian;
				if(defined($mut_counts[$i])&&defined($mut_counts[$i]->{$a_idx})){
					if($mut_counts[$i]->{$a_idx}>=1){
						my $r=$mut_counts[$i]->{$a_idx}/$mut_subtree_length[$i]->{$a_idx};
						#my $r=$mut_counts[$i]->{$a_idx};
						#allele dependent site specific mutation rate
						$self->{MUT_RATES}->[$i]->{$a_idx}=$r;
					}
				}
			}
		}
	}
	$self->{MAX_SITE_IDX}=$max_site_idx;
	$self->{TREE_LENGTH}=$tree_length;
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

sub get_mutation_rate{
	my $self=shift;
	my ($site_idx,$allele_idx)=@_;
	die "\nThe accounted site index is wrong or it is the index of a constant site!" unless defined $self->{MUT_RATES}->[$site_idx];
	die "\nThe allele with index $allele_idx has not been found in the site $site_idx!" unless defined $self->{MUT_RATES}->[$site_idx]->{$allele_idx};
	return $self->{MUT_RATES}->[$site_idx]->{$allele_idx};
}

sub print{
	my $self=shift;
	my $fh=shift;
	$fh=*STDOUT unless defined $fh;
	my $nd=int(log($self->{TREE_LENGTH})/log(10));
	$nd=$nd<0?0:$nd;
	$nd++;
	print $fh "site_idx\tsite\tmut_rates_from_alleles\n";
	for(my $i=0;$i<=$self->{MAX_SITE_IDX};$i++){
		next unless defined $self->{MUT_RATES}->[$i];
		my $site=$self->{IDX2SITE}->[$i];
		print $fh ($i+1)."\t$site";
		for(my $j=0;$j<$self->{ALPHA_SIZE};$j++){
			print $fh "\t";
			print $fh sprintf("%.".$nd."f", $self->{MUT_RATES}->[$i]->{$j}) if defined $self->{MUT_RATES}->[$i]->{$j};
		}
		print $fh "\n";
	}
}

1;