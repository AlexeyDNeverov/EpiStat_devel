package ParallelSubstWeights;
use Bio::Phylo::Forest::Tree;
use Class::Struct;
use strict;

sub is_confident_branch{
	my ($rh_branch_supports,$brname,$threshold)=@_;
	return $rh_branch_supports->{$brname}>=$threshold;
}

#struct SubstInfo => {
#	site => '$',
#	bases => '@'
#};

sub unparse_subst_abbr{
	my $si=shift;
	my @bases=@{$si->bases};
	return $bases[0].$si->site.$bases[1];
}

sub calculate_mutations_confidences{
	my ($tree,$rh_mutation_map,$rh_branch_supports,$confidence_threshold)=@_;
	my %closest_conf_branch;
	my %mutation_counts;
	my %mutation_weights;
	
	$tree->visit_breadth_first(
		-in => sub{
			my $node=shift;
			my $name=$node->get_name;
			$mutation_counts{$name}={};
			if($node->is_root){
				$closest_conf_branch{$name}=$node;
			}else{
				if((!$node->is_terminal) && is_confident_branch($rh_branch_supports,$name,$confidence_threshold)){
					$closest_conf_branch{$name}=$node;
				}else{
					my $pnode=$node->get_parent;
					my $pname=$pnode->get_name;
					$closest_conf_branch{$name}=$closest_conf_branch{$pname};
				}
			}
		}
	);

	$tree->visit_depth_first(
		-in => sub{
			my $node=shift;
			if(!$node->is_root){
				my $name=$node->get_name;
				my $pnode=$node->get_parent;
				my $pname=$pnode->get_name;
				if($node->is_terminal){
					if(defined $rh_mutation_map->{$name}){
						foreach my $site(keys %{$rh_mutation_map->{$name}}){
							#my $str=unparse_subst_abbr($rh_mutation_map->{$name}->{$site});
							my $str=$site;
							$mutation_counts{$pname}->{$str}=[(0,0,$site)] unless defined $mutation_counts{$pname}->{$str};
							$mutation_counts{$pname}->{$str}->[0]++;
						}
					}
				}else{
					foreach my $str(keys %{$mutation_counts{$name}}){
						my $site=$mutation_counts{$name}->{$str}->[2];
						$mutation_counts{$pname}->{$str}=[(0,0,$site)] unless defined $mutation_counts{$pname}->{$str};
						unless(defined($rh_mutation_map->{$name})&&defined($rh_mutation_map->{$name}->{$site})){
							$mutation_counts{$pname}->{$str}->[0]+=$mutation_counts{$name}->{$str}->[0];
							if(is_confident_branch($rh_branch_supports,$name,$confidence_threshold)){
								$mutation_counts{$pname}->{$str}->[1]+=$mutation_counts{$name}->{$str}->[0];
							}else{
								$mutation_counts{$pname}->{$str}->[1]+=$mutation_counts{$name}->{$str}->[1];
							}
						}
					}
					if(defined $rh_mutation_map->{$name}){
						foreach my $site(keys %{$rh_mutation_map->{$name}}){
							#my $str=unparse_subst_abbr($rh_mutation_map->{$name}->{$site});
							my $str=$site;
							$mutation_counts{$pname}->{$str}=[(0,0,$site)] unless defined $mutation_counts{$pname}->{$str};
							$mutation_counts{$pname}->{$str}->[0]++;
							if(is_confident_branch($rh_branch_supports,$name,$confidence_threshold)){
								$mutation_counts{$pname}->{$str}->[1]++;
							}
						}
					}
				}
			}
		}
	);
	
	$tree->visit_depth_first(
		-in => sub{
			my $node=shift;
			if(!$node->is_root){
				my $name=$node->get_name;
				my $bref=$closest_conf_branch{$name};
				if(defined $rh_mutation_map->{$name}){
					foreach my $site(keys %{$rh_mutation_map->{$name}}){
						#my $str=unparse_subst_abbr($rh_mutation_map->{$name}->{$site});
						my $str=$site;
						if($bref==$node){
							$mutation_weights{$name}->{$site}=1;
						}else{
							my $bref_name=$bref->get_name;
							$mutation_weights{$name}->{$site}=1.0/($mutation_counts{$bref_name}->{$str}->[0]-$mutation_counts{$bref_name}->{$str}->[1]);
						}
					}
				}
			}
		}
	);
	return %mutation_weights;
}

1;

