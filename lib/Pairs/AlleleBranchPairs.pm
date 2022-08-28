package Pairs::AlleleBranchPairs;
@ISA = ("Pairs::BranchPairs","Pairs::AlleleConsecutivePairs",);
use Bio::Phylo::Forest::Tree;
use SitePairInfo;
use EpistatAnalysis::PhyloDistances qw(calc_subst_distance);
use strict;
#constants
use constant {
	ALLELE_BRANCH_PAIRS =>7,
};

#interface
sub new{
	my $class=shift;
	my $self={};
	bless($self,$class);
	my ($f_intragene,$tree,$f_ignore_terminals,$hr_bg_subst_map,$hr_tg_subst_map,$hr_bg_site_idx,$hr_tg_site_idx,$ra_alphabet)=@_;
	$self->Pairs::BranchPairs::_init(@_);
	$self->Pairs::AlleleConsecutivePairs::_init($f_intragene,$tree,$hr_bg_subst_map,$hr_tg_subst_map,$hr_bg_site_idx,$hr_tg_site_idx,$ra_alphabet);
	return $self;
}

sub calc_statistics{
	my $self=shift;
	my ($tau,$stat_type,$hr_bg_subst_map,$hr_tg_subst_map)=@_;
	die "\nAlleleBranchPairs::calc_statistics() Error: unsupported statistics type $stat_type. Allowed value is 7!" unless $stat_type==ALLELE_BRANCH_PAIRS;
	die "\nAlleleBranchPairs::calc_statistics() Error: Both substitution maps must be provided!" unless defined($hr_bg_subst_map)&&defined($hr_tg_subst_map);
	my $bg_site_prof=$self->{BGR_PROFILE};
	my $alpha_size=$bg_site_prof->get_alphabet_size();
	my $tg_site_prof;
	my $f_intragene=$self->{F_INTRAGENE};
	if($f_intragene){
		$tg_site_prof=$bg_site_prof;
	}else{
		$tg_site_prof=$self->{FGR_PROFILE};
	}
	my $hr_site_pairs;
	my $N=1;
	if($f_intragene){
		$hr_site_pairs=$self->{UNORD_PAIRS};
		$N=2;
	}else{
		$hr_site_pairs=$self->{SITE_PAIRS};
	}
	foreach my $p_idx (keys %{$hr_site_pairs}){
		my %allele_site_pairs;
		my @pdistr;
		my $epistat=0;
		my ($s1_idx,$s2_idx)=$self->_line2idx_pair($p_idx);
		for(my $I=0;$I<$N;$I++){
			my $site_pair_info;
			if($f_intragene){
				$site_pair_info=$hr_site_pairs->{$p_idx}->[$I];
				next unless defined $site_pair_info;
			}else{
				$site_pair_info=$hr_site_pairs->{$p_idx};
			}
			my $bg_site=$site_pair_info->bg_site;
			my $tg_site=$site_pair_info->tg_site;
			my %mcounts_all; #counts of mutations in the target site
			my %bg_allele2tg_mutation; #counts of mutations in the target site on the specific background
			my %bg_allele2tg_acounts; #counts ancestral states of mutations in the target site on the specific background
			my %esi;
			foreach my $si(@{$site_pair_info->subst_info}){
				my ($bg_node,$tg_node)=($si->bg_branch,$si->tg_branch);
				my $t=calc_subst_distance($bg_node,$tg_node);
				my $count=$si->count();
				my $bg_node_name=$bg_node->get_name;
				my $tg_node_name=$tg_node->get_name;
				my @base_idxs1=@{$hr_bg_subst_map->{$bg_node->get_name}->{$bg_site}->bases};
				die "\nAlleleBranchPairs::calc_statistics() Error: undefined alleles in the site $bg_site on the branch $bg_node_name" unless @base_idxs1==2;
				my @base_idxs2=@{$hr_tg_subst_map->{$tg_node->get_name}->{$tg_site}->bases};
				die "\nAlleleBranchPairs::calc_statistics() Error: undefined alleles in the site $tg_site on the branch $tg_node_name" unless @base_idxs2==2;
				my $w1=$hr_bg_subst_map->{$bg_node_name}->{$bg_site}->weight;
				my $w2=$hr_tg_subst_map->{$tg_node_name}->{$tg_site}->weight;
				my $w=$count;
				if(defined($w1)&&defined($w2)){
					$w*=$w1*$w2;
				}
				$w*=exp(-$t/$tau) if defined $tau;
				my $alleles=$base_idxs2[0]*$alpha_size+$base_idxs2[1];
				my $tg_anc=$base_idxs2[0];
				my $bg_drv=$base_idxs1[1];
				$bg_allele2tg_mutation{$bg_drv}={} unless defined $bg_allele2tg_mutation{$bg_drv};
				$bg_allele2tg_acounts{$bg_drv}={} unless defined $bg_allele2tg_acounts{$bg_drv};
				$bg_allele2tg_mutation{$bg_drv}->{$alleles}++;
				$bg_allele2tg_acounts{$bg_drv}->{$tg_anc}++;
				unless(defined $mcounts_all{$alleles}){
					if($I==0){
						$mcounts_all{$alleles}=$tg_site_prof->get_mutation_count($tg_anc,$base_idxs2[1],$s2_idx);
					}else{
						$mcounts_all{$alleles}=$tg_site_prof->get_mutation_count($tg_anc,$base_idxs2[1],$s1_idx);
					}
				}
				$esi{$alleles}={} unless defined $esi{$alleles};
				$esi{$alleles}->{$bg_drv}+=$w;
			}
			
			foreach my $bg_allele(keys %bg_allele2tg_mutation){
				foreach my $alleles(keys %{$bg_allele2tg_mutation{$bg_allele}}){
					#my $anc_idx=int($alleles/$alpha_size);
					#my $p=$bg_allele2tg_mutation{$bg_allele}->{$alleles}/($mcounts_all{$alleles}+$bg_allele2tg_acounts{$bg_allele}->{$anc_idx}-$bg_allele2tg_mutation{$bg_allele}->{$alleles});
					#P(s1,A1|b2,s2,B2)
					my $p=$bg_allele2tg_mutation{$bg_allele}->{$alleles}/$mcounts_all{$alleles};
					$epistat+=$p*$esi{$alleles}->{$bg_allele};
				}
			}
		}
		for(my $I=0;$I<$N;$I++){
			my $site_pair_info;
			if($f_intragene){
				$site_pair_info=$hr_site_pairs->{$p_idx}->[$I];
				next unless defined $site_pair_info;
			}else{
				$site_pair_info=$hr_site_pairs->{$p_idx};
			}
			$site_pair_info->epistat($epistat);
		}
	}
	my @sp_info;
	foreach my $str(keys %{$self->{SITE_PAIRS}}){
		push @sp_info, $self->{SITE_PAIRS}->{$str};
	}
	return sort {$a->bg_site<=>$b->bg_site||$a->tg_site<=>$b->tg_site} @sp_info;
}

1;