package Pairs::AlleleConsecutivePairs;	#idea 82
@ISA = ("Pairs::ConsecutivePairs");
use Bio::Phylo::Forest::Tree;
use SitePairInfo;
use EpistatAnalysis::PhyloDistances qw(calc_subst_distance);
use SiteModel::EmpiricalProfile;
use strict;
#constants
use constant {
	ALL_ALLELE_PAIRS_A1 =>3,
	ALL_ALLELE_PAIRS_B2 =>4,
	INDEP_ALLELE_PAIRS_A1 =>5,
	INDEP_ALLELE_PAIRS_B2 =>6,
};

sub _init{
	my $self=shift;
	my ($f_intragene,$tree,$hr_bg_subst_map,$hr_tg_subst_map,$hr_bg_site_idx,$hr_tg_site_idx,$ra_alphabet)=@_;
	my $alpha_size=@{$ra_alphabet};
	$self->{F_INTRAGENE}=$f_intragene;
	$self->{BGR_PROFILE}=SiteModel::EmpiricalProfile->new($tree,$hr_bg_subst_map,$alpha_size,$hr_bg_site_idx,1);
	$self->{FGR_PROFILE};
	if($f_intragene){
		if($hr_bg_site_idx!=$hr_tg_site_idx){
			foreach my $site(keys %{$hr_bg_site_idx}){
				die "\nError Pairs::AlleleConsequtivePairs::_init(): different site indexation is not allowed for intragenic pairs! site=$site ".$hr_bg_site_idx->{$site}."\t".$hr_tg_site_idx->{$site}
					unless defined($hr_tg_site_idx->{$site})&&($hr_bg_site_idx->{$site}==$hr_tg_site_idx->{$site});
			}
			foreach my $site(keys %{$hr_tg_site_idx}){
				die "\nError Pairs::AlleleConsequtivePairs::_init(): different site indexation is not allowed for intragenic pairs! site=$site" 
					unless defined($hr_bg_site_idx->{$site});
			}
		}else{
			$self->{FGR_PROFILE}=$self->{BGR_PROFILE};
		}
	}
	$self->{FGR_PROFILE}=SiteModel::EmpiricalProfile->new($tree,$hr_tg_subst_map,$alpha_size,$hr_tg_site_idx,1) unless defined $self->{FGR_PROFILE};
	$self->{UNORD_PAIRS};
	if($f_intragene){
		$self->{UNORD_PAIRS}={};
		foreach my $p_idx (keys %{$self->{SITE_PAIRS}}){
			my ($s1_idx,$s2_idx)=$self->_line2idx_pair($p_idx);
			if($s1_idx<$s2_idx){
				$self->{UNORD_PAIRS}->{$p_idx}=[] unless defined $self->{UNORD_PAIRS}->{$p_idx};
				$self->{UNORD_PAIRS}->{$p_idx}->[0]=$self->{SITE_PAIRS}->{$p_idx};
			}else{
				my $cp_idx=$self->_idx_pair2line($s2_idx,$s1_idx);
				$self->{UNORD_PAIRS}->{$cp_idx}=[] unless defined $self->{UNORD_PAIRS}->{$cp_idx};
				$self->{UNORD_PAIRS}->{$cp_idx}->[1]=$self->{SITE_PAIRS}->{$p_idx};
			}
		}
	}
}

#interface
sub new{
	my $class=shift;
	my $self={};
	bless($self,$class);
	my ($f_intragene,$tree,$f_ignore_terminals,$hr_bg_subst_map,$hr_tg_subst_map,$hr_bg_site_idx,$hr_tg_site_idx,$ra_bg_site2cord,$f_branch_pairs,$ra_alphabet)=@_;
	$self->SUPER::_init(@_);
	$self->_init($f_intragene,$tree,$hr_bg_subst_map,$hr_tg_subst_map,$hr_bg_site_idx,$hr_tg_site_idx,$ra_alphabet);
	return $self;
}

sub _account_for_recent_tg_mutations{
	my $self=shift;
	my ($bg_site_idx,$tg_site_idx,$tg_site,$hr_tg_subst_map,$rh_tg_mutation_nodes,$rh_bg_allele2tg_mutation,$rh_mcounts_all,$rh_acounts_all)=@_;
	my $tg_site_prof=$self->{FGR_PROFILE};
	my $bg_site_prof=$self->{BGR_PROFILE};
	my $alpha_size=$tg_site_prof->get_alphabet_size();
	foreach my $tg_node_name(@{$tg_site_prof->get_root_recent_mutations($tg_site_idx)}){
		unless(defined $rh_tg_mutation_nodes->{$tg_node_name}){
			my $bg_allele_idx=$bg_site_prof->get_root_seq_ref->[$bg_site_idx];
			die "\nThere isn't a mutation in the foreground site $tg_site on the branch $tg_node_name!" unless defined($hr_tg_subst_map->{$tg_node_name}->{$tg_site});
			my @base_idxs2=@{$hr_tg_subst_map->{$tg_node_name}->{$tg_site}->bases};
			die "\n\nAlleleConsequtivePairs::calc_statistics() Error: undefined alleles in the site $tg_site on the branch $tg_node_name" unless @base_idxs2==2;
			my $tg_anc=$base_idxs2[0];
			my $alleles=$base_idxs2[0]*$alpha_size+$base_idxs2[1];
			$rh_bg_allele2tg_mutation->{$bg_allele_idx}->{$alleles}++;
			unless(defined $rh_mcounts_all->{$alleles}){
				$rh_mcounts_all->{$alleles}=$tg_site_prof->get_mutation_count($tg_anc,$base_idxs2[1],$tg_site_idx);
			}
			unless(defined $rh_acounts_all->{$tg_anc}){
				$rh_acounts_all->{$tg_anc}=$tg_site_prof->get_ancestral_state_count($tg_anc,$tg_site_idx);
			}
		}
	}
}

sub calc_statistics{
	my $self=shift;
	my ($tau,$stat_type,$hr_bg_subst_map,$hr_tg_subst_map)=@_;
	die "\nAlleleConsequtivePairs::calc_statistics() Error: unsupported statistics type $stat_type. Allowed values from 3 to 6" 
		unless $stat_type==ALL_ALLELE_PAIRS_A1||$stat_type==ALL_ALLELE_PAIRS_B2||
		$stat_type==INDEP_ALLELE_PAIRS_A1||$stat_type==INDEP_ALLELE_PAIRS_B2;
	my $f_independent_pairs=($stat_type==INDEP_ALLELE_PAIRS_A1||$stat_type==INDEP_ALLELE_PAIRS_B2);
	my $f_A1=($stat_type==ALL_ALLELE_PAIRS_A1||$stat_type==INDEP_ALLELE_PAIRS_A1);
	die "\nAlleleConsequtivePairs::calc_statistics() Error: Both substitution maps must be provided!" unless defined($hr_bg_subst_map)&&defined($hr_tg_subst_map);
	my $bg_site_prof=$self->{BGR_PROFILE};
	my $alpha_size=$bg_site_prof->get_alphabet_size();
	my $tg_site_prof=$self->{FGR_PROFILE};
	my $f_intragene=$self->{F_INTRAGENE};
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
			my ($bg_site_idx,$tg_site_idx);
			if($I==0){
				($bg_site_idx,$tg_site_idx)=($s1_idx,$s2_idx);
			}else{
				#Note that site indices correspond to the first ordered pair (first site cordinate less than the second),
				# thus if $I==1 the foreground site has the $s1_idx site index
				($bg_site_idx,$tg_site_idx)=($s2_idx,$s1_idx);
			}
			my %mcounts_all; #counts of mutations in the target site
			my %acounts_all; #counts of ancestral states of mutations in the target site
			my %bg_allele2tg_mutation; #counts of mutations in the target site on the specific background
			my %bg_allele2tg_acounts; #counts ancestral states of mutations in the target site on the specific background
			my %bg_node2esi;
			my %esi;
			my %tg_mutation_nodes;
			foreach my $si(@{$site_pair_info->subst_info}){
				my ($bg_node,$tg_node)=($si->bg_branch,$si->tg_branch);
				my $t=calc_subst_distance($bg_node,$tg_node);
				my $count=$si->count();
				my $bg_node_name=$bg_node->get_name;
				my $tg_node_name=$tg_node->get_name;
				$tg_mutation_nodes{$tg_node_name}=1;
				my @base_idxs1=@{$hr_bg_subst_map->{$bg_node_name}->{$bg_site}->bases};
				die "\n\nAlleleConsequtivePairs::calc_statistics() Error: undefined alleles in the site $bg_site on the branch $bg_node_name" unless @base_idxs1==2;
				my @base_idxs2=@{$hr_tg_subst_map->{$tg_node_name}->{$tg_site}->bases};
				die "\n\nAlleleConsequtivePairs::calc_statistics() Error: undefined alleles in the site $tg_site on the branch $tg_node_name" unless @base_idxs2==2;
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
				$bg_allele2tg_acounts{$bg_drv}={} unless defined  $bg_allele2tg_acounts{$bg_drv};
				$bg_allele2tg_mutation{$bg_drv}->{$alleles}++;
				$bg_allele2tg_acounts{$bg_drv}->{$tg_anc}++;
				unless(defined $acounts_all{$tg_anc}){
					$acounts_all{$tg_anc}=$tg_site_prof->get_ancestral_state_count($tg_anc,$tg_site_idx);
				}
				unless(defined $mcounts_all{$alleles}){
					$mcounts_all{$alleles}=$tg_site_prof->get_mutation_count($tg_anc,$base_idxs2[1],$tg_site_idx);
				}
				$esi{$alleles}={} unless defined $esi{$alleles};
				if($f_independent_pairs){
					$bg_node2esi{$bg_node_name}=[({},0)] unless defined $bg_node2esi{$bg_node_name};
					$bg_node2esi{$bg_node_name}->[0]->{$alleles}={} unless defined $bg_node2esi{$bg_node_name}->[0]->{$alleles};
					$bg_node2esi{$bg_node_name}->[0]->{$alleles}->{$bg_drv}+=$w;
					$bg_node2esi{$bg_node_name}->[1]++;
				}else{
					$esi{$alleles}->{$bg_drv}+=$w;
				}
			}
			if($f_independent_pairs){
				foreach my $name(keys %bg_node2esi){
					foreach my $alleles(keys %{$bg_node2esi{$name}->[0]}){
						foreach my $bg_drv(keys %{$bg_node2esi{$name}->[0]->{$alleles}}){
							$esi{$alleles}->{$bg_drv}+=$bg_node2esi{$name}->[0]->{$alleles}->{$bg_drv}/$bg_node2esi{$name}->[1];
						}
					}
				}
			}
			$self->_account_for_recent_tg_mutations($bg_site_idx,$tg_site_idx,$tg_site,$hr_tg_subst_map,\%tg_mutation_nodes,\%bg_allele2tg_mutation,\%mcounts_all,\%acounts_all);
			foreach my $bg_allele(keys %bg_allele2tg_mutation){
				foreach my $alleles(keys %{$bg_allele2tg_mutation{$bg_allele}}){
					if($f_A1){
						#P((s1,A1)|(b2,s2,B2))
						my $p=$bg_allele2tg_mutation{$bg_allele}->{$alleles}/$mcounts_all{$alleles};
#DEBUG
die "\nMutation count error. Check code of _account_for_recent_tg_mutations()!" unless $p>=0&&$p<=1.0;
##
						$epistat+=$p*$esi{$alleles}->{$bg_allele};
					}else{
						my $anc_idx=int($alleles/$alpha_size);
						#P(B2|(s1,A1),(b2,s2))
						my $p_bg=$bg_allele2tg_mutation{$bg_allele}->{$alleles}/$bg_allele2tg_acounts{$bg_allele}->{$anc_idx};
						my $p_other=0;
						if($acounts_all{$anc_idx}-$bg_allele2tg_acounts{$bg_allele}->{$anc_idx}){
							#P(!B2|!(s1,A1),(b2,s2))
							$p_other=1-($mcounts_all{$alleles}-$bg_allele2tg_mutation{$bg_allele}->{$alleles})/($acounts_all{$anc_idx}-$bg_allele2tg_acounts{$bg_allele}->{$anc_idx});
						}
#DEBUG
die "\nMutation count error. Check code of _account_for_recent_tg_mutations()!" unless $p_bg>=0&&$p_bg<=1.0&&$p_other>=0&&$p_other<=1.0;
##
						$epistat+=($p_bg+$p_other)*$esi{$alleles}->{$bg_allele}/2;
					}
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

sub print_align_profile{
	my $self=shift;
	my ($fh,$seg_no,$hr_site2idx)=@_;
	die "\nError print_algn_profile(): Undefined segnment number 1 - bgr, 2 - fgr!" unless $seg_no==1||$seg_no==2;
	if($seg_no==1){
		$self->{BGR_PROFILE}->print_align_prof($fh,$hr_site2idx);
	}else{
		$self->{FGR_PROFILE}->print_align_prof($fh,$hr_site2idx);
	}
}

sub print_mutation_profile{
	my $self=shift;
	my ($fh,$seg_no,$hr_site2idx)=@_;
	die "\nError print_mut_profile(): Undefined segnment number 1 - bgr, 2 - fgr!" unless $seg_no==1||$seg_no==2;
	if($seg_no==1){
		$self->{BGR_PROFILE}->print_mut_prof($fh,$hr_site2idx);
	}else{
		$self->{FGR_PROFILE}->print_mut_prof($fh,$hr_site2idx);
	}
}

1;