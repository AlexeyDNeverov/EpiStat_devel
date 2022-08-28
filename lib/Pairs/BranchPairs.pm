package Pairs::BranchPairs;
use Bio::Phylo::Forest::Tree;
use SitePairInfo;
use EpistatAnalysis::PhyloDistances qw(calc_subst_distance);
use strict;
use constant BRANCH_PAIRS =>2;

sub _idx_pair2line{
	my $self=shift;
	my ($bg_site_idx,$tg_site_idx)=@_;
	return $bg_site_idx*$self->{FGR_NSITES}+$tg_site_idx;
}

sub _line2idx_pair{
	my $self=shift;
	my ($line)=@_;
	return (int($line/$self->{FGR_NSITES}),$line%$self->{FGR_NSITES});
}

sub _compl_line{
	my $self=shift;
	die "\nError in Pairs::ConsequtivePairs::_compl_line(): This function is appliable to symmetric matrices only!" unless $self->{FGR_NSITES}==$self->{BGR_NSITES};
	my ($line)=@_;
	my ($i,$j)=(int($line/$self->{FGR_NSITES}),$line%$self->{FGR_NSITES});
	return $j*$self->{FGR_NSITES}+$i;
}

sub _init{
	my $self=shift;
	my ($f_intragene,$tree,$f_ignore_terminals,$hr_bg_subst_map,$hr_tg_subst_map,$hr_bg_site_idx,$hr_tg_site_idx)=@_;
	$self->{SITE_PAIRS}={};
	$self->{BGR_NSITES}=keys %{$hr_bg_site_idx};
	$self->{FGR_NSITES}=keys %{$hr_tg_site_idx};
	$tree->visit_breadth_first(
		-in   => sub {
			my $node=shift;
			if(!($node->is_root||($f_ignore_terminals&&$node->is_terminal))){
				my %bg_subst;
				if(defined $hr_bg_subst_map->{$node->get_name}){
					foreach my $i(keys %{$hr_bg_subst_map->{$node->get_name}}){
						$bg_subst{$i}=1;
					}
				}
				my %tg_subst;
				if(defined $hr_tg_subst_map->{$node->get_name}){
					foreach my $i(keys %{$hr_tg_subst_map->{$node->get_name}}){
						$tg_subst{$i}=1;
					}
				}
				foreach my $i(keys %bg_subst){
					my $bg_site_idx=$hr_bg_site_idx->{$i};
					foreach my $j(keys %tg_subst){
						next if($f_intragene&&$i==$j);
						my $tg_site_idx=$hr_tg_site_idx->{$j};
						my $p_idx=$self->_idx_pair2line($bg_site_idx,$tg_site_idx);
						if(!defined $self->{SITE_PAIRS}->{$p_idx}){
							my $spi=SitePairInfo->new();
							$spi->bg_site($i);
							$spi->tg_site($j);
							$self->{SITE_PAIRS}->{$p_idx}=$spi;
						};
						my $p=1;
						$p=0.5 if $f_intragene;
						SitePairInfo::add_subst_info($self->{SITE_PAIRS}->{$p_idx},$node,$node,$p);
					}
				}
			}
		}
	);
}

#interface
sub new{
	my $class=shift;
	my $self={};
	bless($self,$class);
	$self->_init(@_);
	return $self;
}

sub calc_statistics{
	my $self=shift;
	my ($tau,$hr_bg_subst_map,$hr_tg_subst_map)=@_;
	die "\nConsequtivePairs::calc_statistics() Error: Both or none of substitution maps must be provided!" unless (!defined($hr_bg_subst_map)&&!defined($hr_tg_subst_map))||
		(defined($hr_bg_subst_map)&&defined($hr_tg_subst_map));
	foreach my $p_idx (keys %{$self->{SITE_PAIRS}}){
		my $site_pair_info=$self->{SITE_PAIRS}->{$p_idx};
		my $epistat=0;
		my %bgr_epistat;
		my $bg_site=$site_pair_info->bg_site;
		my $tg_site=$site_pair_info->tg_site;
		foreach my $si(@{$site_pair_info->subst_info}){
			my ($bg_node,$tg_node)=($si->bg_branch,$si->tg_branch);
			die "\nError Pairs::BranchPairs::calc_statistics(): Events on different branches!" unless $bg_node==$tg_node;
			my $t=calc_subst_distance($bg_node,$tg_node);
			my $count=$si->count();
			my $esi=$count;
			$esi*=exp(-$t/$tau) if defined $tau;
			if(defined($hr_bg_subst_map)&&defined($hr_tg_subst_map)){
				my $w1=$hr_bg_subst_map->{$bg_node->get_name}->{$bg_site}->weight;
				my $w2=$hr_tg_subst_map->{$tg_node->get_name}->{$tg_site}->weight;
				if(defined($w1)&&defined($w2)){
					$esi*=$w1*$w2;
				}else{
					print STDERR "\nUndefined mutation weight: ";
					if(!defined $w1){
						print STDERR $bg_node->get_name." $bg_site!";
					}else{
						print STDERR $tg_node->get_name." $tg_site!";
					}
					exit 1;
				}
			}
			$epistat+=$esi;
		}
		$site_pair_info->epistat($epistat);
	}
	my @sp_info;
	foreach my $p_idx(keys %{$self->{SITE_PAIRS}}){
		push @sp_info, $self->{SITE_PAIRS}->{$p_idx};
	}
	@sp_info=sort {$a->bg_site<=>$b->bg_site||$a->tg_site<=>$b->tg_site} @sp_info;
	return @sp_info;
}

1;