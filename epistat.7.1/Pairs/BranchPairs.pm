package Pairs::BranchPairs;
use Bio::Phylo::Forest::Tree;
use SitePairInfo;
use Pairs::Stat qw(stat_func);
use strict;

sub find{
	my ($f_intragene,$tree,$f_ignore_terminals,$hr_bg_subst_map,$hr_tg_subst_map)=@_;
	my %site_pairs;
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
					foreach my $j(keys %tg_subst){
						next if($f_intragene&&$i==$j);
						my $str="$i,$j";
						if(!defined $site_pairs{$str}){
							my $spi=SitePairInfo->new();
							$spi->bg_site($i);
							$spi->tg_site($j);
							$site_pairs{$str}=$spi;
						};
						my $p=1;
						$p=0.5 if $f_intragene;
						SitePairInfo::add_subst_info($site_pairs{$str},$node,$node,$p);
					}
				}
			}
		}
	);
	return %site_pairs;
}

sub calc_statistics{
	my ($site_pair_info,$tau,$hr_bg_subst_map,$hr_tg_subst_map)=@_;
	die "\nConsequtivePairs::calc_statistics() Error: Both or none of substitution maps must be provided!" unless (!defined($hr_bg_subst_map)&&!defined($hr_tg_subst_map))||
		(defined($hr_bg_subst_map)&&defined($hr_tg_subst_map));
	my $epistat=0;
	my %bgr_epistat;
	my $bg_site=$site_pair_info->bg_site;
	my $tg_site=$site_pair_info->tg_site;
	foreach my $si(@{$site_pair_info->subst_info}){
		my ($bg_node,$tg_node)=($si->bg_branch,$si->tg_branch);
		my ($w1,$w2);
		if(defined($hr_bg_subst_map)&&defined($hr_tg_subst_map)){
			$w1=$hr_bg_subst_map->{$bg_node->get_name}->{$bg_site}->weight;
			$w2=$hr_tg_subst_map->{$tg_node->get_name}->{$tg_site}->weight;
			unless(defined($w1)&&defined($w2)){
				print STDERR "\nUndefined mutation weight: ";
				if(!defined $w1){
					print STDERR $bg_node->get_name." $bg_site!";
				}
				if(!defined $w2){
					print STDERR $tg_node->get_name." $tg_site!";
				}
				exit 1;
			}
		}
		die "\nError Pairs::BranchPairs::calc_statistics(): Events on different branches!" unless $bg_node==$tg_node;
		my $esi=stat_func(sub {return exp(-$_[0]/$_[1]) if defined($_[1]); return 1;},$si,$tau,$w1,$w2);
		$epistat+=$esi;
	}
	$site_pair_info->epistat($epistat);
	return $epistat;
}

1;