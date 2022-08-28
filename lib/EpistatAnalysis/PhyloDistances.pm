package EpistatAnalysis::PhyloDistances;
#This package is a collection of methods for analysis of mutations in site pairs on a phylogenetic tree
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
# Symbols to autoexport (:DEFAULT tag)
@EXPORT = qw(calc_subst_distance mk_distance_analysis); 
@EXPORT_OK = qw(mk_distance_analysis);
use Bio::Phylo::IO;
use Class::Struct;
use SitePairUtils;

#information about detected pairs
#struct SubstPairInfo =>{
#	bg_branch => '$',
#	#bg_bases => '@',
#	tg_branch => '$',
#	#tg_bases => '@',
#	count => '$'
#};
#
#struct SitePairInfo => {
#	bg_site => '$',
#	tg_site => '$',
#	epistat => '$',
#	npair_subst => '$',
#	subst_info => '@'
#};

sub calc_subst_distance{
	my ($node1,$node2)=@_;
	my $d=$node2->get_branch_length();
	if($node1!=$node2){
		$d/=2;
		my $pnode=$node2->get_parent;
		while($pnode!=$node1){
			$d+=$pnode->get_branch_length();
			$pnode=$pnode->get_parent;
		};
		$d+=$pnode->get_branch_length()/2;
	}else{
		$d/=3;
	};
	return $d;
}

struct BranchDistanceTo =>{
	root => '$',
	closest_leaf => '$',
	#avg_leaf => '$',
	farest_leaf => '$',
	root_indicator => '$'
};

sub init_distances{
	my $di=shift;
	$di->root(0);
	$di->closest_leaf(0);
	$di->farest_leaf(0);
	$di->root_indicator(0);
}

sub mk_distance_analysis{
	my ($tree,$ra_site_pair_info,$rh_bg_subst_branches,$rh_tg_subst_branches)=@_;
	my @sp_distance_info;
	my @bg_free_mutations;
	my @tg_free_mutations;
	SitePairUtils::find_free_mutations($tree,$ra_site_pair_info,$rh_bg_subst_branches,$rh_tg_subst_branches,\@bg_free_mutations,\@tg_free_mutations);
	die "\nError size of free mutation lists!" unless @bg_free_mutations==@{$ra_site_pair_info}&&@bg_free_mutations==@tg_free_mutations;
	my $tree_length=$tree->calc_tree_length();
	#init nodes
	my %tree_dist_info;
	foreach my $node ($tree->get_nodes){
		$tree_dist_info{$node->get_name}=BranchDistanceTo->new();
		init_distances($tree_dist_info{$node->get_name});
	}
	$tree->visit_breadth_first(
		-in => sub {
			my $node=shift;
			if(!$node->is_root){
				my $pnode=$node->get_parent;
				my $d=$tree_dist_info{$pnode->get_name}->root;
				$d+=$node->get_branch_length;
				$tree_dist_info{$node->get_name}->root($d);
			}
		}
	);
	$tree->visit_depth_first(
		-in => sub {
			my $node=shift;
			if(!$node->is_root){
				my $di=$tree_dist_info{$node->get_name};
				my $br_length=$node->get_branch_length;
				my ($min_dist,$max_dist,$sum_dist)=($tree_length,0,0);
				foreach my $chnode(@{$node->get_children}){
					my $ch_di=$tree_dist_info{$chnode->get_name};
					$min_dist=$ch_di->closest_leaf+$br_length if($min_dist>$ch_di->closest_leaf+$br_length);
					$max_dist=$ch_di->farest_leaf+$br_length if($max_dist<$ch_di->farest_leaf+$br_length);
				}
				if(!$node->is_terminal){
					$di->closest_leaf($min_dist);
					$di->farest_leaf($max_dist);
				}
			}
		}
	);
	for(my $i=0;$i<@{$ra_site_pair_info};$i++){
		my $spi=$ra_site_pair_info->[$i];
		my $sp_di=[(0,BranchDistanceTo->new(),BranchDistanceTo->new(),BranchDistanceTo->new())];
		init_distances($sp_di->[1]);
		init_distances($sp_di->[2]);
		init_distances($sp_di->[3]);
		my $d=0;
		my $root_dist=0;
		my $root_ind=0;
		#my %bg_counts;
		foreach my $si(@{$spi->subst_info}){
			my $bg_name=$si->bg_branch->get_name;
			my $tg_name=$si->tg_branch->get_name;
			#$bg_counts{$bg_name}=1;
			my $u=$tree_dist_info{$tg_name}->root;
			my $v=$u;
			$d+=calc_subst_distance($si->bg_branch,$si->tg_branch);
			if(!$si->tg_branch->is_terminal){
				$v+=$tree_dist_info{$tg_name}->farest_leaf;
			}
			$root_dist+=$tree_dist_info{$bg_name}->root;
			$root_ind+=$u/$v;
		}
		$d/=@{$spi->subst_info};
		#$root_dist/=scalar keys(%bg_counts);
		$root_dist/=@{$spi->subst_info};
		$root_ind/=@{$spi->subst_info};
		$sp_di->[0]=$d;
		$sp_di->[1]->root($root_dist);
		$sp_di->[1]->root_indicator($root_ind);
		#free mutations
		$root_dist=0;
		$root_ind=0;
		foreach my $ni(@{$bg_free_mutations[$i]}){
			my $u=$tree_dist_info{$ni->get_name}->root;
			my $v=$u;
			if(!$ni->is_terminal){
				$v+=$tree_dist_info{$ni->get_name}->farest_leaf;
			};
			$root_dist+=$tree_dist_info{$ni->get_name}->root;
			$root_ind+=$u/$v;
		}
		if(@{$bg_free_mutations[$i]}){
			$root_dist/=@{$bg_free_mutations[$i]};
			$root_ind/=@{$bg_free_mutations[$i]};
		}else{
			$root_ind=undef;
		};
		$sp_di->[2]->root($root_dist);
		$sp_di->[2]->root_indicator($root_ind);
		$root_dist=0;
		$root_ind=0;
		foreach my $ni(@{$tg_free_mutations[$i]}){
			my $u=$tree_dist_info{$ni->get_name}->root;
			my $v=$u;
			if(!$ni->is_terminal){
				$v+=$tree_dist_info{$ni->get_name}->farest_leaf;
			};
			$root_dist+=$tree_dist_info{$ni->get_name}->root;
			$root_ind+=$u/$v;
		}
		if(@{$tg_free_mutations[$i]}){
			$root_dist/=@{$tg_free_mutations[$i]};
			$root_ind/=@{$tg_free_mutations[$i]};
		}else{
			$root_ind=undef;
		};
		$sp_di->[3]->root($root_dist);
		$sp_di->[3]->root_indicator($root_ind);
		push @sp_distance_info,$sp_di;
	}
	return @sp_distance_info;
}

1;
