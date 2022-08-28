package EpistatAnalysis::PhyloDrawSites;
#This package is a collection of methods for analysis of mutations in site pairs on a phylogenetic tree
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
# Symbols to autoexport (:DEFAULT tag)
@EXPORT = qw(site_pairs2figtree_str sites2figtree_str subst_distrib2figtree_str); 
@EXPORT_OK = qw(site_pairs2figtree_str sites2figtree_str subst_distrib2figtree_str);
use Bio::Phylo::IO;
use Class::Struct;
use List::BinarySearch qw(binsearch binsearch_pos);
use Math::Gradient qw(multi_gradient multi_array_gradient);
use Color::Rgb;
use SitePairUtils;
use TreeUtils::Phylo::FigTree qw(tree2str prunned_tree2str);
die "Set up the system variable \$EPISTAT_HOME" unless defined $ENV{EPISTAT_HOME};
my $rgb_txt=$ENV{EPISTAT_HOME}."/rgb.txt";
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

my $bgr_color="00ffff";
my $fgr_color="ff00ff";
my $free_color="ff0033";
my $mix_color="00cc00";

sub unparse_subst_abbr{
	my ($rh_subst_map,$node_name,$spos,$ra_idx2allele)=@_;
	my $si=$rh_subst_map->{$node_name}->{$spos};
	my @bases=@{$si->bases};
	if(defined $ra_idx2allele){
		$bases[0]=$ra_idx2allele->[$bases[0]];
		$bases[1]=$ra_idx2allele->[$bases[1]];
	}
	my $str=$bases[0].$spos.$bases[1];
	my $w=$si->weight;
	if(defined $w){
		$w=sprintf("%.3f",$w) if($w=~/\d+\.\d+/);
		$str.=":$w";
	}
	return $str;
}

sub _mk_phenotypes_info{
	my ($tree,$rh_phenotypes,$ra_pheno_labels,$npheno)=@_;
	my %phenotypes_info;
	foreach my $node($tree->get_nodes){
		next if $node->is_root;
		my $name=$node->get_name;
		my $pname=$node->get_parent->get_name;
		my $str;
		for(my $i=0;$i<$npheno;$i++){
			unless($rh_phenotypes->{$name}->[$i]==$rh_phenotypes->{$pname}->[$i]){
				my $tmp=$rh_phenotypes->{$pname}->[$i]."_".$ra_pheno_labels->[$i]."_".$rh_phenotypes->{$name}->[$i];
				unless(defined $str){
					$str=$tmp;
				}else{
					$str.=";$tmp";
				}
			}
		}
		$phenotypes_info{$name}=$str if defined $str;
	}
	return %phenotypes_info;
}

#'ra_fig_site_pair_info' - a ref to an array of SitePairInfo objects that are required to be drawn
sub site_pairs2figtree_str{
	my ($ra_fig_site_pair_info,$rh_bg_subst_branches,$rh_tg_subst_branches,$rh_bg_subst_map,$rh_tg_subst_map,$tree,%args)=@_;
	my $ra_idx2allele=$args{'-ra_idx2allele'};
	my $rh_branch_confidences=$args{'-rh_branch_confidences'};
	my $f_prune_tree=0;
	$f_prune_tree=$args{'-f_prune_tree'} if(defined $args{'-f_prune_tree'});
	my $npheno=0;
	my $ra_pheno_labels=$args{'-ra_pheno_labels'};
	$npheno=@{$ra_pheno_labels} if defined $ra_pheno_labels;
	my %phen_change_info;
	my %phenotypes_info;
	my $rh_phenotypes=$args{'-rh_phenotypes'};
	if($npheno){
		$rh_phenotypes=$args{'-rh_phenotypes'};
		die "\nError site_pairs2figtree_str(): The reference for branch phenotypes is expected!" unless defined $rh_phenotypes;
		%phen_change_info=_mk_phenotypes_info($tree,$rh_phenotypes,$ra_pheno_labels,$npheno);
		foreach my $name(keys %{$rh_phenotypes}){
			$phenotypes_info{$name}=join('',@{$rh_phenotypes->{$name}});
		}
	}
	my @fig_site_pairs=sort {$a->bg_site<=>$b->bg_site||$a->tg_site<=>$b->tg_site} @{$ra_fig_site_pair_info};
	my @bg_free_mutations;
	my @tg_free_mutations;
	SitePairUtils::find_free_mutations($tree,\@fig_site_pairs,$rh_bg_subst_branches,$rh_tg_subst_branches,\@bg_free_mutations,\@tg_free_mutations);
	my %mutations_info;
	my %branch_color;
	for(my $i=0;$i<@fig_site_pairs;$i++){
		my %bg_names;
		my %tg_names;
		foreach my $spi(@{$fig_site_pairs[$i]->subst_info}){
			$bg_names{$spi->bg_branch->get_name}=1;
			$tg_names{$spi->tg_branch->get_name}=1;
		}
		my $bg_site=$fig_site_pairs[$i]->bg_site;
		my $tg_site=$fig_site_pairs[$i]->tg_site;
		foreach my $name(keys %bg_names){
			$mutations_info{$name}.=";" if defined $mutations_info{$name};
			my $str=unparse_subst_abbr($rh_bg_subst_map,$name,$bg_site,$ra_idx2allele);
			$mutations_info{$name}.="BGR_".$str;
			if(!defined $branch_color{$name}){
				$branch_color{$name}=$bgr_color;
			}else{
				$branch_color{$name}=$mix_color;
			}
		};
		foreach my $name(keys %tg_names){
			$mutations_info{$name}.=";" if defined $mutations_info{$name};
			my $str=unparse_subst_abbr($rh_tg_subst_map,$name,$tg_site,$ra_idx2allele);
			$mutations_info{$name}.="FGR_".$str;
			if(!defined $branch_color{$name}){
				$branch_color{$name}=$fgr_color;
			}else{
				$branch_color{$name}=$mix_color;
			}
		};
		my $key=SitePairInfo->new();
		$key->bg_site($tg_site);
		$key->tg_site($bg_site);
		my $j=binsearch {$a->bg_site<=>$b->bg_site||$a->tg_site<=>$b->tg_site} $key, @fig_site_pairs;
		if((!defined($j))||$i<$j){
			foreach my $node(@{$bg_free_mutations[$i]}){
				my $name=$node->get_name;
				$mutations_info{$name}.=";" if defined $mutations_info{$name};
				my $str=unparse_subst_abbr($rh_bg_subst_map,$name,$bg_site,$ra_idx2allele);
				$mutations_info{$name}.="FREE_".$str;
				if(!defined $branch_color{$name}){
					$branch_color{$name}=$free_color;
				}else{
					$branch_color{$name}=$mix_color;
				}
			};
			foreach my $node(@{$tg_free_mutations[$i]}){
				my $name=$node->get_name;
				$mutations_info{$name}.=";" if defined $mutations_info{$name};
				my $str=unparse_subst_abbr($rh_tg_subst_map,$name,$tg_site,$ra_idx2allele);
				$mutations_info{$name}.="FREE_".$str;
				if(!defined $branch_color{$name}){
					$branch_color{$name}=$free_color;
				}else{
					$branch_color{$name}=$mix_color;
				}
			}
		}
	}
	my $str;
	my $rf_tree2str=\&tree2str;
	$rf_tree2str=\&pruned_tree2str if($f_prune_tree);
	my %inargs=(mutations => \%mutations_info, color => \%branch_color);
	$inargs{'confidence'}=$rh_branch_confidences if(defined $rh_branch_confidences);
	if($npheno){
		$inargs{'phen_change'}=\%phen_change_info;
		$inargs{'phen'}=\%phenotypes_info;
	}
	$str=&$rf_tree2str($tree,%inargs);
	return $str;
}

#'ra_fig_site_pair_pos' - a ref to an array of array references contained a pair of site positions "[($bg_site_pos,$fg_site_pos)]" that are required to be drawn
sub sites2figtree_str{
	my ($ra_fig_site_pair_pos,$rh_bg_subst_branches,$rh_tg_subst_branches,$rh_bg_subst_map,$rh_tg_subst_map,$tree,%args)=@_;
	my $ra_idx2allele=$args{'-ra_idx2allele'};
	my $rh_branch_confidences=$args{'-rh_branch_confidences'};
	my $f_prune_tree=0;
	$f_prune_tree=$args{'-f_prune_tree'} if(defined $args{'-f_prune_tree'});
	my @bg_site_pos;
	my @tg_site_pos;
	foreach my $spi(@{$ra_fig_site_pair_pos}){
		push @bg_site_pos,$spi->[0] if defined $spi->[0];
		push @tg_site_pos,$spi->[1] if defined $spi->[1];
	}
	my %mutations_info;
	my %branch_color;
	foreach my $spos(@bg_site_pos){
		foreach my $name(@{$rh_bg_subst_branches->{$spos}}){
			$mutations_info{$name}.=";" if defined $mutations_info{$name};
			my $str=unparse_subst_abbr($rh_bg_subst_map,$name,$spos,$ra_idx2allele);
			$mutations_info{$name}.="BGR_".$str;
			if(!defined $branch_color{$name}){
				$branch_color{$name}=$bgr_color;
			}else{
				$branch_color{$name}=$mix_color;
			}
		}
	}
	foreach my $spos(@tg_site_pos){
		foreach my $name(@{$rh_tg_subst_branches->{$spos}}){
			$mutations_info{$name}.=";" if defined $mutations_info{$name};
			my $str=unparse_subst_abbr($rh_tg_subst_map,$name,$spos,$ra_idx2allele);
			$mutations_info{$name}.="FGR_".$str;
			if(!defined $branch_color{$name}){
				$branch_color{$name}=$fgr_color;
			}else{
				$branch_color{$name}=$mix_color;
			}
		}
	}
	my $str;
	my $rf_tree2str=\&tree2str;
	$rf_tree2str=\&pruned_tree2str if($f_prune_tree);
	my %inargs=(mutations => \%mutations_info, color => \%branch_color);
	$inargs{'confidence'}=$rh_branch_confidences if(defined $rh_branch_confidences);
	$str=&$rf_tree2str($tree,%inargs);
	return $str;
}

sub _diff_mut_counts{
	my ($obs,$exp)=@_;
	return $obs - $exp;
}

sub subst_distrib2figtree_str{
	my ($ra_sites,$rh_subst_map,$tree)=@_;
	#calculating number of substitutions in subtrees of branches
	my %subst_counts;
	$tree->visit_depth_first(
		-pre => sub {
				my $node=shift;
				my $name=$node->get_name;
				die "\nUnnamed nodes in the tree are not allowed!" unless defined $name;
				$subst_counts{$name}=[(0) x 2];
			},
		-in => sub{
				my $node=shift;
				my $name=$node->get_name;
				foreach my $chnode(@{$node->get_children}){
					my $chname=$chnode->get_name;
					$subst_counts{$name}->[0]+=$subst_counts{$chname}->[0];
					$subst_counts{$name}->[1]+=$subst_counts{$chname}->[1]+$chnode->get_branch_length;
					if(defined $rh_subst_map->{$chname}){
						foreach my $site(@{$ra_sites}){
							if(defined $rh_subst_map->{$chname}->{$site}){
								my $w=$rh_subst_map->{$chname}->{$site}->weight;
								$w=1 unless defined $w;
								$subst_counts{$name}->[0]+=$w;
							}
						}
					}
				}
			}
	);
	my $rname=$tree->get_root->get_name;
	return undef unless($subst_counts{$rname}->[0]>0);
	my $lambda=$subst_counts{$rname}->[0]/$subst_counts{$rname}->[1];
	my @stat_range=($subst_counts{$rname}->[0],-$subst_counts{$rname}->[0]);
	for my $node(@{$tree->get_internals}){
		my $name=$node->get_name;
		my $str=_diff_mut_counts($subst_counts{$name}->[0],$subst_counts{$name}->[1]*$lambda);
		$stat_range[0]=$str if $stat_range[0]>$str;
		$stat_range[1]=$str if $stat_range[1]<$str;
	}
	my $ncat=10;
	my @stat_cats=multi_gradient($ncat,$stat_range[0],0,$stat_range[1]);
	#my(@hot_spots) = ([ 0, 255, 0 ], [ 255, 255, 0 ], [ 127, 127, 127 ], [ 0, 0, 255 ], [ 127, 0, 0 ], [ 255, 255, 255 ]);
	my(@hot_spots) = ([ 0, 0, 255 ], [ 255, 0, 0 ]);
	my @color_cats=multi_array_gradient($ncat, @hot_spots);
	my %branch_color;
	my $rgb = new Color::Rgb(rgb_txt=>$rgb_txt);
	my %obs_mut_counts;
	my %exp_nut_counts;
	for my $node(@{$tree->get_internals}){
		my $name=$node->get_name;
		my $str=_diff_mut_counts($subst_counts{$name}->[0],$subst_counts{$name}->[1]*$lambda);
		my $idx=binsearch_pos {$a<=>$b} $str, @stat_cats;
		if($idx<@stat_cats){
			$branch_color{$name}=$rgb->rgb2hex(@{$color_cats[$idx]});
		}else{
			$branch_color{$name}=$rgb->rgb2hex(@{$hot_spots[-1]});
		}
		$obs_mut_counts{$name}=$subst_counts{$name}->[0];
		$exp_nut_counts{$name}=sprintf("%.1f",$subst_counts{$name}->[1]*$lambda);
	}
	return tree2str($tree,obs_nmut => \%obs_mut_counts,exp_nmut => \%exp_nut_counts, color => \%branch_color);
}

1;
