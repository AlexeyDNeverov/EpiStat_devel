package Pairs::ConsecutivePairs;
use Bio::Phylo::Forest::Tree;
use SitePairInfo;
use EpistatAnalysis::PhyloDistances qw(calc_subst_distance);
use strict;
#constants
use constant {
	ALL_PAIRS =>0,
	INDEP_PAIRS =>1,
};

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
	my ($f_intragene,$tree,$f_ignore_terminals,$hr_bg_subst_map,$hr_tg_subst_map,$hr_bg_site_idx,$hr_tg_site_idx,$ra_bg_site2cord,$f_branch_pairs)=@_;
	$f_branch_pairs=1 unless defined $f_branch_pairs;
	$self->{SITE_PAIRS}={};
	$self->{BGR_NSITES}=keys %{$hr_bg_site_idx};
	$self->{FGR_NSITES}=keys %{$hr_tg_site_idx};
	my %bg_subst_upref;
	my %tg_subst_upref;
	$tree->visit_breadth_first(
		-pre => sub {
			my $node=shift;
			$bg_subst_upref{$node->get_name}=[];
			$tg_subst_upref{$node->get_name}=[];
			$node->set_generic("level" => 0);
			if(!($node->is_root||($f_ignore_terminals&&$node->is_terminal))){
				my $pnode=$node->get_parent();
				@{$bg_subst_upref{$node->get_name}}=@{$bg_subst_upref{$pnode->get_name}};
				@{$tg_subst_upref{$node->get_name}}=@{$tg_subst_upref{$pnode->get_name}};
			}
		},
		-in   => sub {
			my $node=shift;
			if(!($node->is_root||($f_ignore_terminals&&$node->is_terminal))){
				$node->set_generic("level" => $node->get_parent()->get_generic("level")+1);		
				my %bg_subst;
				if(defined $hr_bg_subst_map->{$node->get_name}){
					foreach my $i(keys %{$hr_bg_subst_map->{$node->get_name}}){
						$bg_subst{$i}=1;
					}
				};
				if(defined $hr_tg_subst_map->{$node->get_name}){
					#generate possible mutation pairs
					foreach my $tg_spos(keys %{$hr_tg_subst_map->{$node->get_name}}){
						my $tg_site_idx=$hr_tg_site_idx->{$tg_spos};
						my $ptg_node=$tg_subst_upref{$node->get_name}->[$tg_site_idx];
						for(my $i=0;$i<@{$ra_bg_site2cord};$i++){
							my $bg_spos=$ra_bg_site2cord->[$i];
							next if($f_intragene&&$bg_spos==$tg_spos);
							my $pbg_node=$bg_subst_upref{$node->get_name}->[$i];
							next unless defined($pbg_node)||defined($bg_subst{$bg_spos}); #no substitution in background for the site $i
							if(defined($ptg_node)&&defined($pbg_node)){
								#There is a substitution in the same site in the target followed the substitution in the background 
								next unless ($ptg_node->get_generic("level")<=$pbg_node->get_generic("level"))||defined($bg_subst{$bg_spos});
							};
							my $p_idx=$self->_idx_pair2line($i,$tg_site_idx);
							my $spi=$self->{SITE_PAIRS}->{$p_idx};
							if(!defined $spi){
								$spi=SitePairInfo->new();
								$spi->bg_site($bg_spos);
								$spi->tg_site($tg_spos);
							}
	#my $bg_name=$pbg_node->get_name;
	#my $tg_name=$node->get_name;
	#print "\n$bg_spos\t$tg_spos\t$bg_name\t$tg_name\t$i\t$tg_site_idx";
							#the site pair has calculated epistatic statistics
							my $p=1;
							if(defined $bg_subst{$bg_spos}){
								$p=0.5;
								SitePairInfo::add_subst_info($spi,$node,$node,$p) if $f_branch_pairs;
							};
							if(defined $pbg_node){
								if((!defined($ptg_node)) || ($ptg_node->get_generic("level")<=$pbg_node->get_generic("level"))){
									$p*=0.5 if($pbg_node==$ptg_node);
									SitePairInfo::add_subst_info($spi,$pbg_node,$node,$p);
								}
							};
							unless(defined $self->{SITE_PAIRS}->{$p_idx}){
								$self->{SITE_PAIRS}->{$p_idx}=$spi if @{$spi->subst_info};
							}
						};
						#update the target site references
						$tg_subst_upref{$node->get_name}->[$tg_site_idx]=$node;
					}
				};
				if(defined $hr_bg_subst_map->{$node->get_name}){
					#update the background site references
					foreach my $bg_spos(keys %{$hr_bg_subst_map->{$node->get_name}}){
						my $bg_idx=$hr_bg_site_idx->{$bg_spos};
						$bg_subst_upref{$node->get_name}->[$bg_idx]=$node;
					}
				};
			};	
		},
		-post => sub{
			my $node=shift;
			if(!($node->is_root||$node->is_terminal)){
				delete $tg_subst_upref{$node->get_parent->get_name};
				delete $bg_subst_upref{$node->get_parent->get_name};
			};
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
	my ($tau,$f_stat_type,$hr_bg_subst_map,$hr_tg_subst_map)=@_;
	die "\nConsequtivePairs::calc_statistics() Error: Undefined scale parameter tau!" unless defined $tau;
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
			my $t=calc_subst_distance($bg_node,$tg_node);
			my $count=$si->count();
			my $esi=$count*exp(-$t/$tau);
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
			if($f_stat_type==ALL_PAIRS){
				$epistat+=$esi;
			}elsif($f_stat_type==INDEP_PAIRS){
				my $name=$bg_node->get_name;
				$bgr_epistat{$name}=[(0,0)] unless defined $bgr_epistat{$name};
				$bgr_epistat{$name}->[0]+=$esi;
				$bgr_epistat{$name}->[1]+=1;
			}else{
				die "\nError calc_epistat(): Unknown type of statistics!";
			}
		}
		if($f_stat_type==INDEP_PAIRS){
			foreach my $name(keys %bgr_epistat){
				$epistat+=$bgr_epistat{$name}->[0]/$bgr_epistat{$name}->[1];
			}
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