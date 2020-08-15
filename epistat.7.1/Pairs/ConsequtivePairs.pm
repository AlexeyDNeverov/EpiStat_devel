package Pairs::ConsequtivePairs;
use Bio::Phylo::Forest::Tree;
use SitePairInfo;
use Pairs::Stat qw(stat_func);
use strict;

sub find{
	my ($f_intragene,$tree,$f_ignore_terminals,$hr_bg_subst_map,$hr_tg_subst_map,$hr_bg_site_idx,$hr_tg_site_idx,$ra_bg_site2cord)=@_;
	my %bg_subst_upref;
	my %tg_subst_upref;
	my %site_pairs;
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
							my $str="$i,$tg_site_idx";
							if(!defined $site_pairs{$str}){
								my $spi=SitePairInfo->new();
								$spi->bg_site($bg_spos);
								$spi->tg_site($tg_spos);
								$site_pairs{$str}=$spi;
							};
	#my $bg_name=$pbg_node->get_name;
	#my $tg_name=$node->get_name;
	#print "\n$bg_spos\t$tg_spos\t$bg_name\t$tg_name\t$i\t$tg_site_idx";
							#the site pair has calculated epistatic statistics
							my $p=1;
							if(defined $bg_subst{$bg_spos}){
								$p=0.5;
								SitePairInfo::add_subst_info($site_pairs{$str},$node,$node,$p);
							};
							if(defined $pbg_node){
								if((!defined($ptg_node)) || ($ptg_node->get_generic("level")<=$pbg_node->get_generic("level"))){
									$p*=0.5 if($pbg_node==$ptg_node);
									SitePairInfo::add_subst_info($site_pairs{$str},$pbg_node,$node,$p);
								}
							};
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
	return %site_pairs;
}

sub calc_statistics{
	my ($site_pair_info,$f_stat_type,$f_stat_func,$tau,$hr_bg_subst_map,$hr_tg_subst_map)=@_;
	die "\nConsequtivePairs::calc_statistics() Error: Undefined scale parameter tau!" unless defined($tau)||($f_stat_func==0);
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
		my $esi;
		if($f_stat_func==1){
			$esi=stat_func(sub {exp(-$_[0]/$_[1]);},$si,$tau,$w1,$w2);
		}else{
			$esi=stat_func(sub {$_[0];},$si,$tau,$w1,$w2);
		}
		if($f_stat_type==0){
			$epistat+=$esi;
		}elsif($f_stat_type==1){
			my $name=$bg_node->get_name;
			$bgr_epistat{$name}=[(0,0)] unless defined $bgr_epistat{$name};
			$bgr_epistat{$name}->[0]+=$esi;
			$bgr_epistat{$name}->[1]+=1;
		}elsif($f_stat_type==3){
			$epistat+=$esi unless $bg_node==$tg_node;
		}else{
			die "\nError calc_epistat(): Unknown type of statistics!";
		}
	}
	if($f_stat_type==1){
		foreach my $name(keys %bgr_epistat){
			$epistat+=$bgr_epistat{$name}->[0]/$bgr_epistat{$name}->[1];
		}
	}
	$site_pair_info->epistat($epistat);
	return $epistat;
}

1;