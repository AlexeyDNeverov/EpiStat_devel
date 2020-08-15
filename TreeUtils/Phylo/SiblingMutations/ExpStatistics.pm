package TreeUtils::Phylo::SiblingMutations::ExpStatistics;
#This module calculates exponential distance statistics for sibling mutations
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
@ISA = qw(TreeUtils::Phylo::SiblingMutations::BasicVisitor);

use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
@EXPORT = qw(); # Symbols to autoexport (:DEFAULT tag)
@EXPORT_OK = qw(SitePairExpStatistics);

use strict;
use Bio::Phylo::IO;
use Class::Struct;
use Bit::Vector;
use TreeUtils::Phylo::SiblingMutations::Find;

#matrix of phylogenetic distances between sites
my $mtx_init_value=0;
my $tau=50;

struct SitePairMatrix => {
	site_set1 => '$', #instance of Bit::Vector
	site_set2 => '$', #instance of Bit::Vector
	gene_index1 => '$',
	gene_index2 => '$', #gene_index2 == gene_index1 in the case of both sites within the same gene
	#Data matrix representations:
	rh_matrix_val => '$', #for sparse matrices
	ra_matrix_val => '$' #for dence matrices
};

struct SitePairExpStatistics => {
	bg_site => '$',
	tg_site => '$',
	value => '$'
};

#index access
sub _site_pair{
	my $self=shift;
	if(@_){
		my ($i,$j,$v)=@_;
		my $mtx=$self->{DATA};
		my $tmp;
		if(defined $mtx->rh_matrix_val){
			$tmp=$mtx->rh_matrix_val->{"$i,$j"};
			if(defined $tmp){
				$mtx->rh_matrix_val->{"$i,$j"}=$v if defined $v;
			};
		}elsif(defined $mtx->ra_matrix_val){
			$tmp=$mtx->ra_matrix_val->[$i]->[$j];
			$mtx->ra_matrix_val->[$i]->[$j]=$v if defined $v;
		};
		return $tmp;
	};
	return undef;
};

#interface declaration
#constructor
sub new{
	my $class=shift;
	my $self={};
	bless($self,$class);
	$self->_init(@_);
	return $self;
};

sub _init{
	my $self=shift;
	if(@_){
		my ($in_tree,$rh_mutmap,$rh_bgr_mutmap,$ra_site_pairs);
		my %args = @_;
		$self->{TAU}=$tau;
		$self->{F_TERMINALS}=1;
		while(my ($k,$v)=each %args){
			if($k eq "-tree"){
				$in_tree=$v;
			}elsif($k eq "-mutation_map"){
				$rh_mutmap=$v;
			}elsif($k eq "-bgr_mutation_map"){
				$rh_bgr_mutmap=$v;
			}elsif($k eq "-site_pairs"){
				$ra_site_pairs=$v;
			}elsif($k eq "-tau"){
				$self->{TAU}=$v;
			}elsif($k eq "-terminals"){
				$self->{F_TERMINALS}=0 if $v=~/^(false|FALSE|0)/;
			}else{
				die "\nUnknown parameter: $k!";
			};
			$self->{$k}=$v;
		};
		die "\nInput tree was not specified!" unless defined $in_tree;
		die "\nNot defined mutations on the input tree branches!" unless defined $rh_mutmap;
		#$self->{TREE}=$in_tree->clone();
		$self->{TREE}=$in_tree;
		$self->{TREE}->set_generic();
		my $is_intergen=0;
		$is_intergen=1 if(defined $rh_bgr_mutmap);
		my $ra_site2idx=[];
		$ra_site2idx->[0]={};
		$ra_site2idx->[1]={} if $is_intergen;
		$self->{IDX2SITE}=[];
		$self->{IDX2SITE}->[0]=[];
		$self->{IDX2SITE}->[1]=[] if $is_intergen;
		$self->{DATA}=SitePairMatrix->new();
		$self->{DATA}->gene_index1(0);
		if(defined $rh_bgr_mutmap){
			$self->{DATA}->gene_index2(1);
		}else{
			$self->{DATA}->gene_index2(0);
		};
		if(defined $ra_site_pairs){
			$self->{DATA}->rh_matrix_val({});
			my @si=(0,0); #site counters for two genes
			my %sites1;
			my %sites2;
			foreach my $spi(@{$ra_site_pairs}){
				if(!defined $ra_site2idx->[0]->{$spi->bg_site}){
					$ra_site2idx->[0]->{$spi->bg_site}=$si[0]++;
					push @{$self->{IDX2SITE}->[0]},$spi->bg_site;
				};
				if(!defined $ra_site2idx->[$is_intergen]->{$spi->tg_site}){
					$ra_site2idx->[$is_intergen]->{$spi->tg_site}=$si[$is_intergen]++;
					push @{$self->{IDX2SITE}->[$is_intergen]},$spi->tg_site;
				};
				my ($s1,$s2)=($ra_site2idx->[0]->{$spi->bg_site},$ra_site2idx->[$is_intergen]->{$spi->tg_site});
				$sites1{$s1}=1;
				$sites2{$s2}=1;
				$self->{DATA}->rh_matrix_val->{"$s1,$s2"}=$mtx_init_value;
			};
			my $n=@{$self->{IDX2SITE}->[0]};
			my @tmp=keys %sites1;
			$self->{DATA}->site_set1(Bit::Vector->new($n));
			foreach my $i(@tmp){
				$self->{DATA}->site_set1->Bit_On($i);
			};
			$n=@{$self->{IDX2SITE}->[1]} if $is_intergen;
			@tmp=keys %sites2;
			$self->{DATA}->site_set2(Bit::Vector->new($n));
			foreach my $i(@tmp){
				$self->{DATA}->site_set2->Bit_On($i);
			};
		}else{#no restriction for site pairs, use dence matrix data representation
			my $si=0;
			$self->{DATA}->ra_matrix_val=[];
			if(defined $rh_bgr_mutmap){
				foreach my $node($self->{TREE}->get_nodes){
					my $name=$node->get_name();
					foreach my $site(@{$rh_bgr_mutmap->{$name}}){
						if(!defined $ra_site2idx->[0]->{$site}){
							$ra_site2idx->[0]->{$site}=$si++;
							push @{$self->{IDX2SITE}->[0]},$site;
						};
					}
				};
				$si=0;
			}
			foreach my $node($self->{TREE}->get_nodes){
				my $name=$node->get_name();
				foreach my $site(@{$rh_mutmap->{$name}}){
					if(!defined $ra_site2idx->[$is_intergen]->{$site}){
						$ra_site2idx->[$is_intergen]->{$site}=$si++;
						push @{$self->{IDX2SITE}->[$is_intergen]},$site;
					};
				}
			};
			my $n=@{$self->{IDX2SITE}->[0]};
			$self->{DATA}->site_set1(Bit::Vector->new($n));
			$n=@{$self->{IDX2SITE}->[$is_intergen]};
			$self->{DATA}->site_set2(Bit::Vector->new($n));
			for(my $i=0;$i<@{$self->{IDX2SITE}->[0]};$i++){
				$self->{DATA}->ra_matrix_val->[$i]=[];
				$self->{DATA}->site_set1->Bit_On($i);
				for(my $j=0;$j<@{$self->{IDX2SITE}->[$is_intergen]};$j++){
					$self->{DATA}->ra_matrix_val->[$i]->[$j]=$mtx_init_value;
					$self->{DATA}->site_set2->Bit_On($j);
				}
			}
		};
		$self->{SITE2IDX}=$ra_site2idx;
		$self->_init_tree();
	};
	$self->_calculate();
};

struct MutationInfo =>{
	up_count => '$',
	down_dist => '@',
	down_ref => '@'
};

sub _init_tree{
	my $self=shift;
	foreach my $node($self->{TREE}->get_nodes){
		$node->set_generic("-sites_on_branch" => []);
		$node->get_generic("-sites_on_branch")->[0]=[];
		$node->get_generic("-sites_on_branch")->[1]=[] if $is_intergen;
		my $name=$node->get_name();
		if(defined $rh_bgr_mutmap){
			foreach my $site(@{$rh_bgr_mutmap->{$name}}){
				my $sidx=$ra_site2idx->[0]->{$site};
				push @{$node->get_generic("-sites_on_branch")->[0]},$sidx if(defined $sidx);
			};
		};
		foreach my $site(@{$rh_mutmap->{$name}}){
			my $sidx=$ra_site2idx->[$is_intergen]->{$site};
			push @{$node->get_generic("-sites_on_branch")->[$is_intergen]},$sidx if(defined $sidx);
		};
		$node->set_generic("-sites_access" => []);
		$node->get_generic("-sites_access")->[0]=[];
		$node->get_generic("-sites_access")->[1]=[] if $is_intergen;
		for(my $i=0;$i<@{$self->{IDX2SITE}->[0]};$i++){
			my $mi=MutationInfo->new();
			$mi->up_count(0);
			$node->get_generic("-sites_access")->[0]->[$i]=$mi;
		};
		if($is_intergen){
			for(my $i=0;$i<@{$self->{IDX2SITE}->[1]};$i++){
				my $mi=MutationInfo->new();
				$mi->up_count(0);
				$node->get_generic("-sites_access")->[1]->[$i]=$mi;
			}
		};
	};
	#calculate upward site occurrences
	$self->{TREE}->visit_breadth_first(
		-in => sub {
			my $node=shift;
			return if $node->is_root;
			my $pnode=$node->get_parent;
			my $ra_mut_info=$node->get_generic("-sites_access");
			my $ra_mut_info_pnode=$pnode->get_generic("-sites_access");
			my $ra_mut_on_branche=$node->get_generic("-sites_on_branch");
			my $ngenes=@{$ra_mut_info};
			my $dist=$node->get_branch_length;
			$dist+=$pnode->get_branch_length unless $pnode->is_root;
			$dist/=2;
			for(my $I=0;$I<$ngenes;$I++){
				for(my $i=0;$i<@{$ra_mut_info->[$I]};$i++){
					$ra_mut_info->[$I]->[$i]->up_count($ra_mut_info_pnode->[$I]->[$i]->up_count);
				};
				foreach my $sidx(@{$ra_mut_on_branche->[$I]}){
					$ra_mut_info->[$I]->[$sidx]->up_count($ra_mut_info->[$I]->[$sidx]->up_count+1);
				};
			};
		}
	);
	#find closest downward mutations on independent phylogenetic lineages
	$self->{TREE}->visit_depth_first(
		-in => sub {
			my $node=shift;
			return if $node->is_root;
			my $pnode=$node->get_parent;
			my $ra_mut_info=$node->get_generic("-sites_access");
			my $ra_mut_info_pnode=$pnode->get_generic("-sites_access");
			my $ra_mut_on_branche=$node->get_generic("-sites_on_branch");
			my $ngenes=@{$ra_mut_info};
			my $dist=$node->get_branch_length;
			$dist+=$pnode->get_branch_length unless $pnode->is_root;
			$dist/=2;
			for(my $I=0;$I<$ngenes;$I++){
				my %sites_on_branch;
				foreach my $sidx(@{$ra_mut_on_branche->[$I]}){
					$sites_on_branch{$sidx}=1;
					push @{$ra_mut_info_pnode->[$I]->[$sidx]->down_ref},$node;
					push @{$ra_mut_info_pnode->[$I]->[$sidx]->down_dist},$dist;
				};
				for(my $i=0;$i<@{$ra_mut_info->[$I]};$i++){
					if(!exists $sites_on_branch{$i}){
						for(my $j=0;$j<@{$ra_mut_info->[$I]->[$i]->down_ref};$j++){
							my $d=$dist+$ra_mut_info->[$I]->[$i]->down_dist($j);
							my $r=$ra_mut_info->[$I]->[$i]->down_ref($j);
							push @{$ra_mut_info_pnode->[$I]->[$i]->down_ref},$r;
							push @{$ra_mut_info_pnode->[$I]->[$i]->down_dist},$d;
						}
					}
				}
			}	
		}
	);
};

sub _is_accessible_site_pair{
	my ($src_node,$mca_node,$node,$s_from,$s_to)=@_;
	my $is_intergenic=@{$src_node->get_generic("-sites_access")}-1;
	my $s_ok=1;
	my $pnode=$src_node->get_parent;
	if($pnode!=$mca_node){
		$s_ok=0 if $mca_node->get_generic("-sites_access")->[0]->[$s_from]->up_count-$pnode->get_generic("-sites_access")->[0]->[$s_from]->up_count;
		$s_ok=0 if $mca_node->get_generic("-sites_access")->[$is_intergenic]->[$s_to]->up_count-$pnode->get_generic("-sites_access")->[$is_intergenic]->[$s_to]->up_count;
	};
	$pnode=$node->get_parent;
	if($pnode!=$mca_node){
		$s_ok=0 if $mca_node->get_generic("-sites_access")->[0]->[$s_from]->up_count-$pnode->get_generic("-sites_access")->[0]->[$s_from]->up_count;
		$s_ok=0 if $mca_node->get_generic("-sites_access")->[$is_intergenic]->[$s_to]->up_count-$pnode->get_generic("-sites_access")->[$is_intergenic]->[$s_to]->up_count;
	};
	return $s_ok;
};

sub gen_iterator{
	my $self=shift;
	my ($s1_idx,$s2_idx);
	my ($s1_last,$s2_last);
	my $is_intergen=@{$self->{IDX2SITE}}-1;
	my @mtx_keys;
	if(defined $self->{DATA}->rh_matrix_val){
		@mtx_keys=keys %{$self->{DATA}->rh_matrix_val};
		$s1_idx=0;
		$s1_last=@mtx_keys;
	}elsif(defined $self->{DATA}->ra_matrix_val){
		$s1_idx=0;
		$s1_last=@{$self->{IDX2SITE}->[0]};
		$s2_idx=0;
		$s2_last=@{$self->{IDX2SITE}->[$is_intergen]};
	};
	return sub {
		my $v=SitePairExpStatistics->new();
		if(defined $self->{DATA}->rh_matrix_val){
			return undef if $s1_idx == $s1_last;
			my @idx=split ',', $mtx_keys[$s1_idx];
			$v->bg_site($self->{IDX2SITE}->[0]->[$idx[0]]);
			$v->tg_site($self->{IDX2SITE}->[$is_intergen]->[$idx[1]]);
			$v->value($self->{DATA}->rh_matrix_val->{$mtx_keys[$s1_idx]});
			$s1_idx++;
			return $v; 
		}elsif(defined $self->{DATA}->ra_matrix_val){
			return undef if $s1_idx == $s1_last || $s2_idx == $s2_last;
			$v->bg_site($self->{DATA}->{IDX2SITE}->[0]->[$s1_idx]);
			$v->tg_site($self->{DATA}->{IDX2SITE}->[$is_intergen]->[$s2_idx]);
			$v->value($self->{DATA}->ra_matrix_val->[$s1_idx]->[$s2_idx]);
			$s2_idx++;
			if($s2_idx==$s2_last){
				$s1_idx++;
				$s2_idx=0;
			};
			return $v;
		};
		return undef;
	};
};
#Private

sub _calculate{
	my $self=shift;
	my $n=$self->{TREE}->get_nodes();
	my $n0=$n;
	foreach my $node($self->{TREE}->get_nodes()){
		next if $node->is_root;
		findSiblings($node,$self);
		$n--;
#last if $n0-$n > 10;
		print STDERR "\nLeft to process $n nodes!";
	}
};

sub update_node{
	my $self=shift;
	my ($node_source,$node_branch,$node_target)=@_;
	my $l=$node_source->get_score()+$node_branch->get_branch_length();
	$node_target->set_score($l);
	return 1;
};

sub update_visitor{
	my $self=shift;
	my $site_pair_mtx=$self->{DATA};
	my ($src_node,$mca_node,$node)=@_;
	return 0 if $node->is_terminal && !$self->{F_TERMINALS};
	my $pnode=$node->get_parent;
	my $l=$mca_node->get_score()+$src_node->get_branch_length()/2;
###HERE
	my %both_sites1;
	my %both_sites2;
	my $n=0;
	foreach my $s2_idx(@{$node->get_generic("-sites_on_branch")->[$site_pair_mtx->gene_index2]}){
		next unless $sites2_bvec->bit_test($s2_idx);
		my $p=1.0;
		$p=0.5 if $both_sites2->bit_test($s2_idx);
		foreach my $s1_idx(@{$src_node->get_generic("-sites_on_branch")->[$site_pair_mtx->gene_index1]}){
			next unless $sites1_bvec->bit_test($s1_idx);
			$p*=0.5 if $both_sites1->bit_test($s1_idx);
			my $v=$self->_site_pair($s1_idx,$s2_idx);
			if(defined $v){
				$v+=$p*exp(-$l/$tau);
				$self->_site_pair($s1_idx,$s2_idx,$v);
				$n++;
			};
		}
	}
	return $n;
};

sub init_search{
	my $self=shift;
	my $node=shift;
	return if $node->is_root;
	my $dist_mtx=$self->{DATA};
	my $pnode=$node->get_parent;
	$pnode->set_generic("-sites_access" => []);
	if($dist_mtx->gene_index1==$dist_mtx->gene_index2){
		my $bvec=$dist_mtx->site_set1->Shadow();
		$bvec->Or($dist_mtx->site_set1,$dist_mtx->site_set2);
		$pnode->get_generic("-sites_access")->[0]=$bvec;
	}else{
		$pnode->get_generic("-sites_access")->[0]=$dist_mtx->site_set1->Clone();
		$pnode->get_generic("-sites_access")->[1]=$dist_mtx->site_set2->Clone();
	};
	$pnode->set_score(0);
};

sub clean_node{
	my $self=shift;
	my $node=shift;
	$node->set_generic("-sites_access" => undef);
	$node->set_score(0);
};

1;