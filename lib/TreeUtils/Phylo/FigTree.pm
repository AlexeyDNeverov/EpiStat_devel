#This package provides utilities for Bio::Phylo::Forest::Tree object conversion into FigTree nexus format
#The functions in the package assume that a tree object has names for internal nodes
use strict;
use Bio::Phylo::IO;

#This function prints a tree with node tags in FigTree nexus format
#Usage:
#tree2str($tree, dNdS => \%dNdS, reassortment => \%reassortments,...)
sub tree2str{
	my $tree=shift;
	my %args=@_;
	my %clade_names;
	foreach my $label(%args){
		foreach my $name(keys %{$args{$label}}){
			$clade_names{$name}=$label;
		};
	};
	my $str="";
	$tree->visit_depth_first(
        '-pre_daughter'   => sub { $str.='('             },     
        '-post_daughter'  => sub { $str.=')'			 },     
        '-in'             => sub { 
			my $node=shift;
			return if(!$node->get_parent);
			my $name=$node->get_name;
			if($node->is_terminal){
				$str.=$name;
			};
			$str.="[&" if (!$node->is_terminal)||defined($clade_names{$name});
			my $n=0;
			if(defined($args{color}) && defined($args{color}->{$name})){
				$str.="!color=#$args{color}->{$name}";
				$n++;
			};
			if(!$node->is_terminal){
				$str.=',' if $n;
				$str.="Name=\"$name\"";
				$n++;
			};
			foreach my $param (keys %args){
				next if $param eq 'color';
				if(defined $args{$param}->{$name}){
					$str.=',' if $n;
					$str.="$param=\"";
					#$str.="$param=";
					my $val;
					my $ref=ref($args{$param}->{$name});
					if($ref eq "ARRAY"){
						$val=$args{$param}->{$name}->[0];
					}elsif($ref eq "SCALAR"){
						$val=${$args{$param}->{$name}};
					}elsif($ref eq ""){
						$val=$args{$param}->{$name};
					}else{
						warn "\nWarning tree2str(): The type $ref is not supported as branch labels!";
					};
					$str.="$val\"" ;
					#$str.=$val;
					$n++;
				};
			};
			$str.="]" if (!$node->is_terminal)||defined($clade_names{$name});
			my $len=$node->get_branch_length;
			$str.=":$len";
		},
        '-pre_sister'     => sub { $str.=','             },     
	);
	$str.= ';';
}

sub pruned_tree2str{
	my $tree_=shift;
	my %args=@_;
	my %clade_names;
	foreach my $label(%args){
		foreach my $name(keys %{$args{$label}}){
			$clade_names{$name}=1;
		}
	}
	my $str = Bio::Phylo::IO->unparse(
              '-format' => 'adjacency',
			  '-phylo' => $tree_);
	my $tree=Bio::Phylo::IO->parse(
		'-format' => 'adjacency',
		'-string' =>$str);
	$tree=$tree->[0]->first;
	$tree->visit_depth_first(  
        '-post' => sub {
			my $node=shift;
			my $name=$node->get_name;
			unless($node->is_terminal()||defined($clade_names{$name})){
				my $n=0;
				foreach my $chnode(@{$node->get_children}){
					if(defined $clade_names{$chnode->get_name}){
						$n=1;
					}
				}
				if($n){
					$clade_names{$name}=1;
				}else{
					foreach my $chnode(@{$node->get_children}){
						$node->prune_child($chnode);
					}
				}
			}
		},
	);
	$tree->remove_orphans();
	$tree->remove_unbranched_internals();
	return tree2str($tree,@_);
}

1;