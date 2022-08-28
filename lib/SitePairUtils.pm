package SitePairUtils;
#This package is a collection of methods for analysis of mutations in site pairs on a phylogenetic tree
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
# Symbols to autoexport (:DEFAULT tag)
@EXPORT = qw(find_free_mutations); 
@EXPORT_OK = qw(find_free_mutations);
use Bio::Phylo::IO;
use Class::Struct;
use List::BinarySearch qw(binsearch);

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

#Method searches for each site pair mutations in background and foreground sites which are not formed consequtitive pairs
#Args:
#$tree - an Bio::Phylo::Tree object
#$ra_site_pair_info - a reference on array of SitePairInfo structure objects
#$rh_bg_subst_branches - a reference on a hash indexed by site positions which stores a reference on array of names of nodes in the $tree containing mutations in a corresponding BACKGROUND site
#$rh_bg_subst_branches - a reference on a hash indexed by site positions which stores a reference on array of names of nodes in the $tree containing mutations in a corresponding FOREGROUND site
#Out args:
#$ra_bg_free_mut
#$ra_tg_free_mut
sub find_free_mutations{
	my ($tree,$ra_site_pair_info,$rh_bg_subst_branches,$rh_tg_subst_branches,
		#out:
		$ra_bg_free_mut,$ra_tg_free_mut)=@_;
	my @done;
	@{$ra_bg_free_mut}=();
	@{$ra_tg_free_mut}=();
	my %name2branch;
	foreach my $node($tree->get_nodes){
		my $name=$node->get_name;
		$name2branch{$name}=$node;
	}
	for(my $i=0;$i<@{$ra_site_pair_info};$i++){
		next if $done[$i];
		$ra_bg_free_mut->[$i]=[];
		$ra_tg_free_mut->[$i]=[];
		my $bg_site=$ra_site_pair_info->[$i]->bg_site;
		my $tg_site=$ra_site_pair_info->[$i]->tg_site;
		my %exclude_bg_branches;
		my %exclude_tg_branches;
		foreach my $spi(@{$ra_site_pair_info->[$i]->subst_info}){
			my ($bg_name,$tg_name)=($spi->bg_branch->get_name,$spi->tg_branch->get_name);
			$exclude_bg_branches{$bg_name}=1;
			$exclude_tg_branches{$tg_name}=1;
		}
		my $key=SitePairInfo->new();
		$key->bg_site($tg_site);
		$key->tg_site($bg_site);
		my $j=binsearch {$a->bg_site<=>$b->bg_site||$a->tg_site<=>$b->tg_site} $key, @{$ra_site_pair_info};
		if(defined $j){
			foreach my $spi(@{$ra_site_pair_info->[$j]->subst_info}){
				my ($bg_name,$tg_name)=($spi->bg_branch->get_name,$spi->tg_branch->get_name);
				$exclude_tg_branches{$bg_name}=1;
				$exclude_bg_branches{$tg_name}=1;
			}
		}
		foreach my $name(@{$rh_bg_subst_branches->{$bg_site}}){
			my $node=$name2branch{$name};
			push @{$ra_bg_free_mut->[$i]},$node unless defined $exclude_bg_branches{$name};
		}
		foreach my $name(@{$rh_tg_subst_branches->{$tg_site}}){
			my $node=$name2branch{$name};
			push @{$ra_tg_free_mut->[$i]},$node unless defined $exclude_tg_branches{$name};
		}
		if(defined $j){
			$ra_tg_free_mut->[$j]=$ra_bg_free_mut->[$i];
			$ra_bg_free_mut->[$j]=$ra_tg_free_mut->[$i];
			$done[$j]=1;
		}
		$done[$i]=1;
	}
}

1;