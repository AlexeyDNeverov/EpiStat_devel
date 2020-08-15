package Pairs::Stat;
use SubstPairInfo;
use EpistatAnalysis::PhyloDistances qw(calc_subst_distance);
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
# Symbols to autoexport (:DEFAULT tag)
@EXPORT = qw(stat_func); 
@EXPORT_OK = qw(stat_func);

use strict;

sub stat_func {
    my $code = \&{shift @_}; # ensure we have something like CODE
	my ($subst_pair_info,$tau,$w1,$w2)=@_;
	my ($bg_node,$tg_node)=($subst_pair_info->bg_branch,$subst_pair_info->tg_branch);
	my $t=calc_subst_distance($bg_node,$tg_node);
	my $count=$subst_pair_info->count();
	my $esi=$count;
	$esi*=$code->($t,$tau);
	if(defined($w1)&&defined($w2)){
		$esi*=$w1*$w2;
	}
	return $esi;
}

1;