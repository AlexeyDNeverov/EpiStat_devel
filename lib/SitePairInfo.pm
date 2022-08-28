package SitePairInfo;
use Class::Struct;
use SubstPairInfo;

struct SitePairInfo => {
		bg_site => '$',
		tg_site => '$',
		epistat => '$',
		npair_subst => '$',
		subst_info => '@'
};

sub add_subst_info{
	my $site_pair_info=shift;
	my($bg_node,$tg_node,$count)=@_;
	my $subst_info=SubstPairInfo->new();
	$subst_info->tg_branch($tg_node);
	$subst_info->bg_branch($bg_node);
	$subst_info->count($count);
	$site_pair_info->npair_subst($site_pair_info->npair_subst+$count);
	push @{$site_pair_info->subst_info},$subst_info;
}
	

1;