package BranchSiteMutationsMatrix;
#This module provides class for representing a phylogenetic tree with mutations on branches (XPARR) in the binary matrix form. 
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
# Symbols to autoexport (:DEFAULT tag)
@EXPORT = qw(print_branch_site_matrix); 
@EXPORT_OK = qw(print_branch_site_matrix);

sub print_branch_site_matrix{
	my ($out_fn_prefix,$f_mtx_mode,$f_intragene,$ra_bg_sites,$ra_tg_sites,$ra_bg_branches,$ra_tg_branches,
		$rh_bg_subst_map, $rh_tg_subst_map)=@_;
	my %bg_branches;
	my @bg_idx2branches;
	my %tg_branches;
	my @tg_idx2branches;
	my %bg_sites;
	my @bg_idx2sites;
	my %tg_sites;
	my @tg_idx2sites;
	my $mmtx_ext=".mmtx";
	my $idx2site_ext=".idx2site";
	my $idx2branch_ext=".idx2branch";
	if($f_mtx_mode!=2){
		foreach my $bname(@{$ra_bg_branches}){
			$bg_branches{$bname}=1;
		}
		for(my $i=0;$i<@{$ra_bg_sites};$i++){
			$bg_sites{$ra_bg_sites->[$i]}=$i;
		}
		unless($f_mtx_mode==3&&$f_intragene){
			@bg_idx2sites=sort @{$ra_bg_sites}
		}
	}
	if($f_mtx_mode!=1){
		if($f_mtx_mode==3&&$f_intragene){
			foreach my $spos(@{$ra_tg_sites}){
				$bg_sites{$spos}=1;
			}
			foreach my $bname(@{$ra_tg_branches}){
				$bg_branches{$bname}=1;
			}
			my $i=0;
			foreach my $spos(sort keys %bg_sites){
				$bg_idx2sites[$i]=$spos;
				$bg_sites{$spos}=$i++;
			}
		}else{
			for(my $i=0;$i<@{$ra_tg_sites};$i++){
				$tg_sites{$ra_tg_sites->[$i]}=$i;
			}
			@tg_idx2sites=sort @{$ra_tg_sites};
			foreach my $bname(@{$ra_tg_branches}){
				$tg_branches{$bname}=1;
			}
		}
	}
	if($f_mtx_mode==1||$f_mtx_mode==3){
		my $i=0;
		foreach my $bname(sort keys %bg_branches){
			$bg_idx2branches[$i]=$bname;
			$bg_branches{$bname}=$i++;
		}
		my $fn=$out_fn_prefix;
		if($f_intragene){
			$fn.=".intra";
		}else{
			$fn.=".inter";
		}
		my $str;
		my @branch_margins;
		my @site_margins;
		if($f_intragene&&$f_mtx_mode==3){
			$fn.=".both";
			$str=$fn.$mmtx_ext;
			_print_matrix($str,\@bg_idx2branches,\@bg_idx2sites,$rh_bg_subst_map,$rh_tg_subst_map,
				\@branch_margins,\@site_margins);
		}else{
			$fn.=".bgr";
			$str=$fn.$mmtx_ext;
			_print_matrix($str,\@bg_idx2branches,\@bg_idx2sites,$rh_bg_subst_map,undef,
				\@branch_margins,\@site_margins);
		}
		$str=$fn.$idx2branch_ext;
		_print_indices($str,\@bg_idx2branches,\@branch_margins);
		$str=$fn.$idx2site_ext;
		_print_indices($str,\@bg_idx2sites,\@site_margins);
	}
	if($f_mtx_mode==2||($f_mtx_mode==3&&$f_intragene==0)){
		my $i=0;
		foreach my $bname(sort keys %tg_branches){
			$tg_idx2branches[$i]=$bname;
			$tg_branches{$bname}=$i++;
		};
		my $fn=$out_fn_prefix;
		if($f_intragene){
			$fn.=".intra";
		}else{
			$fn.=".inter";
		}
		$fn.=".fgr";
		my @branch_margins;
		my @site_margins;
		my $str=$fn.$mmtx_ext;
		_print_matrix($str,\@tg_idx2branches,\@tg_idx2sites,undef,$rh_tg_subst_map,
			\@branch_margins,\@site_margins);
		$str=$fn.$idx2branch_ext;
		_print_indices($str,\@tg_idx2branches,\@branch_margins);
		$str=$fn.$idx2site_ext;
		_print_indices($str,\@tg_idx2sites,\@site_margins);
	}
}

sub _print_matrix{
	my ($fn,$ra_idx2branches,$ra_idx2sites,$rh_bg_subst_map,$rh_tg_subst_map,
		$ra_branch_margins,$ra_site_margins)=@_;
	open OPF, ">$fn" or die "\nUnable to open output file $fn!";
	for(my $i=0;$i<@{$ra_idx2branches};$i++){
		my $str;
		my $bname=$ra_idx2branches->[$i];
		for (my $j=0;$j<@{$ra_idx2sites};$j++){
			my $t=0;
			my $spos=$ra_idx2sites->[$j];
			if(defined $rh_bg_subst_map){
				if(defined $rh_bg_subst_map->{$bname}){
					$t=1 if defined $rh_bg_subst_map->{$bname}->{$spos};
				};
			};
			if(defined $rh_tg_subst_map){
				if(defined $rh_tg_subst_map->{$bname}){
					$t=1 if defined $rh_tg_subst_map->{$bname}->{$spos};
				};
			};		
			$str.=$t;
			$ra_site_margins->[$j]+=$t if defined $ra_site_margins;
			$ra_branch_margins->[$i]+=$t if defined $ra_branch_margins;
		};
		if(defined $ra_branch_margins){
			warn "\nIn matrix $fn\n\trow $i is empty!" unless $ra_branch_margins->[$i];
		};
		print OPF "$str\n";
	};
	if(defined $ra_site_margins){
		for (my $j=0;$j<@{$ra_idx2sites};$j++){
			warn "\nIn matrix $fn\n\tcolumn $j is empty!" unless $ra_site_margins->[$j];
		}
	};
	close OPF;
}

sub _print_indices{
	my ($fn,$ra_indices,$ra_margins)=@_;
	open OPF, ">$fn" or die "\nUnable to open output file $fn!";
	for(my $i=0;$i<@{$ra_indices};$i++){
		my $str=$ra_indices->[$i];
		print OPF "$str";
		if(defined $ra_margins){
			$str=$ra_margins->[$i];
			print OPF "\t$str";
		};
		print OPF "\n";
	};
	close OPF;
}

1;