package IO::EmpiricalMutationModel;
#This module provides class for representing a phylogenetic tree with mutations on branches (XPARR) in the binary matrix form. 
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
# Symbols to autoexport (:DEFAULT tag)
@EXPORT = qw(print_branch_site_mutation_matrix 
		print_allele_mutation_matrix read_allele_mutation_matrix
		print_align_profile read_align_profile
		print_mutation_rate_matrix); 
@EXPORT_OK = qw(print_branch_site_mutation_matrix 
		print_allele_mutation_matrix read_allele_mutation_matrix
		print_align_profile read_align_profile
		print_mutation_rate_matrix);

use SiteModel::AlleleMutRates;

my $mmtx_ext=".mmtx";
my $idx2site_ext=".idx2site";
my $idx2branch_ext=".idx2branch";
my $amfreq_ext=".amfreq";
my $alnprof_ext=".aln_prof";
my $amrates_ext=".amrates";

sub print_branch_site_mutation_matrix{
	my ($out_fn_prefix,$f_mtx_mode,$f_intragene,$ra_bg_sites,$ra_tg_sites,$ra_bg_branches,$ra_tg_branches,
		$rh_bg_subst_map, $rh_tg_subst_map)=@_;
	my %bg_branches;
	my @bg_idx2branches;
	my %tg_branches;
	my @tg_idx2branches;
	my %bg_sites;
	my @bg_idx2sites;
	my @tg_idx2sites;
	if($f_mtx_mode!=2){
		foreach my $bname(@{$ra_bg_branches}){
			$bg_branches{$bname}=1;
		}
		unless($f_mtx_mode==3&&$f_intragene){
			@bg_idx2sites=sort @{$ra_bg_sites}
		}else{
			foreach my $spos(@{$ra_bg_sites}){
				$bg_sites{$spos}=1;
			}
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
				$bg_idx2sites[$i++]=$spos;
			}
		}else{
			@tg_idx2sites=sort @{$ra_tg_sites};
			foreach my $bname(@{$ra_tg_branches}){
				$tg_branches{$bname}=1;
			}
		}
	}
	my $fn=$out_fn_prefix;
	if($f_intragene){
		$fn.=".intra";
	}else{
		$fn.=".inter";
	}
	if($f_mtx_mode==1||$f_mtx_mode==3){
		my $i=0;
		foreach my $bname(sort keys %bg_branches){
			$bg_idx2branches[$i]=$bname;
			$bg_branches{$bname}=$i++;
		}
		my $str;
		my @branch_margins;
		my @site_margins;
		my $fn1=$fn;
		if($f_intragene&&$f_mtx_mode==3){
			$fn1.=".both";
			$str=$fn1.$mmtx_ext;
			_print_matrix($str,\@bg_idx2branches,\@bg_idx2sites,$rh_bg_subst_map,$rh_tg_subst_map,
				\@branch_margins,\@site_margins);
		}else{
			$fn1.=".bgr";
			$str=$fn1.$mmtx_ext;
			_print_matrix($str,\@bg_idx2branches,\@bg_idx2sites,$rh_bg_subst_map,undef,
				\@branch_margins,\@site_margins);
		}
		$str=$fn1.$idx2branch_ext;
		_print_indices($str,\@bg_idx2branches,\@branch_margins);
		$str=$fn1.$idx2site_ext;
		_print_indices($str,\@bg_idx2sites,\@site_margins);
	}
	if($f_mtx_mode==2||($f_mtx_mode==3&&$f_intragene==0)){
		my $i=0;
		foreach my $bname(sort keys %tg_branches){
			$tg_idx2branches[$i]=$bname;
			$tg_branches{$bname}=$i++;
		};
		my $fn2=$fn;
		$fn2.=".fgr";
		my @branch_margins;
		my @site_margins;
		my $str=$fn2.$mmtx_ext;
		_print_matrix($str,\@tg_idx2branches,\@tg_idx2sites,undef,$rh_tg_subst_map,
			\@branch_margins,\@site_margins);
		$str=$fn2.$idx2branch_ext;
		_print_indices($str,\@tg_idx2branches,\@branch_margins);
		$str=$fn2.$idx2site_ext;
		_print_indices($str,\@tg_idx2sites,\@site_margins);
	}
}

sub print_mutation_rate_matrix{
	my ($out_fn_prefix,$tree,$alpha_size,$f_mtx_mode,$f_intragene,$ra_bg_sites,$ra_tg_sites,$rh_bg_subst_map,$rh_tg_subst_map)=@_;
	my $fn=$out_fn_prefix;
	if($f_intragene){
		$fn.=".intra";
	}else{
		$fn.=".inter";
	}
	if($f_mtx_mode==1||$f_mtx_mode==3){
		my $fn1=$fn;
		if($f_intragene&&$f_mtx_mode==3){
			$fn1.=".both";
		}else{
			$fn1.=".bgr";
		}
		$fn1.=$amrates_ext;
		open OPF, ">$fn1" or die "\nUnable to open output file $fn1!";
		my $amrates=SiteModel::AlleleMutRates->new($tree,$ra_bg_sites,$rh_bg_subst_map,$alpha_size,1);
		$amrates->print(*OPF);
		close OPF;
	}
	if($f_mtx_mode==2||($f_mtx_mode==3&&$f_intragene==0)){
		my $fn2=$fn;
		$fn2.=".fgr";
		$fn2.=$amrates_ext;
		open OPF, ">$fn2" or die "\nUnable to open output file $fn2!";
		my $amrates=SiteModel::AlleleMutRates->new($tree,$ra_tg_sites,$rh_tg_subst_map,$alpha_size,1);
		$amrates->print(*OPF);
		close OPF;
	}
}

sub print_allele_mutation_matrix{
	my ($out_fn_prefix,$f_mtx_mode,$f_intragene,$allele_site_model,$hr_bg_site2idx,$hr_tg_site2idx)=@_;
	my $fn=$out_fn_prefix;
	if($f_intragene){
		$fn.=".intra";
	}else{
		$fn.=".inter";
	}
	if($f_mtx_mode==1||$f_mtx_mode==3){
		my $fn1=$fn;
		if($f_intragene&&$f_mtx_mode==3){
			$fn1.=".both";
		}else{
			$fn1.=".bgr";
		}
		$fn1.=$amfreq_ext;
		open OPF, ">$fn1" or die "\nUnable to open output file $fn1!";
		$allele_site_model->print_mutation_profile(*OPF,1,$hr_bg_site2idx);
		close OPF;
	}
	if($f_mtx_mode==2||($f_mtx_mode==3&&$f_intragene==0)){
		my $fn2=$fn;
		$fn2.=".fgr";
		$fn2.=$amfreq_ext;
		open OPF, ">$fn2" or die "\nUnable to open output file $fn2!";
		$allele_site_model->print_mutation_profile(*OPF,2,$hr_tg_site2idx);
		close OPF;
	}
}

#site and allele indices are numbered from 1.
sub read_allele_mutation_matrix{
	my ($amfreg_fn,$rao_mut_profile,$rao_idx2site)=@_;
	@{$rao_mut_profile}=() if defined $rao_mut_profile;
	@{$rao_idx2site}=() if defined $rao_idx2site;
	open INPF, "<$amfreg_fn" or die "\nUnable to open input file $amfreg_fn!";
	my $sites_col;
	while(<INPF>){
		chomp;
		next unless /\S+/;
		my @line=split "\t";
		if($line[0] eq "site_idx"){
			$sites_col=1 if $line[1] eq "site";
			die "\nThe input file $amfreg_fn is not a valid site-specific aminoacid mutation matrix!" unless ($line[-2] eq "anc_allele_idx")&&
					($line[-1] eq "mutation_freqs");
			last;
		}
	}
	die "\nWrong format of allele mutation matrix: $amfreg_fn!" if eof INPF;
	my $n=0;
	while(<INPF>){
		chomp;
		next unless /\S+/;
		my @line=split "\t",$_,-1;
		die "\nUnexpected value type!" unless $line[0]=~/^\d+/;
		my $i=1;
		my $site_idx=$line[0]-1;
		if(defined $sites_col){
			$i++;
			$rao_idx2site->[$site_idx]=$line[$sites_col] if defined $rao_idx2site;
		}
		my $anc_allele_idx=$line[$i]-1;
		$rao_mut_profile->[$site_idx]={} unless defined $rao_mut_profile->[$site_idx];
		$rao_mut_profile->[$site_idx]->{$anc_allele_idx}={};
		for(my $j=$i+1;$j<@line;$j++){
			my $der_allele_idx=$j-($i+1);
			$rao_mut_profile->[$site_idx]->{$anc_allele_idx}->{$der_allele_idx}=$line[$j] if $line[$j] ne "";
		}
		$n++;
	}
	close INPF;
	return $n;
}

#site indices are numbered from 1.
sub print_align_profile{
	my ($out_fn_prefix,$f_mtx_mode,$f_intragene,$allele_site_model,$hr_bg_site2idx,$hr_tg_site2idx)=@_;
	my $fn=$out_fn_prefix;
	if($f_intragene){
		$fn.=".intra";
	}else{
		$fn.=".inter";
	}
	if($f_mtx_mode==1||$f_mtx_mode==3){
		my $fn1=$fn;
		if($f_intragene&&$f_mtx_mode==3){
			$fn1.=".both";
		}else{
			$fn1.=".bgr";
		}
		$fn1.=$alnprof_ext;
		open OPF, ">$fn1" or die "\nUnable to open output file $fn1!";
		$allele_site_model->print_align_profile(*OPF,1,$hr_bg_site2idx);
		close OPF;
	}
	if($f_mtx_mode==2||($f_mtx_mode==3&&$f_intragene==0)){
		my $fn2=$fn;
		$fn2.=".fgr";
		$fn2.=$alnprof_ext;
		open OPF, ">$fn2" or die "\nUnable to open output file $fn2!";
		$allele_site_model->print_align_profile(*OPF,2,$hr_tg_site2idx);
		close OPF;
	}
}

#site indices are numbered from 1.
sub read_align_profile{
	my ($alnprof_fn,$rao_aln_profile,$rao_idx2site)=@_;
	@{$rao_aln_profile}=() if defined $rao_aln_profile;
	@{$rao_idx2site}=() if defined $rao_idx2site;
	open INPF, "<$alnprof_fn" or die "\nUnable to open input file $alnprof_fn!";
	my $sites_col;
	while(<INPF>){
		chomp;
		next unless /\S+/;
		my @line=split "\t";
		if($line[0] eq "site_idx"){
			$sites_col=1 if $line[1] eq "site";
			die "\nThe input file $alnprof_fn is not a valid alignment profile!" unless $line[-1] eq "allele_freqs";
			last;
		}
	}
	die "\nWrong format of allele mutation matrix: $alnprof_fn!" if eof INPF;
	my $n=0;
	while(<INPF>){
		chomp;
		next unless /\S+/;
		my @line=split "\t",$_,-1;
		die "\nUnexpected value type!" unless $line[0]=~/^\d+/;
		my $i=1;
		my $site_idx=$line[0]-1;
		if(defined $sites_col){
			$i++;
			$rao_idx2site->[$site_idx]=$line[$sites_col] if defined $rao_idx2site;
		}
		$rao_aln_profile->[$site_idx]={} unless defined $rao_aln_profile->[$site_idx];
		for(my $j=$i;$j<@line;$j++){
			my $allele_idx=$j-$i;
			$rao_aln_profile->[$site_idx]->{$allele_idx}=$line[$j] if $line[$j] ne "";
		}
		$n++;
	}
	close INPF;
	return $n;
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