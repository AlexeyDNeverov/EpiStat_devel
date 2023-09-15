#!/usr/bin/env perl
#This script renerates a random distribution of mutations on the tree branches according to allele- and site-specific substitition rates.
#Params:
#<parameters> - The file with script params.
#		Contents:
#		XPAR="FN" - input file for EpiStat application [Kryazhimsky11]. Contains tree in the adjacent file format
#		PairsType=("intra[gene]"|"inter[gene]")
#		MatrixMode="1|2|3" - mode of mutation matrix  (see below): 1 - background;2 - target; 3 - both
#		Alphabet="STR" - the string of ordered one letter symbols for allele states in mutation codes: <state1><site><state2>.
#			!!!This parameter is required to generate alleles for mutations in the output XPARR file
#		BranchDistUnits="(DEFAULT||SYN|NSYN|TOTAL)" - Units for branch length distances
#		[BgrTreatProfile="FN"] - The profile is used to generate unidentified alleles which codes in mutations from the XPARR file are absent in the alphabet.
#			The IQTree *.sitefreq format is expected!
#			If 'TreatProfile' is not specified, the program will be aborted if unidentified alleles were in the XPARR file.
#		[BgrRootSequence="PROTEIN_SEQ"] - the string with protein sequence in the root of the tree specified in XPARR file.
#			This sequence is used to initialize alleles. If not specified the root sequence is generated from the XPARR file.
#		[FgrTreatProfile="FN"] - The profile is used to generate unidentified alleles which codes in mutations from the XPARR file are absent in the alphabet.
#			The IQTree *.sitefreq format is expected!
#			If 'TreatProfile' is not specified, the program will be aborted if unidentified alleles were in the XPARR file.
#		[FgrRootSequence="PROTEIN_SEQ"] - the string with protein sequence in the root of the tree specified in XPARR file.
#			This sequence is used to initialize alleles. If not specified the root sequence is generated from the XPARR file.

use strict;
#server inst
use lib "$ENV{EPISTAT_LIB}";
use Bio::Phylo::IO;
use Getopt::Std;
use File::Basename;
use SiteModel::ProteinMixtureModelProfile;
use SiteModel::EmpiricalProfile;
use SiteModel::AlleleMutRates;
use IO::XPARR qw(parse_xparr);

my $xpar_fn;
my $f_intragene;
my $f_mtx_mode;
my $distance_mode;
my $xparr_header_str;
#####################
my $brsubsets_num=1;
my @subset_folders;
#allele distribution
my $alphabet_str;
my @alphabet;
my $ra_alphabet;
my $bgr_root_seq;
my $bgr_treat_prof_fn;
my $fgr_root_seq;
my $fgr_treat_prof_fn;
my $rh_allele2idx;
my $f_treat_allleles=2; #ignore mutations with undefined alleles
my $fgr_treat_prof;
my $bgr_treat_prof;

open INFILE, "$ARGV[0]" or die "\nUnable to open input file: $ARGV[0]!";
while(<INFILE>){
	if(/#/){
		my $str1=$`;
		my $str2=$';
		$_=$str1 unless($str1=~/\"[^\"]*$/&&$str2=~/^[^\"]*\"/);
	}
	if(/(\w+)\s*=\s*\"\s*(.+?)\s*\"/){
		my ($key,$value)=($1,$2);
		if($key eq "XPAR"){
			$xpar_fn=$value;
		}elsif($key eq "MatrixMode"){
			$f_mtx_mode=$value;
			die "\nUnknown value ($value) for the PrintMatrixMode parameter in file $ARGV[0]!\n\tUINT in range is [1..3] expected!" unless $value=~/^[123]\s*$/;
		}elsif($key eq "PairsType"){
			if($value=~m/^intra(gene)*/i){
				$f_intragene=1;
			}elsif($value=~m/^inter(gene)*/i){
				$f_intragene=0;
			}else{
				die "\nUnknown value ($value) for the PairsType parameter in file $ARGV[0]!"; 
			}
		}elsif($key eq "Alphabet"){
			$alphabet_str= uc $value;
			$alphabet_str=~s/\s//;
		}elsif($key eq "BranchDistUnits"){
			if($value=~m/^default/i){
				$distance_mode=0;
			}elsif($value=~m/^syn/i){
				$distance_mode=1;
			}elsif($value=~m/^nsyn/i){
				$distance_mode=2;
			}elsif($value=~m/^total/i){
				$distance_mode=3;
			}else{
				die "\nUnknown value for 'BranchDistUnits'!";
			}
		}elsif($key eq "BgrTreatProfile"){
			$bgr_treat_prof_fn=$value;
		}elsif($key eq "BgrRootSequence"){
			$bgr_root_seq=$value;
			$bgr_root_seq=~s/\s//;
		}elsif($key eq "FgrTreatProfile"){
			$fgr_treat_prof_fn=$value;
		}elsif($key eq "FgrRootSequence"){
			$fgr_root_seq=$value;
			$fgr_root_seq=~s/\s//;
		}else{
			die "\nWrong parameter: $key!";
		}
	}
}
close INFILE;
die "\nError: Parameter PairsType is undefined in the file $ARGV[0]!" unless defined $f_intragene;
die "\nThe alphabet string isn't defined!" unless defined $alphabet_str;
die "\nThe units of the branch length measure are not specified!" unless defined $distance_mode;
if($f_intragene){
	if($f_mtx_mode==3){
		if(defined($bgr_treat_prof_fn)&&defined($fgr_treat_prof_fn)){
			die "\nThe same profiles are expected for the intragene matrix mode=3!" unless $bgr_treat_prof_fn eq $fgr_treat_prof_fn;
			$fgr_treat_prof_fn=undef;
		}elsif(defined($fgr_treat_prof_fn)){
			$bgr_treat_prof_fn=$fgr_treat_prof_fn;
			$fgr_treat_prof_fn=undef;
		}
		if(defined($bgr_root_seq)&&defined($fgr_root_seq)){
			die "\nThe same root sequences are expected for the intragene matrix mode=3!" unless length($bgr_root_seq)==length($fgr_root_seq);
			warn "\nYou've specified root sequences for the background as well as for the foreground. The sequence for the background will be ignored!";
			$fgr_root_seq=undef;
		}elsif(defined($fgr_root_seq)){
			$bgr_root_seq=$fgr_root_seq;
			$fgr_root_seq=undef;
		}
	}else{
		if(defined($bgr_treat_prof_fn)||defined($fgr_treat_prof_fn)){
			if($f_mtx_mode==2){
				die "\nThe profile is required to be specified either for both of background and foreground or don't specified at all!" unless defined $fgr_treat_prof_fn;
			}else{
				die "\nThe profile is required to be specified either for both of background and foreground or don't specified at all!" unless defined $bgr_treat_prof_fn;
			}
		}
		if(defined($bgr_root_seq)||defined($fgr_root_seq)){
			if($f_mtx_mode==2){
				die "\nThe root sequence is required to be specified either for foreground or don't specified at all!" unless defined $fgr_root_seq;
				warn "\nYou've specified root sequences for the background as well as for the foreground. The sequence for the background will be ignored!" if defined $bgr_root_seq;
				$bgr_root_seq=undef;
			}else{
				die "\nThe root sequence is required to be specified either for background or don't specified at all!" unless defined $bgr_root_seq;
				warn "\nYou've specified root sequences for the background as well as for the foreground. The sequence for the foreground will be ignored!" if defined $fgr_root_seq;
				$fgr_root_seq=undef;
			}
		}
	}
}
@alphabet=split "", $alphabet_str;
$ra_alphabet=\@alphabet;
if(defined($bgr_treat_prof_fn)||defined($fgr_treat_prof_fn)){
	$f_treat_allleles=1;
	$bgr_treat_prof=SiteModel::ProteinMixtureModelProfile->new($bgr_treat_prof_fn,$ra_alphabet) if defined $bgr_treat_prof_fn;
	$fgr_treat_prof=SiteModel::ProteinMixtureModelProfile->new($fgr_treat_prof_fn,$ra_alphabet) if defined $fgr_treat_prof_fn;
	$fgr_treat_prof=$bgr_treat_prof if($f_intragene&&$f_mtx_mode==3);
}
my %allele2idx;
for(my $i=0;$i<@{$ra_alphabet};$i++){
	if(defined $allele2idx{$ra_alphabet->[$i]}){
		die "\nError: The alphabet string contains the character ".$ra_alphabet->[$i]." several times!";
	}else{
		$allele2idx{$ra_alphabet->[$i]}=$i;
	}
}

my $tree=Bio::Phylo::IO->parse(
    -file => $xpar_fn,
    -format => 'adjacency');
#!!!The return value of Adjacency parser not corresponded to interface of the Bio::Phylo::IO->parse
$tree=$tree->[0]->first;

#my $str=Bio::Phylo::IO->unparse(-phylo => $tree,-format => 'newick');
#print $str;
my %branch_length;
#substitution counts
my %bg_nsubst;
my %tg_nsubst;
 
#Original distribution of substitutions on branches
my %tg_subst_map;
my %bg_subst_map;
my %syn_subst_map; #synonimous substitutions
my %phenotype_map;
my $f_phenotypes;
my %bg_site2idx; #converts background site coordinate into array index
my @bg_idx2site;
my %tg_site2idx; #converts target site coordinate into array index
my @tg_idx2site;

#New distribution of substitutions on branches
my $rh_smpl_bg_subst_map;
my $rh_smpl_tg_subst_map;

sub sample_base_from_profile{
	my ($ra_prof)=@_;
	my @prof=@{$ra_prof};
	my $n=@prof;
	for(my $j=1;$j<$n;$j++){
		$prof[$j]+=$prof[$j-1];
	}
	my $p=rand;
	my $k=0;
	for(;$k<$n;$k++){
		last if $prof[$k]>=$p;
	}
	return $k;
}

#translates characters into allele indices
sub init_sequence{
	my %args=@_;
	my $ra_out_seq=$args{'-result'};
	die "\nThe array reference to output sequence is required: specify '-result'!" unless defined $ra_out_seq;
	@{$ra_out_seq}=();
	if(defined $args{'-sequence_str'}){
		my $ra_alphabet=$args{'-alphabet'};
		die "\nThe reference to an array of ordered states is required: specify '-alphabet'!" unless defined $ra_alphabet;
		my %state2idx;
		for(my $i=0;$i<@{$ra_alphabet};$i++){
			$state2idx{$ra_alphabet->[$i]}=$i;
		}
		my @sequence=split "",$args{'-sequence_str'};
		foreach my $base(@sequence){
			my $base_idx=$state2idx{$base};
			die "\nThe symbol $base is not in the alphabet:\n\t".join "",@{$ra_alphabet} unless defined $base_idx;
			push @{$ra_out_seq},$base_idx;
		}
		return 1;
	}elsif(defined $args{'-profile'}){
		my $profile=$args{'-profile'};
		@{$ra_out_seq}=$profile->get_root_seq;
		return 1;
	}
	return 0;
}

sub init_alignment{
	my ($tree,$ra_root_seq,$allele_mut_rates,$profile,$ra_idx2site,$orh_align,$orh_subst_map)=@_;
	%{$orh_align}=();
	%{$orh_subst_map}=();
	$orh_align->{$tree->get_root()->get_name}=$ra_root_seq;
	my $alen=@{$ra_root_seq};
	$tree->visit_breadth_first(
		-in   => sub {
			my $node=shift;
			if(!$node->is_root){
				my $name=$node->get_name;
				my $pnode=$node->get_parent;
				my $pname=$pnode->get_name;
				my $ra_seq=[@{$orh_align->{$pname}}];
				my $l=$node->get_branch_length;
				for(my $i=0;$i<$alen;$i++){
					my $r=$allele_mut_rates->get_mutation_rate($i,$ra_seq->[$i]);
					my $p=exp(-$r*$l);
					if(rand()>=$p){
						my @prof=$profile->get_site_profile($i,$ra_seq->[$i]);
						my $k=sample_base_from_profile(\@prof);
						$ra_seq->[$i]=$k;
						my $site=$ra_idx2site->[$i];
						$orh_subst_map->{$name}->{$site}=1;
					}
				}
				$orh_align->{$name}=$ra_seq;
			}
		}
	);
	return 1;
}

sub unparse_subst1{
	my ($ra_sites,$ra_alphabet)=@_;
	my @sites=sort {$a->site()<=>$b->site()} @{$ra_sites};
	my @bases=@{$sites[0]->bases};
	@bases=($ra_alphabet->[$bases[0]],$ra_alphabet->[$bases[1]]) if(defined $ra_alphabet);
	my $str=$bases[0].$sites[0]->site().$bases[1];
	for(my $i=1;$i<@sites;$i++){
		@bases=@{$sites[$i]->bases};
		@bases=($ra_alphabet->[$bases[0]],$ra_alphabet->[$bases[1]]) if(defined $ra_alphabet);
		$str.=";".$bases[0].$sites[$i]->site().$bases[1];
	}
	return $str;
}

sub unparse_subst2{
	my $rh_align=shift;
	my $name=shift;
	my $pname=shift;
	my $rh_site2idx=shift;
	my $ra_alphabet=shift;
	my $ra_sites=shift;
	my @sites=sort {$a<=>$b} @{$ra_sites};
	my $str="";
	for(my $i=0;$i<@sites;$i++){
		my $site=$sites[$i];
		my $site_idx=$rh_site2idx->{$site};
		my @bases=($rh_align->{$pname}->[$site_idx],$rh_align->{$name}->[$site_idx]);
		@bases=($ra_alphabet->[$bases[0]],$ra_alphabet->[$bases[1]]);
		die "\nError unparse_subst(): for the node $name the mutation in the site $site_idx is declared but not found:!" if $rh_align->{$pname}->[$site_idx]==$rh_align->{$name}->[$site_idx];
		$str.=$bases[0].$site.$bases[1];
		$str.=";" unless $i==$#sites;
	}
	return $str;
}	
	
#begin script
my $rh_bgr_align;
my $rh_fgr_align;
my $rh_bg_site2idx;
my $rh_tg_site2idx;

{
	my $rh_syn_ncounts;
	my $rh_nsyn_ncounts;
	if($distance_mode==1||$distance_mode==3){
		$rh_syn_ncounts={};
	}
	if($distance_mode==2||$distance_mode==3){
		$rh_nsyn_ncounts={};
	}
	$xparr_header_str=parse_xparr($xpar_fn,"-f_treat_alleles" => $f_treat_allleles,"-bgr_treat_profile" => $bgr_treat_prof,"-fgr_treat_profile" => $fgr_treat_prof,
		"-allele_indices" => \%allele2idx,
		"-orh_bgr_site_nsubst" => \%bg_nsubst,"-orh_fgr_site_nsubst" => \%tg_nsubst,
		"-orh_bgr_subst_map" => \%bg_subst_map, "-orh_fgr_subst_map" => \%tg_subst_map,
		"-orh_bgr_site2idx" => \%bg_site2idx, "-orh_fgr_site2idx" => \%tg_site2idx,
		"-ora_bgr_idx2site" => \@bg_idx2site, "-ora_fgr_idx2site" => \@tg_idx2site,
		"-orh_syn_counts" => $rh_syn_ncounts, "-orh_fgr_nsyn_counts" => $rh_nsyn_ncounts,
		"-orh_column_str=3" => \%syn_subst_map, "-orh_column_str=6" => \%phenotype_map
	);
	$f_phenotypes=($xparr_header_str=~m/phenotypes:/);
	if($distance_mode){
		foreach my $node($tree->get_nodes){
			my $name=$node->get_name();
			my $l=$node->get_branch_length;
			$branch_length{$name}=$l;
			$l=0;
			if($distance_mode==1||$distance_mode==3){
				$l+=$rh_syn_ncounts->{$name};
			}
			if($distance_mode==2||$distance_mode==3){
				$l+=$rh_nsyn_ncounts->{$name};
			}
			$node->set_branch_length($l);
		}
	}
}

if($f_mtx_mode!=2){
	$rh_smpl_bg_subst_map={};
	my $ra_bgr_root_seq=[];
	my $bgr_profile=SiteModel::EmpiricalProfile->new($tree,\%bg_subst_map,scalar(@alphabet),\%bg_site2idx,1);
	my $bgr_arates=SiteModel::AlleleMutRates->new($tree,\%bg_subst_map,scalar(@alphabet),\%bg_site2idx);
#DEBUG: Print allele rate matrix
#$bgr_arates->print;
#exit;
##############
	if(defined $bgr_root_seq){
		init_sequence('-result' => $ra_bgr_root_seq,'-sequence_str' =>$bgr_root_seq,'-alphabet'=>$ra_alphabet);
	}else{
		init_sequence('-result' => $ra_bgr_root_seq,'-profile' => $bgr_arates);
	}
	$rh_bgr_align={};
	init_alignment($tree,$ra_bgr_root_seq,$bgr_arates,$bgr_profile,\@bg_idx2site,$rh_bgr_align,$rh_smpl_bg_subst_map);
}

if($f_mtx_mode!=1){
	unless($f_intragene&&$f_mtx_mode==3){
		$rh_smpl_tg_subst_map={};
		my $ra_fgr_root_seq=[];
		my $fgr_profile=SiteModel::EmpiricalProfile->new($tree,\%tg_subst_map,scalar(@alphabet),\%tg_site2idx,1);
		my $fgr_arates=SiteModel::AlleleMutRates->new($tree,\%tg_subst_map,scalar(@alphabet),\%tg_site2idx);
		if(defined $fgr_root_seq){
			init_sequence('-result' => $ra_fgr_root_seq,'-sequence_str' =>$fgr_root_seq,'-alphabet'=>$ra_alphabet);
		}else{
			init_sequence('-result' => $ra_fgr_root_seq,'-profile' => $fgr_arates);
		}
		$rh_fgr_align={};
		init_alignment($tree,$ra_fgr_root_seq,$fgr_arates,$fgr_profile,\@tg_idx2site,$rh_fgr_align,$rh_smpl_tg_subst_map);
	}else{
		$rh_fgr_align=$rh_bgr_align;
	}
}

#print XPARR
print $xparr_header_str;
foreach my $node($tree->get_nodes){
	my $name=$node->get_name;
	print "\n$name";
	if(!$node->is_root){
		my $pname=$node->get_parent->get_name;
		my $length=$node->get_branch_length;
		$length=$branch_length{$name} if($distance_mode);
		print "\t$pname\t$length\t$syn_subst_map{$name}\t";
		my @sites1;
		my @sites2;
		my @tmp;
		my $str="";
		{
			if(defined $rh_smpl_tg_subst_map){
				@sites2=keys %{$rh_smpl_tg_subst_map->{$name}};
			}elsif(defined $tg_subst_map{$name}){
				@tmp=keys %{$tg_subst_map{$name}};
			}
			if(!$f_intragene){
				foreach my $site(@tmp){
					push @sites1,$tg_subst_map{$name}->{$site};
				}
			}elsif(!defined($rh_smpl_tg_subst_map)){
				foreach my $site(@tmp){
					push @sites1,$tg_subst_map{$name}->{$site} unless defined $bg_nsubst{$site};
				}
				@tmp=keys %{$rh_smpl_bg_subst_map->{$name}};
				for(my $i=0;$i<@tmp;$i++){
					push @sites2,$tmp[$i] if defined $tg_nsubst{$tmp[$i]};
				}
			}
			$str.=unparse_subst1(\@sites1,$ra_alphabet) if @sites1;
			if(@sites2){
				$str.=";" if @sites1;
				if(defined $rh_fgr_align){
					$str.=unparse_subst2($rh_fgr_align,$name,$pname,\%tg_site2idx,$ra_alphabet,\@sites2);
				}else{
					$str.=join ";",sort {$a<=>$b} @sites2;
				}
			}
		}
		print "$str\t";
		@sites1=();
		@sites2=();
		@tmp=();
		$str="";
		{
			if(defined $rh_smpl_bg_subst_map){
				@sites2=keys %{$rh_smpl_bg_subst_map->{$name}};
			}else{
				@tmp=keys %{$bg_subst_map{$name}} if(defined $bg_subst_map{$name});
			}
			if(!$f_intragene){
				foreach my $site(@tmp){
					push @sites1,$bg_subst_map{$name}->{$site};
				}
			}elsif(!defined $rh_smpl_bg_subst_map){
				for(my $i=0;$i<@tmp;$i++){
					push @sites1,$tmp[$i] unless defined $tg_nsubst{$tmp[$i]};
				}
				@tmp=keys %{$rh_smpl_tg_subst_map->{$name}};
				for(my $i=0;$i<@tmp;$i++){
					push @sites2,$tmp[$i] if(defined $bg_nsubst{$tmp[$i]});
				}
			}
			$str.=unparse_subst1(\@sites1,$ra_alphabet) if @sites1;
			if(@sites2){
				$str.=";" if @sites1;
				if(defined $rh_bgr_align){
					$str.=unparse_subst2($rh_bgr_align,$name,$pname,\%bg_site2idx,$ra_alphabet,\@sites2);
				}else{
					$str.=join ";",sort {$a<=>$b} @sites2;
				}
			}
			print $str;
		}
		if($f_phenotypes){
			print "\t";
			print $phenotype_map{$name} if defined $phenotype_map{$name};
		}
	}
}
			