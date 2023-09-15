#!/usr/bin/env perl
#This script converts matrix containing distribution of mutations on a tree branches into XPARR format.
#[options:]
#	-b <file name> - Sets BgrMutMatrix parameter
#	-t <file name> - Sets FgrMutMatrix parameter
#Params:
#<parameters> - The file with script params.
#		Contents:
#		XPAR="FN" - input file for EpiStat application [Kryazhimsky11]. Contains tree in the adjacent file format
#		PairsType=("intra[gene]"|"inter[gene]")
#		MatrixMode="1|2|3" - mode of mutation matrix  (see below): 1 - background;2 - target; 3 - both
#		IgnoreTerminals="0|1" - Ignore terminal branches
#			Default=1 - Ignore
###SECTION null model with alleles{
#		[Alphabet="STR"] - the string of ordered one letter symbols for allele states in mutation codes: <state1><site><state2>.
#			!!!This parameter is required to generate alleles for mutations in the output XPARR file
#		[BgrTreatProfile="FN"] - The profile is used to generate unidentified alleles which codes in mutations from the XPARR file are absent in the alphabet.
#			The IQTree *.sitefreq format is expected!
#			If 'TreatProfile' is not specified, the program will be aborted if unidentified alleles were in the XPARR file.
#		[BgrRootSequence="PROTEIN_SEQ"] - the string with protein sequence in the root of the tree that is specified in XPARR file.
#			This sequence is used to initialize alleles.
#		[FgrTreatProfile="FN"] - The profile is used to generate unidentified alleles which codes in mutations from the XPARR file are absent in the alphabet.
#			The IQTree *.sitefreq format is expected!
#			If 'TreatProfile' is not specified, the program will be aborted if unidentified alleles were in the XPARR file.
#		[FgrRootSequence="PROTEIN_SEQ"] - the string with protein sequence in the root of the tree that is specified in XPARR file.
#			This sequence is used to initialize alleles.
#		[MkRootSequenceFrom=("xparr"|"align(ment)?")] - specifies how to generate root sequences if either of 'BgrRootSequence' or 'FgrRootSequence' parameters has not been specified.
#			"xparr" - root sequences are generated from the mutation distribution on the tree branches.
#			"align" - root sequences are sampled from the alignment profile.
###}END SECTION
#		[NBranchSubsets="1|2"] - Number of subsets of branches for which distributions of mutations were defined separately.
#			!!!Note that for each subset binary matrices of mutations in the background or in the foreground are required
#			Default="1"
#		[SubsetFolder="FOLDER"] - The folder containing data of a subset. If specified used as a path prefix for the corresponding mutations matrix.
#			!!!The number of specifications is required to be equal to the number of subsets or to the zero.
#			Default=""
#		[BgrMutMatrix="FN"] - Distribution mutations in background. Binary matrix: raws - branches; columns - sites
#			!!!If defined, the number of specifications is required to be equal to the number of subsets or to the one.
#			!!!In the case of single specification the matrix with the same file name is expected in each subset folder.
#			Default=undef
#			!!!Option -b has a priority
#		BgrIdx2Site="FN" - Distribution mutations in background. Table to convert columns' indexes of MutMatrix into site positions
#			!!!If defined, the number of specifications is required to be equal to the number of subsets
#			Default=undef
#		BgrIdx2Branch="FN" - Distribution mutations in background. Table to convert raws' indexes of MutMatrix into branch names
#			!!!If defined, the number of specifications is required to be equal to the number of subsets
#			Default=undef
#		[FgrMutMatrix="FN"] - Distribution mutations in foreground. Binary matrix: raws - branches; columns - sites
#			!!!If defined, the number of specifications is required to be equal to the number of subsets or to the one.
#			!!!In the case of a single specification the matrix with the same file name is expected in each subset folder.
#			Default=undef
#			!!!Option -t has a priority
#		FgrIdx2Site="FN" - Distribution mutations in foreground. Table to convert columns' indexes of MutMatrix into site positions
#			!!!If defined, the number of specifications is required to be equal to the number of subsets
#			Default=undef
#		FgrIdx2Branch="FN" - Distribution mutations in foreground. Table to convert raws' indexes of MutMatrix in
#			!!!If defined, the number of specifications is required to be equal to the number of subsets
#			Default=undef

use strict;
use lib "$ENV{EPISTAT_LIB}";
use Bio::Phylo::IO;
use Class::Struct;
use Getopt::Std;
use File::Basename;
use SiteModel::ProteinMixtureModelProfile;
use SiteModel::EmpiricalProfile;
use IO::EmpiricalMutationModel qw(read_allele_mutation_matrix read_align_profile);
use IO::XPARR qw(parse_xparr);

my %args;
if(!getopts('b:t:',\%args)){
	die "\nError in option string!";
}

my $xpar_fn;
my @bg_mmtx_fn;
push @bg_mmtx_fn,split(',',$args{b}) if defined $args{b};
my @tg_mmtx_fn;
push @tg_mmtx_fn,split(',',$args{t}) if defined $args{t};
my @bg_idx2site_fn;
my @bg_idx2branch_fn;
my @tg_idx2site_fn;
my @tg_idx2branch_fn;
my $f_intragene;
my $f_mtx_mode;
my $f_ignore_terminals=1;
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
my $f_mk_root_seq;

my $rh_allele2idx;
my $f_treat_allleles=0;
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
		}elsif($key eq "IgnoreTerminals"){
			$f_ignore_terminals=$value;
			die "\nUnknown value ($value) for the IgnoreTerminals parameter in file $ARGV[0]!\n\tUINT in range is [0..1] expected!" unless $value=~/^[01]\s*$/;
		}elsif($key eq "BgrMutMatrix"){
			$value=~s/^\s+//;
			push @bg_mmtx_fn,$value unless defined $args{b};
		}elsif($key eq "BgrIdx2Site"){
			push @bg_idx2site_fn,$value;
		}elsif($key eq "BgrIdx2Branch"){
			push @bg_idx2branch_fn,$value;
		}elsif($key eq "FgrMutMatrix"){
			$value=~s/^\s+//;
			push @tg_mmtx_fn,$value unless defined $args{t};
		}elsif($key eq "FgrIdx2Site"){
			push @tg_idx2site_fn,$value;
		}elsif($key eq "FgrIdx2Branch"){
			push @tg_idx2branch_fn,$value;
		}elsif($key eq "PairsType"){
			if($value=~m/^intra(gene)*/i){
				$f_intragene=1;
			}elsif($value=~m/^inter(gene)*/i){
				$f_intragene=0;
			}else{
				die "\nUnknown value ($value) for the PairsType parameter in file $ARGV[0]!"; 
			}
		}elsif($key eq "NBranchSubsets"){
			$value=~s/^\s+//;
			if($value=~/^\d+\s*/){
				$brsubsets_num=$value;
			}else{
				die "\nUnpropper value for the 'NBranchSubsets' parameter: $value! UINT is expected.";
			}
		}elsif($key eq "SubsetFolder"){
			$value=~s/\s+$//;
			$value=~s/\/$//;
			$value.="/";
			push @subset_folders,$value;
		}elsif($key eq "Alphabet"){
			$alphabet_str= uc $value;
			$alphabet_str=~s/\s//;
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
		}elsif($key eq "MkRootSequenceFrom"){
			if($value=~/^xparr/i){
				$f_mk_root_seq="xparr";
			}elsif($value=~/^align/i){
				$f_mk_root_seq="align";
			}else{
				die "\nUnknown method for generation of root sequences: $value!";
			}
		}else{
			die "\nWrong parameter: $key!";
		}
	}
}
close INFILE;
die "\nError: Parameter PairsType is undefined in the file $ARGV[0]!" unless defined $f_intragene;
die "\nThe number of 'BgrIdx2Site' records must be equal to the number of branch subsets: $brsubsets_num!" unless @bg_idx2site_fn==0||$brsubsets_num==@bg_idx2site_fn;
die "\nThe number of 'BgrIdx2Branch' records must be equal to the number of branch subsets: $brsubsets_num!" unless @bg_idx2branch_fn==0||$brsubsets_num==@bg_idx2branch_fn;
die "\nThe number of 'FgrIdx2Site' records must be equal to the number of branch subsets: $brsubsets_num!" unless @tg_idx2site_fn==0||$brsubsets_num==@tg_idx2site_fn;
die "\nThe number of 'FgrIdx2Branch' records must be equal to the number of branch subsets: $brsubsets_num!" unless @tg_idx2branch_fn==0||$brsubsets_num==@tg_idx2branch_fn;
die "\nUndefined where has the mutations been perturbed in background or foreground?" unless @bg_mmtx_fn||@tg_mmtx_fn;
die "\nThe number of 'SubsetFolder' records must be equal to the number of branch subsets: $brsubsets_num!" unless @subset_folders==0||$brsubsets_num==@subset_folders;
if(defined($alphabet_str)){
	unless(defined($bgr_root_seq)&&defined($fgr_root_seq)){
		die "\nThe method to generate root sequence is not specified, use the 'MkRootSequenceFrom' parameter!" unless defined $f_mk_root_seq;
	}elsif(defined $f_mk_root_seq){
		warn "\nThe method of generating of the root seqequence is specified but will be ignored!";
		$f_mk_root_seq=undef;
	}
	@alphabet=split "", $alphabet_str;
	$ra_alphabet=\@alphabet;
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
				warn "\nYou specified root sequences for the background as well as for the foreground. The sequence for the background will be ignored!";
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
	$f_treat_allleles=2;
	if(defined($bgr_treat_prof_fn)||defined($fgr_treat_prof_fn)){
		$f_treat_allleles=1;
		$bgr_treat_prof=SiteModel::ProteinMixtureModelProfile->new($bgr_treat_prof_fn,$ra_alphabet) if defined $bgr_treat_prof_fn;
		$fgr_treat_prof=SiteModel::ProteinMixtureModelProfile->new($fgr_treat_prof_fn,$ra_alphabet) if defined $fgr_treat_prof_fn;
		$fgr_treat_prof=$bgr_treat_prof if $f_intragene&&$f_mtx_mode==3;
	}
	$rh_allele2idx={};
	for(my $i=0;$i<@{$ra_alphabet};$i++){
		if(defined $rh_allele2idx->{$ra_alphabet->[$i]}){
			die "\nError: The alphabet string contains the character ".$ra_alphabet->[$i]." several times!";
		}else{
			$rh_allele2idx->{$ra_alphabet->[$i]}=$i;
		}
	}
}else{
	if(defined $bgr_treat_prof_fn){
		warn "\nThe site profile $bgr_treat_prof_fn is specified but will be ignored!";
		$bgr_treat_prof_fn=undef;
	}
	if(defined $fgr_treat_prof_fn){
		warn "\nThe site profile $fgr_treat_prof_fn is specified but will be ignored!";
		$fgr_treat_prof_fn=undef;
	}
	if(defined $bgr_root_seq){
		warn "\nThe root sequence is specified but will be ignored!";
		$bgr_root_seq=undef;
	}
	if(defined $fgr_root_seq){
		warn "\nThe root sequence is specified but will be ignored!";
		$fgr_root_seq=undef;
	}
	if(defined $f_mk_root_seq){
		warn "\nThe generating method for a root seqequence is specified but will be ignored!";
		$f_mk_root_seq=undef;
	}
}

if(@bg_mmtx_fn){
	if(@bg_mmtx_fn==1){
		for(my $i=1;$i<$brsubsets_num;$i++){
			push @bg_mmtx_fn,$bg_mmtx_fn[0];
		}
	}elsif(@bg_mmtx_fn!=$brsubsets_num){
		die "\nThe number of background mutation matrices must be equal to the number of branch subsets: $brsubsets_num!";
	}
}
if(@tg_mmtx_fn){
	if(@tg_mmtx_fn==1){
		for(my $i=1;$i<$brsubsets_num;$i++){
			push @tg_mmtx_fn,$tg_mmtx_fn[0];
		}
	}elsif(@tg_mmtx_fn!=$brsubsets_num){
		die "\nThe number of foreground mutation matrices must be equal to the number of branch subsets: $brsubsets_num!";
	}
}
if(@subset_folders){
	die "\nThe number of 'SubsetFolder' records must be equal to the number of branch subsets: $brsubsets_num!" unless $brsubsets_num==@subset_folders;
	for(my $i=0;$i<$brsubsets_num;$i++){
		$bg_mmtx_fn[$i]=$subset_folders[$i].$bg_mmtx_fn[$i] if @bg_mmtx_fn;
		$tg_mmtx_fn[$i]=$subset_folders[$i].$tg_mmtx_fn[$i] if @tg_mmtx_fn;
	}
}

my $tree=Bio::Phylo::IO->parse(
    -file => $xpar_fn,
    -format => 'adjacency');
#!!!The return value of Adjacency parser not corresponded to interface of the Bio::Phylo::IO->parse
$tree=$tree->[0]->first;
my $rh_terminals;
my $rh_xparr_term_lines;
if($f_ignore_terminals){
	$rh_xparr_term_lines={};
	$rh_terminals={};
	foreach my $node(@{$tree->get_terminals}){
		my $name=$node->get_name;
		$rh_terminals->{$name}=1;
	}
}

#my $str=Bio::Phylo::IO->unparse(-phylo => $tree,-format => 'newick');
#print $str;

#substitution counts
my %bg_nsubst;
my %tg_nsubst;
 
#Original distribution of substitutions on branches
my %tg_subst_map;
my %bg_subst_map;
my %syn_subst_map; #synonimous substitutions
my %phenotype_map;
my $f_phenotypes=0;
my %bg_site2idx; #converts background site coordinate into array index
my @bg_idx2site;
my %tg_site2idx; #converts target site coordinate into array index
my @tg_idx2site;

#New distribution of substitutions on branches
my $hr_smpl_bg_subst;
my $hr_smpl_tg_subst;
my $hr_smpl_bg_subst_map;
my $hr_smpl_tg_subst_map;

sub get_site_profile{
	my ($ra_profile,$alpha_size,$site_idx,$anc_allele_idx)=@_;
	die "\nThe site with index $site_idx is constant or index is wrong!" unless defined $ra_profile->[$site_idx];
	my @tmp=(0) x $alpha_size;
	unless(defined $anc_allele_idx){
		foreach my $state(keys %{$ra_profile->[$site_idx]}){
			$tmp[$state]=$ra_profile->[$site_idx]->{$state};
		}
	}else{
		if(defined $ra_profile->[$site_idx]->{$anc_allele_idx}){
			foreach my $state(keys %{$ra_profile->[$site_idx]->{$anc_allele_idx}}){
				$tmp[$state]=$ra_profile->[$site_idx]->{$anc_allele_idx}->{$state};
			}
		}else{
			die "\nIn the site with index $site_idx no ancestral allele $anc_allele_idx has been accounted!";
		}
	}
	return @tmp;
}

sub init_subst_distribution{
	my ($nsubsets,$ra_mmtx_fn,$ra_branch_fn,$ra_site_fn,$rh_map,$rh_site)=@_;
	%{$rh_site}=();
	%{$rh_map}=();
	for(my $k=0;$k<$nsubsets;$k++){
		my $mmtx_fn=$ra_mmtx_fn->[$k];
		my $branch_fn=$ra_branch_fn->[$k];
		my $site_fn=$ra_site_fn->[$k];
		die "\nError init_mutation_distribution(): Undefined parameter(s)!" unless defined($mmtx_fn)&&defined($branch_fn)&&defined($site_fn);
		
		my @idx2site;
		my @idx2branch;
		my $i=0;
		open INPF, "<$branch_fn" or die "\nUnable to open input file: $branch_fn!";
		while(<INPF>){
			chomp;
			if($_ ne ""){
				my @line=split "\t";
				$idx2branch[$i++]=$line[0];
			};
		};
		close INPF;
		$i=0;
		open INPF, "<$site_fn" or die "\nUnable to open input file: $site_fn!";
		while(<INPF>){
			chomp;
			if($_ ne ""){
				my @line=split "\t";
				$rh_site->{$line[0]}=0 unless defined $rh_site->{$line[0]};
				$idx2site[$i++]=$line[0];
			};
		};
		close INPF;
		open INPF, "<$mmtx_fn" or die "\nUnable to open input file: $mmtx_fn!";
		$i=0;
		#DEBUG
		#while(<INPF>){
		#	chomp;
		#	my $bname=$idx2branch[$i++];
		#	next if $_ eq "";
		#	$rh_map->{$bname}={};
		#	while(/(1+)/g){
		#		my $idx=pos;
		#		my $n=length $1;
		#		for(my $j=$idx-$n;$j<$idx;$j++){
		#			my $spos=$idx2site[$j];
		#			$rh_site->{$spos}++;
		#			$rh_map->{$bname}->{$spos}=1;
		#		}
		#	}
		#}
		#END DEBUG			
		while(<INPF>){
			$_=~s/\s+$//;
			chomp;
			my $bname=$idx2branch[$i++];
			next if $_ eq "";
			$rh_map->{$bname}={};
			my @sites=split(/;/, $_);
			foreach my $j(@sites){
				my $spos=$idx2site[$j-1];
				$rh_site->{$spos}++;
				$rh_map->{$bname}->{$spos}=1;
			}
		}
		die "\nThe number of lines in the $mmtx_fn\n\tisn't equal to the number of branches in the $branch_fn!" unless $i==@idx2branch;
		close INPF;
	}
}

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
	my $ra_alphabet=$args{'-alphabet'};
	my $alpha_size;
	if(defined $ra_alphabet){
		$alpha_size=@{$ra_alphabet};
	}
	if(defined $args{'-sequence_str'}){
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
	}elsif(defined($args{'-profile'})||defined('-profile_fn')){
		my $nsites=$args{'-length'};
		die "\nUndefined number of sites to generate the sequence!" unless defined $nsites;
		my $profile=$args{'-profile'};
		my $alnprof_fn=$args{'-profile_fn'};
		my @alnprof;
		my $f_mk_seq_from;
		if(defined $alnprof_fn){
			my @idx2site;
			$alpha_size=$args{'-alpha_size'} unless defined $alpha_size;
			die "\nThe alphabet size is requred, use '-alphabet_size' parameter to specify!" unless defined $alpha_size;
			read_align_profile($alnprof_fn,\@alnprof,\@idx2site);
			$f_mk_seq_from="align";
			###DEBUG
			#if(defined $args{'-rh_site2index'}){
			#	my $rh_site2idx=$args{'-rh_site2index'};
			#	for(my $i=0;$i<@idx2site;$i++){
			#		next unless defined $idx2site[$i];
			#		die "\nIn file $alnprof_fn wrong site indices!" unless $rh_site2idx->{$idx2site[$i]}==$i;
			#	}
			#}
			#####
		}elsif(!($args{'-seq_from_align'})){
			@{$ra_out_seq}=$profile->get_root_seq;
			return 1;
		}
		my $i=0;
		while($i<$nsites){
			my @prof;
			if(defined $profile){
				@prof=$profile->get_site_profile($i);
			}else{
				@prof=get_site_profile(\@alnprof,$alpha_size,$i);
			}
			$i++;
			my $k=sample_base_from_profile(\@prof);
			push @{$ra_out_seq},$k;
		}
###DEBUG
#print STDERR "\n";
#for(my $i=0;$i<$nsites;$i++){
#	print STDERR $alphabet[$ra_out_seq->[$i]];
#}
#exit;
######
		return 1;
	}
	return 0;
}

sub init_alleles{
	my ($tree,$ra_root_seq,$rh_smpl_subst_map,$rh_site2idx,%args)=@_;
	my $profile;
	my $alpha_size;
	my @amfreq;
	if(defined($args{'-profile'})||defined($args{'-profile_fn'})){
		$profile=$args{'-profile'};
		my $amfreq_fn=$args{'-profile_fn'};
		if(defined $amfreq_fn){
			my @idx2site;
			$alpha_size=$args{'-alpha_size'};
			die "\nThe alphabet size is requred, use '-alphabet_size' parameter to specify!" unless defined $alpha_size;
			read_allele_mutation_matrix($amfreq_fn,\@amfreq,\@idx2site);
			###DEBUG
			#for(my $i=0;$i<@idx2site;$i++){
			#	next unless defined $idx2site[$i];
			#	die "\nIn file $amfreq_fn wrong site indices!" unless $rh_site2idx->{$idx2site[$i]}==$i;
			#}
			#####
		}
	}else{
		die "\nUnable to initilize mutational model!";
	}
	my %align;
	$align{$tree->get_root()->get_name}=$ra_root_seq;
	$tree->visit_breadth_first(
		-in   => sub {
			my $node=shift;
			if(!$node->is_root){
				my $name=$node->get_name;
				my $pnode=$node->get_parent;
				my $pname=$pnode->get_name;
				my $ra_seq=[@{$align{$pname}}];
				if(defined $rh_smpl_subst_map->{$name}){
					foreach my $site(keys %{$rh_smpl_subst_map->{$name}}){
						my $site_idx=$rh_site2idx->{$site};
						die "\nUnable to find index for the site $site" unless defined $site_idx;
						my @prof;
						if(defined $profile){
							@prof=$profile->get_site_profile($site_idx,$ra_seq->[$site_idx]);
						}else{
							@prof=get_site_profile(\@amfreq,$alpha_size,$site_idx,$ra_seq->[$site_idx]);
						}
						my $k=sample_base_from_profile(\@prof);
						$rh_smpl_subst_map->{$name}->{$site}=[($ra_seq->[$site_idx],$k)];
						$ra_seq->[$site_idx]=$k;
					}
				}
				$align{$name}=$ra_seq;
			}
		},
		-post => sub {
			my $node=shift;
			my $name=$node->get_name;
			$align{$name}=undef;
			delete $align{$name};
		}
	);
#!!!DEBUG: print alignment
#foreach my $bname(sort keys %{$rh_align}){
#	print ">$bname\n";
#	for(my $i=0;$i<@{$ra_root_seq};$i++){
#		if(defined $align{$bname}->[$i]){
#			print $alphabet[$align{$bname}->[$i]];
#		}else{
#			print " ";
#		}
#	}
#	print "\n";
#}
#####
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
	my $rh_subst_map=shift;
	my $name=shift;
	my $pname=shift;
	my $ra_alphabet=shift;
	my $ra_sites=shift;
	my @sites=sort {$a<=>$b} @{$ra_sites};
	my $str="";
	my $rh_branch_subst_map=$rh_subst_map->{$name};
	for(my $i=0;$i<@sites;$i++){
		my $site=$sites[$i];
		die "\nError unparse_subst2(): unable to find site $site on the branch $name!" unless defined $rh_branch_subst_map->{$site};
		my @baseids=@{$rh_branch_subst_map->{$site}};
		die "\nError unparse_subst(): for the node $name the mutation in the site $site is declared but not found:!" if $baseids[0]==$baseids[1];
		my @bases=($ra_alphabet->[$baseids[0]],$ra_alphabet->[$baseids[1]]);
		$str.=$bases[0].$site.$bases[1];
		$str.=";" unless $i==$#sites;
	}
	return $str;
}	
	
#begin script
my $bgr_profile;
my $fgr_profile;
my $rh_bgr_align;
my $rh_fgr_align;
my $rh_bg_site2idx;
my $rh_tg_site2idx;

$xparr_header_str=parse_xparr($xpar_fn,"-f_treat_alleles" => $f_treat_allleles,"-bgr_treat_profile" => $bgr_treat_prof,"-fgr_treat_profile" => $fgr_treat_prof,
	"-allele_indices" => $rh_allele2idx,"-f_ignore_terminals" => $f_ignore_terminals, "-term_branches" => $rh_terminals,
	"-orh_bgr_site_nsubst" => \%bg_nsubst,"-orh_fgr_site_nsubst" => \%tg_nsubst,
	"-orh_bgr_subst_map" => \%bg_subst_map, "-orh_fgr_subst_map" => \%tg_subst_map,
	"-orh_bgr_site2idx" => \%bg_site2idx, "-orh_fgr_site2idx" => \%tg_site2idx,
	"-ora_bgr_idx2site" => \@bg_idx2site, "-ora_fgr_idx2site" => \@tg_idx2site,
	"-orh_column_str=3" => \%syn_subst_map, "-orh_column_str=6" => \%phenotype_map,
	"-orh_xparr_terms" => $rh_xparr_term_lines
);
$f_phenotypes=($xparr_header_str=~m/phenotypes:/);

my $f_bgr_eq_fgr=0;
if($f_intragene&&$f_mtx_mode==3){
	if(@bg_mmtx_fn==0){
		$f_bgr_eq_fgr=1;
		@bg_mmtx_fn=@tg_mmtx_fn;
		@bg_idx2site_fn=@tg_idx2site_fn;
		@bg_idx2branch_fn=@tg_idx2branch_fn;
	}elsif(@tg_mmtx_fn==0){
		$f_bgr_eq_fgr=1;
		@tg_mmtx_fn=@bg_mmtx_fn;
		@tg_idx2site_fn=@bg_idx2site_fn;
		@tg_idx2branch_fn=@bg_idx2branch_fn;
	}
}
if($f_mtx_mode!=2){
	die "Error: one or more input files for background are undefined!" unless @bg_mmtx_fn&&@bg_idx2branch_fn&&@bg_idx2site_fn;
	$hr_smpl_bg_subst_map={};
	$hr_smpl_bg_subst={};
	init_subst_distribution($brsubsets_num,\@bg_mmtx_fn,\@bg_idx2branch_fn,\@bg_idx2site_fn,$hr_smpl_bg_subst_map,$hr_smpl_bg_subst);
	###DEBUG
	foreach my $site(keys %bg_nsubst){
		die "\nThe site $site is not found in any of files specified in 'BgrIdx2Site' options!" unless defined $hr_smpl_bg_subst->{$site};
	}
	###
	if(defined($alphabet_str)){
		my $ra_bgr_root_seq=[];
		$bgr_profile=SiteModel::EmpiricalProfile->new($tree,\%bg_subst_map,scalar(@alphabet),\%bg_site2idx,1);
		if(defined $bgr_root_seq){
			init_sequence('-result' => $ra_bgr_root_seq,'-sequence_str' =>$bgr_root_seq,'-alphabet'=>$ra_alphabet);
		}else{
			my $t=1;
			$t=0 if defined($f_mk_root_seq)&&($f_mk_root_seq eq "xparr");
			init_sequence('-result' => $ra_bgr_root_seq,'-profile' => $bgr_profile,'-length' => scalar(@bg_idx2site), '-seq_from_align' => $t);
		}
		init_alleles($tree,$ra_bgr_root_seq,$hr_smpl_bg_subst_map,\%bg_site2idx,'-profile' => $bgr_profile);
	}
}
if($f_mtx_mode!=1){
	die "Error: one or more input files for foreground are undefined!" unless @tg_mmtx_fn&&@tg_idx2branch_fn&&@tg_idx2site_fn;
	unless($f_bgr_eq_fgr){
		$hr_smpl_tg_subst_map={};
		$hr_smpl_tg_subst={};
		init_subst_distribution($brsubsets_num,\@tg_mmtx_fn,\@tg_idx2branch_fn,\@tg_idx2site_fn,$hr_smpl_tg_subst_map,$hr_smpl_tg_subst);
		###DEBUG
		foreach my $site(keys %tg_nsubst){
			die "\nThe site $site is not found in any of files specified in 'FgrIdx2Site' options!" unless defined $hr_smpl_tg_subst->{$site};
		}
		###
		if(defined($alphabet_str)){
			my $ra_fgr_root_seq=[];
			$fgr_profile=SiteModel::EmpiricalProfile->new($tree,\%tg_subst_map,scalar(@alphabet),\%tg_site2idx,1);
			if(defined $fgr_root_seq){
				init_sequence('-result' => $ra_fgr_root_seq,'-sequence_str' =>$fgr_root_seq,'-alphabet'=>$ra_alphabet);
			}else{
				my $t=1;
				$t=0 if defined($f_mk_root_seq)&&($f_mk_root_seq eq "xparr");
				init_sequence('-result' => $ra_fgr_root_seq,'-profile' => $fgr_profile,'-length' => scalar(@tg_idx2site), '-seq_from_align' => $t);
			}
			init_alleles($tree,$ra_fgr_root_seq,$hr_smpl_tg_subst_map,\%tg_site2idx,'-profile' => $fgr_profile);
		}
	}else{
		$hr_smpl_tg_subst_map=$hr_smpl_bg_subst_map;
		$hr_smpl_tg_subst=$hr_smpl_bg_subst;
	}
}

if($f_intragene&&(!$f_bgr_eq_fgr)){
	#check that permutated subsets of sites don't overlap
	foreach my $site(keys %{$hr_smpl_tg_subst}){
		if(exists $hr_smpl_bg_subst->{$site}){
			die "\nError: Independent permutations of both background and foreground overlapped subsets of sites are impossible!";
		}
	}
}

#print XPARR
print $xparr_header_str;;
foreach my $node($tree->get_nodes){
	my $name=$node->get_name;
	if($f_ignore_terminals&&$node->is_terminal){
		print "\n".$rh_xparr_term_lines->{$name};
		next;
	}
	print "\n$name";
	if(!$node->is_root){
		my $pname=$node->get_parent->get_name;
		my $length=$node->get_branch_length;
		print "\t$pname\t$length\t$syn_subst_map{$name}\t";
		my @sites1;
		my @sites2;
		my @tmp;
		my $str="";
		if(defined $tg_subst_map{$name}){
			my $rh_smpl_subst_map;
			if(defined $hr_smpl_tg_subst_map){
				@sites2=keys %{$hr_smpl_tg_subst_map->{$name}};
				$rh_smpl_subst_map=$hr_smpl_tg_subst_map;
			}else{
				@tmp=keys %{$tg_subst_map{$name}};
			}
			if(!$f_intragene){
				foreach my $site(@tmp){
					push @sites1,$tg_subst_map{$name}->{$site};
				}
			}elsif(!defined($hr_smpl_tg_subst_map)){
				foreach my $site(@tmp){
					push @sites1,$tg_subst_map{$name}->{$site} unless defined $bg_nsubst{$site};
				}
				@tmp=keys %{$hr_smpl_bg_subst_map->{$name}};
				for(my $i=0;$i<@tmp;$i++){
					push @sites2,$tmp[$i] if defined $tg_nsubst{$tmp[$i]};
				}
				$rh_smpl_subst_map=$hr_smpl_bg_subst_map;
			}
			die "\nError: no mutations in foreground on the branch: $name!" unless @sites1+@sites2;
			$str.=unparse_subst1(\@sites1,$ra_alphabet) if @sites1;
			if(@sites2){
				$str.=";" if @sites1;
				if(defined $alphabet_str){
					$str.=unparse_subst2($rh_smpl_subst_map,$name,$pname,$ra_alphabet,\@sites2);
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
		if(defined $bg_subst_map{$name}){
			my $rh_smpl_subst_map;
			if(defined $hr_smpl_bg_subst_map){
				@sites2=keys %{$hr_smpl_bg_subst_map->{$name}};
				$rh_smpl_subst_map=$hr_smpl_bg_subst_map;
			}else{
				@tmp=keys %{$bg_subst_map{$name}};
			}
			if(!$f_intragene){
				foreach my $site(@tmp){
					push @sites1,$bg_subst_map{$name}->{$site};
				}
			}elsif(!defined $hr_smpl_bg_subst_map){
				for(my $i=0;$i<@tmp;$i++){
					push @sites1,$tmp[$i] unless defined $tg_nsubst{$tmp[$i]};
				}
				@tmp=keys %{$hr_smpl_tg_subst_map->{$name}};
				for(my $i=0;$i<@tmp;$i++){
					push @sites2,$tmp[$i] if(defined $bg_nsubst{$tmp[$i]});
				}
				$rh_smpl_subst_map=$hr_smpl_tg_subst_map;
			}
			die "\nError: no mutations in foreground on the branch: $name!" unless @sites1+@sites2;
			$str.=unparse_subst1(\@sites1,$ra_alphabet) if @sites1;
			if(@sites2){
				$str.=";" if @sites1;
				if(defined $alphabet_str){
					$str.=unparse_subst2($rh_smpl_subst_map,$name,$pname,$ra_alphabet,\@sites2);
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
			