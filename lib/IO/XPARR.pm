package IO::XPARR;
#This module provides class for representing a phylogenetic tree with mutations on branches (XPARR) in the binary matrix form. 
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
# Symbols to autoexport (:DEFAULT tag)
@EXPORT = qw(parse_xparr parse_subst_abbr); 
@EXPORT_OK = qw(parse_xparr);

use Class::Struct;

struct SubstInfo => {
	site => '$',
	weight => '$',
	bases => '@'
};

sub parse_subst_abbr{
	my $str=shift;
	my $ra_bases=shift;
	my $rh_allele2idx=shift;
	my $spos;
	if($str=~/([A-Za-z-])(\d+)([[A-Za-z-])/){
		$spos=$2;
		if(defined $ra_bases){
			@{$ra_bases}=();
			my ($a1,$a2)=(uc $1,uc $3);
			if(defined $rh_allele2idx){
				die "\nError: The alleles could not be omitted if the alphabet string is specified!" unless defined($1)&&defined($3);
				@{$ra_bases}=($rh_allele2idx->{$a1},$rh_allele2idx->{$a2});
			}else{
				@{$ra_bases}=($a1,$a2);
			}
		}
	}else{
		$str=~s/\D//g;
		$spos=$str;
	}
	return $spos;
}

sub parse_xparr{
	my $xpar_fn=shift;
	my %args=@_;
	#Parameters
	#input:
	my $rh_terminals=$args{'-term_branches'};
	my $rh_allele2idx=$args{'-allele_indices'};
	my $f_treat_allleles=0;
	$f_treat_allleles=$args{'-f_treat_alleles'} if defined $f_treat_allleles;
	my $f_ignore_terminals=0;
	$f_ignore_terminals=$args{'-f_ignore_terminals'} if defined $f_ignore_terminals;
	if($f_ignore_terminals){
		die "\nError parse_xparr(): The map of terminal branch names is required! Please specify '-term_branches'" unless defined $rh_terminals;
	}
	my $fgr_treat_prof=$args{'-fgr_treat_profile'};
	my $bgr_treat_prof=$args{'-bgr_treat_profile'};
	#output:
	my $rh_tg_subst_branches=$args{'-orh_fgr_subst2branch'}; #get branches carrying substitutions in each site
	%{$rh_tg_subst_branches}=() if defined $rh_tg_subst_branches;
	my $rh_bg_subst_branches=$args{'-orh_bgr_subst2branch'}; #get branches carrying substitutions in each site
	%{$rh_bg_subst_branches}=() if defined $rh_bg_subst_branches;
	my $ra_pheno_labels=$args{'-ora_pheno_labels'}; #get names of phenotypes
	@{$ra_pheno_labels}=() if defined $ra_pheno_labels;
	my $rh_terms=$args{'-orh_xparr_terms'}; #get lines from xparr file corresponding to terminal branches
	if(defined $rh_terms){
		%{$rh_terms}=();
		die "\nError parse_xparr(): The map of terminal branch names is required! Please specify '-term_branches'" unless defined $rh_terminals;
	}
	my $rh_tg_nsubst=$args{'-orh_fgr_site_nsubst'}; #get numbers of substitutions in sites
	%{$rh_tg_nsubst}=() if defined $rh_tg_nsubst;
	my $rh_bg_nsubst=$args{'-orh_bgr_site_nsubst'}; #get numbers of substitutions in sites
	%{$rh_bg_nsubst}=() if defined $rh_bg_nsubst;
	my $rh_tg_subst_map=$args{'-orh_fgr_subst_map'}; #get substitutions on the tree branches
	%{$rh_tg_subst_map}=() if defined $rh_tg_subst_map;
	my $rh_bg_subst_map=$args{'-orh_bgr_subst_map'}; #get substitutions on the tree branches
	%{$rh_bg_subst_map}=() if defined $rh_bg_subst_map;
	my $rh_tg_site2idx=$args{'-orh_fgr_site2idx'}; #get site indices
	%{$rh_tg_site2idx}=() if defined $rh_tg_site2idx;
	my $rh_bg_site2idx=$args{'-orh_bgr_site2idx'}; #get site indices
	%{$rh_bg_site2idx}=() if defined $rh_bg_site2idx;
	my $ra_tg_idx2site=$args{'-ora_fgr_idx2site'}; #get sites
	@{$ra_tg_idx2site}=() if defined $ra_tg_idx2site;
	my $ra_bg_idx2site=$args{'-ora_bgr_idx2site'}; #get sites
	@{$ra_bg_idx2site}=() if defined $ra_bg_idx2site;
	my $rh_syn_counts=$args{'-orh_syn_counts'}; #get numbers of synonimous mutations on branches
	%{$rh_syn_counts}=() if defined $rh_syn_counts;
	my $rh_tg_nsyn_counts=$args{'-orh_fgr_nsyn_counts'}; #get numbers of synonimous mutations on branches
	%{$rh_tg_nsyn_counts}=() if defined $rh_tg_nsyn_counts;
	my $rh_bg_nsyn_counts=$args{'-orh_bgr_nsyn_counts'}; #get numbers of synonimous mutations on branches
	%{$rh_bg_nsyn_counts}=() if defined $rh_bg_nsyn_counts;
	my $rh_phenotypes=$args{'-orh_phenotypes'}; #get phenotype values on branches
	%{$rh_phenotypes}=() if defined $rh_phenotypes;
	##############################
	my %column_strings;
	my $ncols=0;
	foreach my $key(keys %args){
		if($key=~/^-orh_column_str=(\d+)/){
			my $ci=$1;
			die "\nThe first two columns in the xparr file could not be requested!" if $ci<2; 
			%{$args{$key}}=();
			$ncols++;
			$column_strings{$ci}=$args{$key}; #get the string corresponding to a column with a specified number
		}
	}
	#############################
	my $npheno;
	my @pheno_labels;
	my %tg_nsubst;
	my %bg_nsubst;
	my %tg_subst_map;
	my %bg_subst_map;
	my %syn_counts;
	my %tg_nsyn_counts;
	my %bg_nsyn_counts;
	my %phenotypes;
	my %bg_site2idx;
	my @bg_idx2site;
	my %tg_site2idx;
	my @tg_idx2site;
	open INPF, "<$xpar_fn" or die "\nUnable to open input file $xpar_fn!";
	$_=<INPF>; #skip header line
	chomp;
	s/\s+$//;
	my $header_str=$_;
	my @line=split "\t";
	if(defined $ra_pheno_labels){
		if($line[6]=~/^phenotypes:/){
			@pheno_labels=split ",", $';
			$npheno=@pheno_labels;
		}else{
			"\nNo phenotypes defined on branches of a tree in the file: $xpar_fn!";
		}
		@{$ra_pheno_labels}=@pheno_labels;
	}
	die "\nHeader line was omitted in file: $xpar_fn!" unless /^child\tparent\tlength/;
	$_=<INPF>; #skip root node statement
	die "\nThe root node statement required in file: $xpar_fn!" unless /^[\w.]+\s*$/;
	while(<INPF>){
		chomp;
		s/\r+$//;
		@line=split "\t",$_,-1;
		if(defined($rh_terms)&&defined($rh_terminals->{$line[0]})){
			$rh_terms->{$line[0]}=$_;
		}
		next if($f_ignore_terminals && defined($rh_terminals->{$line[0]}));
		if($ncols){
			for(my $i=2;$i<@line;$i++){
				if(defined $column_strings{$i}){
					$column_strings{$i}->{$line[0]}=$line[$i] if $line[$i] ne "";
				}
			}
		}
		my $n=0;
		if($line[4] ne ""){
			my @sites=split(/;/, $line[4]);
			$n=@sites;
			$tg_subst_map{$line[0]}={};
			foreach my $str(@sites){
				my @bases;
				my $spos=parse_subst_abbr($str,\@bases,$rh_allele2idx);
				$tg_nsubst{$spos}++;
				my $si=SubstInfo->new();
				$si->site($spos);
				$si->bases([@bases]) if @bases;
				$tg_subst_map{$line[0]}->{$spos}=$si;
				if(defined $rh_tg_subst_branches){
					$rh_tg_subst_branches->{$spos}=[] unless defined $rh_tg_subst_branches->{$spos};
					push @{$rh_tg_subst_branches->{$spos}},$line[0];
				}
			}
		}
		$tg_nsyn_counts{$line[0]}=$n;
		$n=0;
		if($line[5] ne ""){
			my @sites=split(/;/, $line[5]);
			$bg_subst_map{$line[0]}={};
			foreach my $str(@sites){
				my @bases;
				my $spos=parse_subst_abbr($str,\@bases,$rh_allele2idx);
				$bg_nsubst{$spos}++;
				my $si=SubstInfo->new();
				$si->site($spos);
				$si->bases([@bases]) if @bases;
				$bg_subst_map{$line[0]}->{$spos}=$si;
				if(defined $rh_bg_subst_branches){
					$rh_bg_subst_branches->{$spos}=[] unless defined $rh_bg_subst_branches->{$spos};
					push @{$rh_bg_subst_branches->{$spos}},$line[0];
				}
			}
		}
		$bg_nsyn_counts{$line[0]}=$n;
		if($npheno){
			if($line[6]=~/^[10]+$/){
				my @pheno=split "", $line[6];
				if(@pheno==$npheno){
					$phenotypes{$line[0]}=[@pheno];
				}else{
					die "\nNumber of labels doesn't equal to the number of phenotypes on the branch $line[0]!";
				}
			}else{
				die "\nUndefined phenotypes on the branch: $line[0]!";
			}
		}
		$n=0;
		if($line[3] ne ""){
			while($line[3]=~m/;/g){$n++};
			$n++;
		}
		$syn_counts{$line[0]}=$n;
	}
	close INPF;
	#create site indices
	{
		my $I=0;
		foreach my $spos(sort {$a <=> $b} keys %bg_nsubst){
			$bg_site2idx{$spos}=$I;
			$bg_idx2site[$I]=$spos;
			$I++;
		}
		$I=0;
		foreach my $spos(sort {$a <=> $b} keys %tg_nsubst){
			$tg_site2idx{$spos}=$I if defined $rh_tg_site2idx;
			$tg_idx2site[$I]=$spos if defined $ra_tg_idx2site;
			$I++;
		}
	}
	if($f_treat_allleles){
		foreach my $bname(keys %bg_subst_map){
			foreach my $spos(keys %{$bg_subst_map{$bname}}){
				my $si=$bg_subst_map{$bname}->{$spos};
				if(defined $bgr_treat_prof){
					die "\nError: Unable to treat alleles in the background mutation on the branch $bname site $spos!" 
						unless $bgr_treat_prof->treat_alleles($si,\%bg_site2idx);
				}else{
					delete $bg_subst_map{$bname}->{$spos} unless defined($si->bases(0))&&defined($si->bases(1));
					$bg_nsubst{$spos}--;
				}
			}
		}
		foreach my $bname(keys %tg_subst_map){
			foreach my $spos(keys %{$tg_subst_map{$bname}}){
				my $si=$tg_subst_map{$bname}->{$spos};
				if(defined $fgr_treat_prof){
					die "\nError: Unable to treat alleles in the foreground mutation on the branch $bname site $spos!" 
						unless $fgr_treat_prof->treat_alleles($si,\%tg_site2idx);
				}else{
					delete $tg_subst_map{$bname}->{$spos} unless defined($si->bases(0))&&defined($si->bases(1));
					$tg_nsubst{$spos}--;
				}
			}
		}
	}
	%{$rh_tg_nsubst}=%tg_nsubst if defined $rh_tg_nsubst;
	%{$rh_bg_nsubst}=%bg_nsubst if defined $rh_bg_nsubst;
	%{$rh_tg_subst_map}=%tg_subst_map if defined $rh_tg_subst_map;
	%{$rh_bg_subst_map}=%bg_subst_map if defined $rh_bg_subst_map;
	%{$rh_syn_counts}=%syn_counts if defined $rh_syn_counts;
	%{$rh_tg_nsyn_counts}=%tg_nsyn_counts if defined $rh_tg_nsyn_counts;
	%{$rh_bg_nsyn_counts}=%bg_nsyn_counts if defined $rh_bg_nsyn_counts;
	%{$rh_phenotypes}=%phenotypes if defined $rh_phenotypes;
	%{$rh_tg_site2idx}=%tg_site2idx if defined $rh_tg_site2idx;
	@{$ra_tg_idx2site}=@tg_idx2site if defined $ra_tg_idx2site;
	%{$rh_bg_site2idx}=%bg_site2idx if defined $rh_bg_site2idx;
	@{$ra_bg_idx2site}=@bg_idx2site if defined $ra_bg_idx2site;
	return $header_str;
}

1;