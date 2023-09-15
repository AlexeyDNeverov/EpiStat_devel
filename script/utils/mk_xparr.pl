#!/usr/bin/env perl
#Converts tree for "Mitochondria" project into XPARR
#Options:
#	-t <TREE> - a file with tree in newick format
#	[-m <0<MinMuts>] - Minimal number of mutations in a site
#	[-p <'sites'|'mutations'>] - Print mode: sites or mutations
#		Default='sites'
#	[-g <'all'|'internal'|'no'>] - Use gaps. Allows to ignore gaps on the tree's leafs
#		Default='all' - Do not ignore terminal gaps on leafs
#	[-G <max_gaps_frac>] - A filter of alignment sites: ignore sites those fraction of internal gaps is higher than <max_gaps_frac>=[0,1]
#		default=1
#	[-s <"n(ucleotides)?"|"a(mino_acids)?">] -  Type of aligned sequences.
#		The script tries to guess sequence type if it has not been directly specified
#	[-T <0|1>] - Translate nucleotide alignments
#		Default=1
#	[-z] - Remove internal nodes with zero mutations (all mutations are considered)
#Params:
#	<ALIGN> - a file with alignment in FASTA
#	[<ALIGN>] - a file with alignment in FASTA. If it was accounted then the first file is interpreted as background and the second as foreground
use lib "$ENV{EPISTAT_LIB}";
#use warnings;
#no warnings 'recursion';
use strict;
use Bio::Phylo::IO;
use Getopt::Std;
use DnaUtilities::FASTASimple;
use Class::Struct;
use Bio::Tools::CodonTable;

my %args;
if(!getopts('t:m:p:g:G:s:T:z',\%args)){
	die "\nError in option string!";
}
my $f_ignore_gaps=0;
my @gaps_fracts;
my @seqs_counts;
my $MaxGapsFrac;
my $f_nucleotides;
my $f_translate=0;
my $codonTable;
if(defined($args{g})){
	if($args{g}=~/^internal/i){
		$f_ignore_gaps=1;
	}elsif($args{g}=~/^no/i){
		$f_ignore_gaps=2;
	}elsif(!$args{g}=~/^all/i){
		die "\nWrong value for the option -g accounted. Valid values are 'internal' or 'all'!";
	}
}
my $f_intragene=1;
$f_intragene=0 if defined $ARGV[1];
if(defined $args{G}){
	$MaxGapsFrac=$args{G};
	die "\nThe maximal fractions of gaps accounted in the '-G' parameter should be an unsigned floating point!" unless $MaxGapsFrac=~m/\d+\.\d+/;
	die "\nThe maximal fractions of gaps accounted in the '-G' parameter should be within the range [0,1]!" unless $MaxGapsFrac>=0&&$MaxGapsFrac<=1;
	push @gaps_fracts,[];
	push @seqs_counts,[];
	unless($f_intragene){
		push @gaps_fracts,[];
		push @seqs_counts,[];
	}
}
my $MinMuts=1;
$MinMuts=int $args{m} if defined $args{m};
die "\nUnpropper minimal mutations' number - not UINT!" unless $MinMuts>0;
my $f_print_mutations=0;
if($args{p}){
	if($args{p}=~m/sites/i){
		$f_print_mutations=0;
	}elsif($args{p}=~m/mutations/i){
		$f_print_mutations=1;
	}else{
		die "\nUnknown print mode: $args{p}!";
	}
}
my $tree;
if($args{t}){
	$tree = Bio::Phylo::IO->parse(
		'-file' => $args{t},
		'-format' => 'newick',
		)->first;
}else{
	die "\nThe file with phylogeny is not specified! Use '-t' option!";
}
if($args{s}){
	if($args{s}=~/^n(ucleotides)?/i){
		$f_nucleotides=1;
	}elsif($args{s}=~/^a(mino_acids)?/i){
		$f_nucleotides=0;
	}else{
		die "\nUncnown sequence type was specified: -s ".$args{s};
	}
}else{
	#guess the sequence type
	open INPF, "<$ARGV[0]" or die "\nUnable to open input file ".$ARGV[0];
	my $str;
	while(<INPF>){
		chomp;
		if(/^>/){
			last if defined $str;
		}elsif(/\S+/){
			s/\s+//;
			$str.=uc $_;
		}
	}
	close INPF;
	my %counts;
	if(defined $str){
		my @seq=split //,$str;
		my $N=0;
		foreach my $sym(@seq){
			if($sym ne '-'){
				$counts{$sym}++;
				$N++;
			}
		}
		if(($counts{A}+$counts{C}+$counts{G}+$counts{T}+$counts{U})/$N>0.9){
			$f_nucleotides=1;
		}else{
			my @alphabet=split //,"ARNDCQEGHILKMFPSTWYV";
			my $n=0;
			foreach my $a(@alphabet){
				$n+=$counts{$a};
			}
			if($n/$N>0.9){
				$f_nucleotides=0;
			}else{
				die "\nUnable to guess sequence type using the first sequence from $ARGV[0]!";
			}
		}
	}
}
$f_translate=1 if $f_nucleotides;
if($args{T}){
	if($args{T} eq "1"){
		die "\nUnable to translate protein sequences!" unless $f_nucleotides;
	}elsif($args{T} eq "0"){
		$f_translate=0;
	}else{
		die "\nWrong parameter for the option -T ".$args{T}." (only 0 or 1 are allowed).";
	}
}
$codonTable=Bio::Tools::CodonTable->new() if $f_translate;

my @visited_nodes;
my $buff_size=1000;

struct SubstInfo => {
	site => '$',
	bases => '@'
};

sub seq_internals{
	my ($ras)=@_;
	my ($from,$to)=(0,@{$ras}-1);
	while($from<@{$ras}){
		if($ras->[$from] =~ /-/){
			$from++;
		}else{
			last;
		}
	}
	while($to>=0){
		if($ras->[$to] =~ /-/){
			$to--;
		}else{
			last;
		}
	}
	return ($from,$to);
}

sub two_seq_diff{
	my($rh_align,$pnode,$node,$f_ignore_gaps,$f_print_mutations,$myCodonTable,$ra_syn)=@_;
	my $pname=$pnode->get_name;
	my $name=$node->get_name;
	my $ras1=$rh_align->{$pname};
	my $ras2=$rh_align->{$name};
	@{$ra_syn}=() if(defined $ra_syn);
	die "\nError two_seq_diff(): Input arrays have unequal sizes!" unless @{$ras1}==@{$ras2};
	my ($from,$to)=(0,@{$ras1}-1);
	if($f_ignore_gaps==1){
		if($node->is_terminal){
			($from,$to)=seq_internals($ras2);
		}
	}
	my @diff;
	for(my $i=$from;$i<=$to;$i++){
		if($ras1->[$i] ne $ras2->[$i]){
			next if ($f_ignore_gaps==2)&&(($ras2->[$i] =~ /-/)||($ras1->[$i] =~ /-/));
			my $si=SubstInfo->new();
			$si->site($i+1);
			my ($a1,$a2)=($ras1->[$i],$ras2->[$i]);
			if(defined $myCodonTable){
				warn "The codon ".($i+1)."(".$ras1->[$i].") in the sequence $pname possibly is a terminator" if $myCodonTable->is_ter_codon($ras1->[$i]);
				warn "The codon ".($i+1)."(".$ras2->[$i].") in the sequence $name possibly is a terminator" if $myCodonTable->is_ter_codon($ras2->[$i]);
				($a1,$a2)=($myCodonTable->translate($ras1->[$i]),$myCodonTable->translate($ras2->[$i]));
				if($a1 eq $a2){
					push @{$ra_syn},$si if(defined $ra_syn);
				}else{
					$si->bases([($a1,$a2)]) if $f_print_mutations;
					push @diff,$si;
				}
			}else{
				$si->bases([($a1,$a2)]) if $f_print_mutations;
				push @diff,$si;
			}
		}
	}
	return @diff;
}

sub count_gaps{
	my ($rai_seq,$rao_gaps_counts,$rao_seqs_counts)=@_;
	my ($from,$to)=(0,@{$rai_seq}-1);
	for(;$from<@{$rai_seq};$from++){last unless $rai_seq->[$from] =~ /-/};
	die "\nFatal error in count_gaps(): The sequence entierly consisted of gaps!" if $to<$from;
	for(;$to>=0;$to--){last unless $rai_seq->[$to] =~ /-/};
	for(my $i=$from;$i<=$to;$i++){
		$rao_seqs_counts->[$i]++;
		$rao_gaps_counts->[$i]++ if $rai_seq->[$i] =~ /-/;
	}
}

sub print_gaps_stats{
	my ($ra_gaps_counts,$ra_seqs_counts,$alength,$out_file_name)=@_;
	my $fh;
	if(defined $out_file_name){
		$fh=IO::File->new(">$out_file_name");
	}else{
		$fh=*STDOUT;
	}
	print $fh "site\tngaps\tntot";
	for(my $i=0;$i<$alength;$i++){
		print $fh "\n".($i+1)."\t".$ra_gaps_counts->[$i]."\t".$ra_seqs_counts->[$i];
	}
	$fh->close() if(defined $out_file_name);
}

my $ra_align=[{},undef];
my $ra_site_nmuts=[{},undef];
my @alength=(0,0);
my %ndaughters;
if(!$f_intragene){
	$ra_align->[1]={};
	$ra_site_nmuts->[1]={};
}

$tree->visit_breadth_first(
	-in => sub{
		my $node=shift;
		my $name=$node->get_name;
		$ndaughters{$name}=@{$node->get_children} if(!$node->is_terminal);
		if(@visited_nodes==$buff_size){
			my @sids;
			foreach my $nd(@visited_nodes){
				push @sids,$nd->get_name;
			}
			my @seq_names;
			my @seqs;
			my $n=DnaUtilities::FASTASimple::fetchFASTA_by_ids($ARGV[0],\@sids,\@seq_names,\@seqs);
			foreach(my $i=0;$i<$n;$i++){
				$ra_align->[0]->{$seq_names[$i]}=[];
				if($f_translate){
					@{$ra_align->[0]->{$seq_names[$i]}}=unpack '(A3)*',$seqs[$i];
					warn "\nThe length of sequence >$seq_names[$i] from $ARGV[0] is not multiple by tree!" if length($seqs[$i])%3;
				}else{
					@{$ra_align->[0]->{$seq_names[$i]}}=unpack '(A1)*',$seqs[$i];
				}
			}
			if($n!=@sids){
				print STDERR "\nError: Some of requested sequences were not retreived from the file $ARGV[0]:";
				foreach my $name(@sids){
					print STDERR "\n\t$name" unless exists $ra_align->[0]->{$name};
				}
				exit;
			}
			if(!$f_intragene){
				@seq_names=();
				@seqs=();
				$n=DnaUtilities::FASTASimple::fetchFASTA_by_ids($ARGV[1],\@sids,\@seq_names,\@seqs);
				foreach(my $i=0;$i<$n;$i++){
					$ra_align->[1]->{$seq_names[$i]}=[];
					if($f_translate){
						@{$ra_align->[1]->{$seq_names[$i]}}=unpack '(A3)*',$seqs[$i];
						warn "\nThe length of sequence >$seq_names[$i] from $ARGV[1] is not multiple by tree!" if length($seqs[$i])%3;
					}else{
						@{$ra_align->[1]->{$seq_names[$i]}}=unpack '(A1)*',$seqs[$i];
					}
				}
				if($n!=@sids){
					print STDERR "\nError: Some of requested sequences were not retreived from the file $ARGV[1]:";
					foreach my $name(@sids){
						print STDERR "\n\t$name" unless exists $ra_align->[1]->{$name};
					}
					exit;
				}
			}
			foreach my $nd(@visited_nodes){
				if(!$nd->is_root){
					my $pname=$nd->get_parent->get_name;
					my $name=$nd->get_name;
					$nd->set_generic('nsyn' => [(undef,undef)]);
					$nd->set_generic('syn' => [(undef,undef)]);
					my @syn;
					my @nsyn=two_seq_diff($ra_align->[0],$nd->get_parent,$nd,$f_ignore_gaps,$f_print_mutations,$codonTable,\@syn);
					if($MinMuts>1){
						foreach my $site(@nsyn){
							$ra_site_nmuts->[0]->{$site->site}++;
						}
					}
					$nd->get_generic('nsyn')->[0]=[@nsyn];
					if($f_intragene){
						$nd->get_generic('syn')->[0]=[@syn] if defined $codonTable;
					}else{
						my @syn;
						my @nsyn=two_seq_diff($ra_align->[1],$nd->get_parent,$nd,$f_ignore_gaps,$f_print_mutations,$codonTable,\@syn);
						if($MinMuts>1){
							foreach my $site(@nsyn){
								$ra_site_nmuts->[1]->{$site->site}++;
							}
						}
						$nd->get_generic('nsyn')->[1]=[@nsyn];
						$nd->get_generic('syn')->[1]=[@syn] if defined $codonTable;
					}
					$ndaughters{$pname}--;
					if($ndaughters{$pname}==0){
						delete $ra_align->[0]->{$pname};
						delete $ra_align->[1]->{$pname} unless $f_intragene;
						delete $ndaughters{$pname};
					}
					if($nd->is_terminal){
	 					count_gaps($ra_align->[0]->{$name},$gaps_fracts[0],$seqs_counts[0]) if(defined $gaps_fracts[0]);
						$alength[0]=@{$ra_align->[0]->{$name}} if $alength[0]<@{$ra_align->[0]->{$name}};
						delete $ra_align->[0]->{$name};
						unless($f_intragene){
							count_gaps($ra_align->[1]->{$name},$gaps_fracts[1],$seqs_counts[1]) if(defined $gaps_fracts[1]);
							$alength[1]=@{$ra_align->[1]->{$name}} if $alength[1]<@{$ra_align->[1]->{$name}};
							delete $ra_align->[1]->{$name};
						}
					}
				}
			}
			@visited_nodes=();
		}
		push @visited_nodes,$node;
	}
);

if(@visited_nodes){
	my @sids;
	foreach my $nd(@visited_nodes){
		push @sids,$nd->get_name;
	}
	my @seq_names;
	my @seqs;
	my $n=DnaUtilities::FASTASimple::fetchFASTA_by_ids($ARGV[0],\@sids,\@seq_names,\@seqs);
	foreach(my $i=0;$i<$n;$i++){
		$ra_align->[0]->{$seq_names[$i]}=[];
		if($f_translate){
			@{$ra_align->[0]->{$seq_names[$i]}}=unpack '(A3)*',$seqs[$i];
			warn "\nThe length of sequence >$seq_names[$i] from $ARGV[0] is not multiple by tree!" if length($seqs[$i])%3;
		}else{
			@{$ra_align->[0]->{$seq_names[$i]}}=unpack '(A1)*',$seqs[$i];
		}
	}
	if($n!=@sids){
		print STDERR "\nError: Some of requested sequences were not retreived from the file $ARGV[0]:";
		foreach my $name(@sids){
			print STDERR "\n\t$name" unless exists $ra_align->[0]->{$name};
		}
		exit;
	}
	if(!$f_intragene){
		@seq_names=();
		@seqs=();
		$n=DnaUtilities::FASTASimple::fetchFASTA_by_ids($ARGV[1],\@sids,\@seq_names,\@seqs);
		foreach(my $i=0;$i<$n;$i++){
			$ra_align->[1]->{$seq_names[$i]}=[];
			if($f_translate){
				@{$ra_align->[1]->{$seq_names[$i]}}=unpack '(A3)*',$seqs[$i];
				warn "\nThe length of sequence >$seq_names[$i] from $ARGV[1] is not multiple by tree!" if length($seqs[$i])%3;
			}else{
				@{$ra_align->[1]->{$seq_names[$i]}}=unpack '(A1)*',$seqs[$i];
			}
		}
		if($n!=@sids){
			print STDERR "\nError: Some of requested sequences were not retreived from the file $ARGV[1]:";
			foreach my $name(@sids){
				print STDERR "\n\t$name" unless exists $ra_align->[1]->{$name};
			}
			exit;
		}
	}
	foreach my $nd(@visited_nodes){
		if(!$nd->is_root){
			my $pname=$nd->get_parent->get_name;
			my $name=$nd->get_name;
			$nd->set_generic('nsyn' => [(undef,undef)]);
			$nd->set_generic('syn' => [(undef,undef)]);
			my @syn;
			my @nsyn=two_seq_diff($ra_align->[0],$nd->get_parent,$nd,$f_ignore_gaps,$f_print_mutations,$codonTable,\@syn);
			if($MinMuts>1){
				foreach my $site(@nsyn){
					$ra_site_nmuts->[0]->{$site->site}++;
				}
			}
			$nd->get_generic('nsyn')->[0]=[@nsyn];
			if($f_intragene){
				$nd->get_generic('syn')->[0]=[@syn] if defined($codonTable);
			}else{
				my @syn;
				my @nsyn=two_seq_diff($ra_align->[1],$nd->get_parent,$nd,$f_ignore_gaps,$f_print_mutations,$codonTable,\@syn);
				if($MinMuts>1){
					foreach my $site(@nsyn){
						$ra_site_nmuts->[1]->{$site->site}++;
					}
				}
				$nd->get_generic('nsyn')->[1]=[@nsyn];
				$nd->get_generic('syn')->[1]=[@syn] if defined $codonTable;
			}
			$ndaughters{$pname}--;
			if($ndaughters{$pname}==0){
				delete $ra_align->[0]->{$pname};
				delete $ra_align->[1]->{$pname} unless $f_intragene;
				delete $ndaughters{$pname};
			}
			if($nd->is_terminal){
				count_gaps($ra_align->[0]->{$name},$gaps_fracts[0],$seqs_counts[0]) if(defined $gaps_fracts[0]);
				$alength[0]=@{$ra_align->[0]->{$name}} if $alength[0]<@{$ra_align->[0]->{$name}};
				delete $ra_align->[0]->{$name};
				unless($f_intragene){
					count_gaps($ra_align->[1]->{$name},$gaps_fracts[1],$seqs_counts[1]) if(defined $gaps_fracts[1]);
					$alength[1]=@{$ra_align->[1]->{$name}} if $alength[1]<@{$ra_align->[1]->{$name}};
					delete $ra_align->[1]->{$name};
				}
			}
		}
	}
	@visited_nodes=();
}
if(defined $MaxGapsFrac){
	for(my $I=0;$I<2;$I++){
		last if ($I==1)&&($f_intragene);
		my $str=$ARGV[$I];
		$str=~s/\.\S+$//;
		$str.=".gaps.stat";
		print_gaps_stats($gaps_fracts[$I],$seqs_counts[$I],$alength[$I],$str);
		for(my $i=0;$i<$alength[$I];$i++){
			$gaps_fracts[$I]->[$i]=0 unless defined $gaps_fracts[$I]->[$i];
			unless(defined $seqs_counts[$I]->[$i]){
				$seqs_counts[$I]->[$i]=0;
				$gaps_fracts[$I]->[$i]=2;
				my $str=$i+1;
				warn "\nThe position $str in the $ARGV[0] is totally consisted of gaps!";
			}else{
				$gaps_fracts[$I]->[$i]/=$seqs_counts[$I]->[$i];
			}
		}
	}
}
if(exists $args{z}){
	#remove branches without mutations
	$tree->visit_depth_first(  
        '-post' => sub {
			my $node=shift;
			my $name=$node->get_name;
			unless($node->is_terminal()||$node->is_root){
				my $n=0;
				if($f_intragene){
					$n+=@{$node->get_generic('syn')->[0]};
					$n+=@{$node->get_generic('nsyn')->[0]};
				}else{
					$n+=@{$node->get_generic('syn')->[1]};
					$n+=@{$node->get_generic('nsyn')->[1]};
				}
				$node->collapse() unless($n);
			}
		},
	);
}
#print XPARR
sub get_mutation_str{
	my ($node,$idx,$ra_site_nmuts,$MinMuts)=@_;
	$MinMuts=0 unless defined $MinMuts;
	die "\nWrong index value. Only 0 and 1 allowed!" unless $idx==0||$idx==1;
	if(!$node->is_root){
		my $ra_sites=$node->get_generic('nsyn')->[$idx];
		my @nsyn;
		my @mask=(1) x scalar(@{$ra_sites});
		if($MinMuts>1){
			for(my $i=0;$i<@mask;$i++){
				my $site=$ra_sites->[$i]->site;
				$mask[$i]=0 if $ra_site_nmuts->[$idx]->{$site}<$MinMuts
			}
		}
		if(defined($MaxGapsFrac)){
			for(my $i=0;$i<@mask;$i++){
				my $site=$ra_sites->[$i]->site-1;
				$mask[$i]=0 if $gaps_fracts[$idx]->[$site]>$MaxGapsFrac;
			}
		}
		my $j=0;
		for(my $i=0;$i<@mask;$i++){
			$nsyn[$j++]=$ra_sites->[$i] if $mask[$i];
		}
		my @tmp;
		foreach my $si(@nsyn){
			my $str=$si->site;
			if($f_print_mutations){
				my @bases=@{$si->bases};
				$str=$bases[0].$str.$bases[1];
			}
			push @tmp,$str;
		}
		return join ";",@tmp;
	}
	return "";
}

sub get_syn_mutation_str{
	my ($node,$idx)=@_;
	die "\nWrong index value. Only 0 and 1 allowed!" unless $idx==0||$idx==1;
	if(!$node->is_root){
		my $ra_sites=$node->get_generic('syn')->[$idx];
		return "" unless defined $ra_sites;
		my @syn;
		for(my $i=0;$i<@{$ra_sites};$i++){
			my $site=$ra_sites->[$i]->site;
			if(defined($MaxGapsFrac)){
				push @syn,$site unless $gaps_fracts[$idx]->[$site-1]>$MaxGapsFrac;
			}else{
				push @syn,$site;
			}
		}
		return join ";",@syn;
	}
	return "";
}

print "child\tparent\tlength";
foreach my $node($tree->get_nodes){
	my $name=$node->get_name;
	my $length=$node->get_branch_length;
	print "\n$name";
	if(!$node->is_root){
		my $pname=$node->get_parent->get_name;
		my $str="";
		print "\t$pname\t$length\t";
		my $str_bg=get_mutation_str($node,0,$ra_site_nmuts,$MinMuts);
		if($f_intragene){
			$str=get_syn_mutation_str($node,0) if $f_translate;
			print "$str\t";
			print "$str_bg\t$str_bg";
		}else{
			$str=get_syn_mutation_str($node,1) if $f_translate;
			print "$str\t";
			$str=get_mutation_str($node,1,$ra_site_nmuts,$MinMuts);
			print "$str\t$str_bg";
		}
	}
}
