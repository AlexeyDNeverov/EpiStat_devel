#This package implements aggregation methods for epistatic statistics obtained for different null-models conditioned on phenotypes
package Phenotype::AggregateStatistics;
use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK $VERSION);
use Exporter;
$VERSION = 1.00; # Or higher
@ISA = qw(Exporter);
# Symbols to autoexport (:DEFAULT tag)
@EXPORT = qw(); 
@EXPORT_OK = qw();
use File::Path qw(make_path rmtree);
use SitePairMatrix;
use Data::Float;
use Time::Progress;
#sort_order: 0- ascending; 1- descending
sub get_best_stat{
	my ($sort_order,$s1,$s2)=@_;
	if(defined($s1)&&defined($s2)){
		if($sort_order==0){
			return ($s1<=$s2)?$s1:$s2;
		}
		return ($s1>=$s2)?$s1:$s2;
	}
	return defined($s1)?$s1:$s2;
}

sub find_missed_samples{
	my ($path,$ext,$ra_samples_storage_node)=@_;
	my @n;
	$path=~s/\s*$//;
	#$path=~s/\/$//;
	#$path.=$ra_samples_storage_node->[2];
	$path.="/";
	my $from=$ra_samples_storage_node->[0];
	my $to=$ra_samples_storage_node->[1];
	for(my $i=$from;$i<=$to;$i++){
		my $str=$path.$i.$ext;
		if(-e $str){
			my $filesize=-s $str;
			if($filesize){
				last if(defined $n[1]);
				$n[0]++;
			}else{
				$n[1]++;
			}
		}else{
			$n[1]++;
		}
	}
	if (wantarray()) {
		# list context
		$n[0]+=$from;
		$n[1]+=$n[0]-1;
		return @n;
	}elsif (defined wantarray()) {
		return $n[1];
	}
	return undef;
}

sub get_random_index{
	my @probs=@_;
	die "\nThe function expects a vector of cumulative probabilities as an argument!" unless(sprintf("%.6f",$probs[-1])==1);
	my $k=0;
	if(@probs>1){
		my $smpl=rand;
		for(;$k<@probs;$k++){
			last if $smpl<$probs[$k];
		}
	}
	return $k;
}

sub aggregateByBest{
	my ($pairs_fn,$f_intragene,$indir,$ra_pheno_samles_dirs,$ra_samples_storage,$file_ext,
			$rh_bg_sites2pheno,$rh_fg_sites2pheno,$site2pheno_order,$outdir,$ntries)=@_;
	$ntries=1 unless defined $ntries;
	$indir=~s/\s+$//;
	$indir=~s/\/$//;
	$indir.="/";
	$outdir=~s/\s+$//;
	$outdir=~s/\/$//;
	my $sp_matrix=SitePairMatrix->new($pairs_fn,$f_intragene);
	my @bg_sites;
	my @fg_sites;
	$sp_matrix->get_sites(\@bg_sites,\@fg_sites);
	for(my $i=0;$i<@bg_sites;$i++){
		die "\nThe site $bg_sites[$i] is absent in the background site to phenotype map!" unless defined $rh_bg_sites2pheno->{$bg_sites[$i]};
	}
	for(my $i=0;$i<@fg_sites;$i++){
		die "\nThe site $fg_sites[$i] is absent in the background site to phenotype map!" unless defined $rh_fg_sites2pheno->{$fg_sites[$i]};
	}
	my $npairs=$sp_matrix->{NLINES};
	my $npheno=@{$ra_pheno_samles_dirs};
	my $nsamples=$ra_samples_storage->[-1]->[1]-$ra_samples_storage->[0]->[0]+1;
	my $progress_bar = Time::Progress->new(min => 1, max => $nsamples);
	die "\nNot allowed value for the 'site2pheno_order': $site2pheno_order!" unless ($site2pheno_order==0)||($site2pheno_order==1);
	my @pairs2pheno;
	for(my $i=0;$i<$npairs;$i++){
		my ($bgs,$fgs)=$sp_matrix->line2site_pair($i);
		my $p2p_val=-(Data::Float::max_finite-1.0); #-inf
		if(defined($rh_bg_sites2pheno->{$bgs}->[0])&&defined($rh_fg_sites2pheno->{$fgs}->[0])){
			$p2p_val=$rh_bg_sites2pheno->{$bgs}->[0]+$rh_fg_sites2pheno->{$fgs}->[0];
		}
		my @p2p_idx;
		push @p2p_idx,0;
		for(my $j=1;$j<$npheno;$j++){
			my $p2p_val2=$rh_bg_sites2pheno->{$bgs}->[$j]+$rh_fg_sites2pheno->{$fgs}->[$j];
			if(get_best_stat($site2pheno_order,$p2p_val,$p2p_val2)==$p2p_val2){
				if($p2p_val==$p2p_val2){
					push @p2p_idx,$j;
				}else{
					@p2p_idx=($j);
				}
				$p2p_val=$p2p_val2;
			}
		}
		push @pairs2pheno,[@p2p_idx];
	}
	for(my $j=0;$j<@{$ra_samples_storage};$j++){
		my $out_path=$outdir.$ra_samples_storage->[$j]->[2];
		make_path($out_path) unless(-d $out_path);
		for(my $I=0;$I<$ntries;$I++){
			my ($from,$to)=find_missed_samples($out_path,$file_ext,$ra_samples_storage->[$j]);
			last if $to<$from;
			for(my $k=$from;$k<=$to;$k++){
				my %aggregated;
				my @pairs2phen_idx;
				my $infile_header;
				for(my $i=0;$i<$npairs;$i++){
					my $ridx=0;
					my $n=@{$pairs2pheno[$i]};
					$ridx=int rand($n) if $n>1;
					$ridx=$pairs2pheno[$i]->[$ridx];
					push @pairs2phen_idx,$ridx;
				}
				for(my $i=0;$i<$npheno;$i++){
					my $in_path=$indir.$ra_pheno_samles_dirs->[$i].$ra_samples_storage->[$j]->[2];
					my $fn=$in_path."/".$k.$file_ext;
					open INPF, "<$fn" or die "\nUnable to open input file: $fn!";
					while(<INPF>){
						if(/^(\d+)\t/){
							my $pair_idx=$1;
							$aggregated{$pair_idx}=$_ if $pairs2phen_idx[$pair_idx]==$i;
						}
						if($i==0){
							if(/^NMaxDataLines/||/^sp_index/){
								$infile_header.=$_;
							}
						}
					}
					close INPF;
				}
				my $fn=$out_path."/".$k.$file_ext;
				open OPF,">$fn" or die "\nUnable to open output file:$fn!";
				print OPF $infile_header; 
				foreach my $pair_idx(sort {$a<=>$b} keys %aggregated){
					print OPF $aggregated{$pair_idx};
				}
				close OPF;
				print STDERR $progress_bar->report("\r%20b  ETA: %E", $k);
			}
		}
	}
}

#This function makes normalization of a vector of logarithms of values so that the vector exponentiation would be summed up to the one
#undefined values are interpreted as '-inf' values
sub log_of_simplex_closure{
	my $ra_lvec=shift;
	my @ln_vec=@{$ra_lvec};
	my $i=0;
	for(;$i<@ln_vec;$i++){
		last if defined $ln_vec[$i]; #undef means -inf
	}
	return (undef) x @ln_vec unless $i<@ln_vec;
	my $mval=$ln_vec[$i++];
	for(;$i<@ln_vec;$i++){
		next unless defined $ln_vec[$i];
		$mval=$ln_vec[$i] if $mval<$ln_vec[$i];
	}
	my $s=0;
	for(my $i=0;$i<@ln_vec;$i++){
		next unless defined $ln_vec[$i];
		$s+=exp($ln_vec[$i]-$mval);
	}
	my $ldenom=$mval+log($s);
	for(my $i=0;$i<@ln_vec;$i++){
		next unless defined $ln_vec[$i];
		$ln_vec[$i]-=$ldenom;
	}
	return @ln_vec;
}

sub aggregateByMean{
	my ($pairs_fn,$f_intragene,$indir,$ra_pheno_samles_dirs,$ra_samples_storage,$file_ext,
			$rh_bg_sites2pheno,$rh_fg_sites2pheno,$site2pheno_order,$outdir,$ntries)=@_;
	$ntries=1 unless defined $ntries;
	$indir=~s/\s+$//;
	$indir=~s/\/$//;
	$indir.="/";
	$outdir=~s/\s+$//;
	$outdir=~s/\/$//;
	my $sp_matrix=SitePairMatrix->new($pairs_fn,$f_intragene);
	my @bg_sites;
	my @fg_sites;
	$sp_matrix->get_sites(\@bg_sites,\@fg_sites);
	for(my $i=0;$i<@bg_sites;$i++){
		die "\nThe site $bg_sites[$i] is absent in the background site to phenotype map!" unless defined $rh_bg_sites2pheno->{$bg_sites[$i]};
	}
	for(my $i=0;$i<@fg_sites;$i++){
		die "\nThe site $fg_sites[$i] is absent in the background site to phenotype map!" unless defined $rh_fg_sites2pheno->{$fg_sites[$i]};
	}
	my $npairs=$sp_matrix->{NLINES};
	my $npheno=@{$ra_pheno_samles_dirs};
	my $nsamples=$ra_samples_storage->[-1]->[1]-$ra_samples_storage->[0]->[0]+1;
	my $progress_bar = Time::Progress->new(min => 1, max => $nsamples);
	die "\nNot allowed value for the 'site2pheno_order': $site2pheno_order!" unless ($site2pheno_order==0)||($site2pheno_order==1);
	my @pairs2pheno;
	for(my $i=0;$i<$npairs;$i++){
		my ($bgs,$fgs)=$sp_matrix->line2site_pair($i);
		push @pairs2pheno,[];
		my $norm=0;
		for(my $j=0;$j<$npheno;$j++){
			my $val;
			if(defined($rh_bg_sites2pheno->{$bgs}->[$j])&&defined($rh_fg_sites2pheno->{$fgs}->[$j])){
				$val=$rh_bg_sites2pheno->{$bgs}->[$j]+$rh_fg_sites2pheno->{$fgs}->[$j];
			}
			push @{$pairs2pheno[-1]},$val;
		}
		@{$pairs2pheno[-1]}=log_of_simplex_closure($pairs2pheno[-1]);
		if($site2pheno_order==0){
			for(my $j=0;$j<$npheno;$j++){
				if(defined $pairs2pheno[-1]->[$j]){
					$pairs2pheno[-1]->[$j]=log(1-exp($pairs2pheno[-1]->[$j]));
				}else{
					$pairs2pheno[-1]->[$j]=0;
				}
			}
			@{$pairs2pheno[-1]}=log_of_simplex_closure($pairs2pheno[-1]);
		}
	}
	for(my $j=0;$j<@{$ra_samples_storage};$j++){
		my $out_path=$outdir.$ra_samples_storage->[$j]->[2];
		make_path($out_path) unless(-d $out_path);
		for(my $I=0;$I<$ntries;$I++){
			my ($from,$to)=find_missed_samples($out_path,$file_ext,$ra_samples_storage->[$j]);
			last if $to<$from;
			for(my $k=$from;$k<=$to;$k++){
				my %aggregate;
				my $infile_header;
				for(my $i=0;$i<$npheno;$i++){
					my $in_path=$indir.$ra_pheno_samles_dirs->[$i].$ra_samples_storage->[$j]->[2];
					my $fn=$in_path."/".$k.$file_ext;
					open INPF, "<$fn" or die "\nUnable to open input file: $fn!";
					while(<INPF>){
						chomp;
						s/\s+$//;
						if($i==0){
							if(/^NMaxDataLines/||/^sp_index/){
								$infile_header.=$_;
							}
						}
						if(/^\d+\t/){
							my @line=split "\t";
							my $pair_idx=$line[0];
							$aggregate{$pair_idx}=[] unless defined $aggregate{$pair_idx};
							for(my $j=1;$j<@line;$j++){
								if(defined $pairs2pheno[$pair_idx]->[$i]){
									$aggregate{$pair_idx}->[$j-1]+=$line[$j]*exp($pairs2pheno[$pair_idx]->[$i]);
								}
								if($j>1){
									$aggregate{$pair_idx}->[$j-1]=sprintf("%.2f",$aggregate{$pair_idx}->[$j-1]);
								}else{
									$aggregate{$pair_idx}->[0]=sprintf("%.6f",$aggregate{$pair_idx}->[0]);
								}
							}
						}
						
					}
					close INPF;
				}
				my $fn=$out_path."/".$k.$file_ext;
				open OPF,">$fn" or die "\nUnable to open output file:$fn!";
				print OPF $infile_header;
				foreach my $pair_idx(sort {$a<=>$b} keys %aggregate){
					my $line=$pair_idx."\t".join("\t",@{$aggregate{$pair_idx}})."\n";
					print OPF $line;
				}
				close OPF;
				print STDERR $progress_bar->report("\r%20b  ETA: %E", $k);
			}
		}
	}
}

#$site2pheno_order - best to worst site to phenotype association statistics sorting order: 0 - for p-values; 1 -for z-scores
sub aggregateByMixture{
	my ($pairs_fn,$f_intragene,$indir,$ra_pheno_samles_dirs,$ra_samples_storage,$file_ext,
			$rh_bg_sites2pheno,$rh_fg_sites2pheno,$site2pheno_order,$outdir,$ntries)=@_;
	$ntries=1 unless defined $ntries;
	$indir=~s/\s+$//;
	$indir=~s/\/$//;
	$indir.="/";
	$outdir=~s/\s+$//;
	$outdir=~s/\/$//;
	my $sp_matrix=SitePairMatrix->new($pairs_fn,$f_intragene);
	my @bg_sites;
	my @fg_sites;
	$sp_matrix->get_sites(\@bg_sites,\@fg_sites);
	for(my $i=0;$i<@bg_sites;$i++){
		die "\nThe site $bg_sites[$i] is absent in the background site to phenotype map!" unless defined $rh_bg_sites2pheno->{$bg_sites[$i]};
	}
	for(my $i=0;$i<@fg_sites;$i++){
		die "\nThe site $fg_sites[$i] is absent in the background site to phenotype map!" unless defined $rh_fg_sites2pheno->{$fg_sites[$i]};
	}
	my $npairs=$sp_matrix->{NLINES};
	my $npheno=@{$ra_pheno_samles_dirs};
	my $nsamples=$ra_samples_storage->[-1]->[1]-$ra_samples_storage->[0]->[0]+1;
	my $progress_bar = Time::Progress->new(min => 1, max => $nsamples);
	die "\nNot allowed value for the 'site2pheno_order': $site2pheno_order!" unless ($site2pheno_order==0)||($site2pheno_order==1);
	my @pairs2phen_probs;
	for(my $i=0;$i<$npairs;$i++){
		my ($bgs,$fgs)=$sp_matrix->line2site_pair($i);
		my @phen_probs;
		for(my $j=0;$j<$npheno;$j++){
			my $val;
			if(defined($rh_bg_sites2pheno->{$bgs}->[$j])&&defined($rh_fg_sites2pheno->{$fgs}->[$j])){
				$val=$rh_bg_sites2pheno->{$bgs}->[$j]+$rh_fg_sites2pheno->{$fgs}->[$j];
			}
			push @phen_probs,$val;
		}
		@phen_probs=log_of_simplex_closure(\@phen_probs);
		if($site2pheno_order==0){
			for(my $j=0;$j<$npheno;$j++){
				$phen_probs[$j]=-$phen_probs[$j];
			}
			@phen_probs=log_of_simplex_closure(\@phen_probs);
			if($site2pheno_order==0){
				for(my $j=0;$j<$npheno;$j++){
					if(defined $phen_probs[$j]){
						$phen_probs[$j]=exp($phen_probs[$j]);
					}else{
						$phen_probs[$j]=0;
					}
					$phen_probs[$j]+=$phen_probs[$j-1] if $j;
				}
			}
		}
		push @pairs2phen_probs,[@phen_probs];
	}
	for(my $j=0;$j<@{$ra_samples_storage};$j++){
		my $out_path=$outdir.$ra_samples_storage->[$j]->[2];
		make_path($out_path) unless(-d $out_path);
		for(my $I=0;$I<$ntries;$I++){
			my ($from,$to)=find_missed_samples($out_path,$file_ext,$ra_samples_storage->[$j]);
			last if $to<$from;
			for(my $k=$from;$k<=$to;$k++){
				my %aggregated;
				my @pairs2phen_idx;
				my $infile_header;
				for(my $i=0;$i<$npairs;$i++){
					my $ridx=get_random_index(@{$pairs2phen_probs[$i]});
					push @pairs2phen_idx,$ridx;
				}
				for(my $i=0;$i<$npheno;$i++){
					my $in_path=$indir.$ra_pheno_samles_dirs->[$i].$ra_samples_storage->[$j]->[2];
					my $fn=$in_path."/".$k.$file_ext;
					open INPF, "<$fn" or die "\nUnable to open input file: $fn!";
					my $pair_idx=0;
					while(<INPF>){
						if(/^(\d+)\t/){
							my $pair_idx=$1;
							$aggregated{$pair_idx}=$_ if $pairs2phen_idx[$pair_idx]==$i;
						}
						if($i==0){
							if(/^NMaxDataLines/||/^sp_index/){
								$infile_header.=$_;
							}
						}
					}
					close INPF;
				}
				my $fn=$out_path."/".$k.$file_ext;
				open OPF,">$fn" or die "\nUnable to open output file:$fn!";
				print OPF $infile_header; 
				foreach my $pair_idx(sort {$a<=>$b} keys %aggregated){
					print OPF $aggregated{$pair_idx};
				}
				close OPF;
				print STDERR $progress_bar->report("\r%20b  ETA: %E", $k);
			}
		}
	}
}


1;