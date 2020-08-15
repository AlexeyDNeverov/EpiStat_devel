package PValuesPOTapprox;
#This module provides functions for calculating statistics on set of trees with randomized mutation distributions
use strict;
use SitePairMatrix;
use Class::Struct;
use IndexedFakeSamplesHeap;
use Time::Progress;
use Scalar::Util qw(blessed);

my $threshr_R="~/devel/epistat.6/threshr.R";
#interface declaration
#constructor
sub new{
	my $class=shift;
	my $self={};
	bless($self,$class);
	$self->_init(@_);
	return $self;
}
struct SitePairStat =>{
	observed => '$',
	emp_pvalue => '$',
	idx => '$',
	gpd_params => '$',
	app_pvalue => '$'
};

struct GPD_params =>{
	threshold => '$',
	quantile => '$',
	shape => '$',
	scale => '$'
};

sub calc_smiths_app_pval{
	my ($x,$gpd_pars)=@_;
	my $epsilon=0.001;
	my $u=$gpd_pars->threshold;
	my $qu=$gpd_pars->quantile;
	my $ksi=$gpd_pars->shape;
	my $beta=$gpd_pars->scale;
print STDERR "\n\tu=$u\tqu=$qu\tksi=$ksi\tbeta=$beta";
	return undef unless defined($u)&&defined($qu)&&defined($ksi)&&defined($beta);
	die "\ncalc_smiths_app_pval(): the argument must be greater than $u!" unless $x>$u;
	if(abs($ksi)<$epsilon){
		return $qu*exp(-($x-$u)/$beta);
	}
	return $qu*(1+$ksi*($x-$u)/$beta)**(-1/$ksi);
}

sub pot{
	my $self=shift;
	my ($spi,$ra_data)=@_;
	my $nperm=@{$ra_data};
	my ($bgs,$fgs)=$self->{SITE_PAIR_MTX}->line2site_pair($spi->idx);
	my $j_max=$nperm-1;
	while($ra_data->[$j_max]>=$spi->observed){$j_max--;}
	die "\nThe POT approximation is not applicable to high (>0.01) empirical pvalues!" if $j_max/$nperm<=0.95;
	if($j_max>$nperm-50-1){$j_max=$nperm-50-1;}
	my $j=int($nperm/2);
	my $median=$ra_data->[$j];
	if($nperm%2==0){
		$median+=$ra_data->[$j+1];
		$median/=2;
	}
	my @quantiles;
	push @quantiles,[($median,$j,0)];
	$j=int(0.9*$nperm);
	while($ra_data->[$j]==$ra_data->[$j-1]||$ra_data->[$j]==$ra_data->[$j+1]){
		$j++;
		last if $j==$j_max;
	}
	if($j<$j_max){
		push @quantiles,[($ra_data->[$j],$j,0)];
	}
	$j=int(0.95*$nperm);
	while($ra_data->[$j]==$ra_data->[$j-1]||$ra_data->[$j]==$ra_data->[$j+1]){
		$j++;
		last if $j==$j_max;
	}
	if($j<$j_max){
		push @quantiles,[($ra_data->[$j],$j,0)];
	}
	my $q=0.95;
	my $qstep=0.02;
	$q+=$qstep;
	$j=int($q*$nperm);
	while($ra_data->[$j]==$ra_data->[$j-1]||$ra_data->[$j]==$ra_data->[$j+1]){
		$j++;
		last if $j==$j_max;
	}
	while($j<$j_max){
		push @quantiles,[($ra_data->[$j],$j,0)];
		$q+=$qstep;
		$j=int($q*$nperm);
		while($ra_data->[$j]==$ra_data->[$j-1]||$ra_data->[$j]==$ra_data->[$j+1]){
			$j++;
			last if $j>=$j_max;
		}
	}
	push @quantiles,[($ra_data->[$j_max],$j_max,0)];
	my $n=1;
	$j=@quantiles-1;
	for(my $i=$nperm-2;$i>=$quantiles[1]->[1];$i--){
		$n++ if($ra_data->[$i]<$ra_data->[$i+1]);
		if($i==$quantiles[$j]->[1]){
			$quantiles[$j]->[2]=$n;
			$j--;
		}
	}
	my @tmp;
	push @tmp,$quantiles[0]->[0];
	for($j=1;$j<@quantiles;$j++){
		push @tmp,$quantiles[$j]->[0] if ($quantiles[$j]->[0]>$tmp[-1])&&
			($quantiles[$j]->[2]/($nperm-$quantiles[$j]->[1])>0.90);
	}
	@quantiles=@tmp[1 .. $#tmp];
	if(@quantiles>1){
		my $tmp=gen_tempname(5).".$bgs"."_"."$fgs";
		open OPF, ">$tmp.tmp" or die "\nUnable to open output file: $tmp.tmp!";
		print OPF $quantiles[1];
		for(my $i=2;$i<@quantiles;$i++){
			print OPF "\t".$quantiles[$i];
		}
		print OPF "\n";
		print OPF $ra_data->[0];
		for(my $i=1;$i<$nperm;$i++){
			print OPF "\t".$ra_data->[$i];
		}
		print OPF "\n";
		close OPF;
		my $str="$threshr_R $tmp.tmp>$tmp.R.out";
		system($str);
		print_child_termination_status();
		open OPF, "<$tmp.R.out" or die "\nUnable to open input file $tmp.R.out!";
		my ($u_idx,$qu,$ksi,$beta);
		while(<OPF>){
			chomp;
			if(/(\w+)\s*=\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)/){
				my ($key,$value)=($1,$2);
				if($key eq "best_u_idx"){
					$u_idx=$value-1;
				}elsif($key eq "mean_Qu"){
					$qu=$value;
				}elsif($key eq "mean_beta"){
					$beta=$value
				}elsif($key eq "mean_ksi"){
					$ksi=$value;
				}
			}
		}
		close OPF;
		if(defined($u_idx)&&defined($qu)&&defined($beta)&&defined($ksi)){
			my $gpd_pars=GPD_params->new();
			$gpd_pars->threshold($quantiles[$u_idx]) if ($u_idx>=0)&&($u_idx<@quantiles);
			$gpd_pars->quantile($qu);
			$gpd_pars->shape($ksi);
			$gpd_pars->scale($beta);
			my $app_pval=calc_smiths_app_pval($spi->observed,$gpd_pars);
			$spi->gpd_params($gpd_pars);
			$spi->app_pvalue($app_pval) if defined $app_pval;
		}
		unlink "$tmp.tmp";
		unlink "$tmp.R.out";
	}
}

sub _init{
	my $self=shift;
	my %args = @_;
	my $nperm=0;
	my ($samples_dir,$stats_ext,$stats_fn,$pairs_fn,$f_intragene,$pvalues_fn);
	my $max_samples_indir=10000000;
	my $samples_heap;
	my $f_square_mtx;
	my @line2compl_line;
	my @pvalues;
	my $pval_POThreshold;
	my $f_ordered_pairs;
	while(my ($k,$v)=each %args){
		$v=~s/^\s*//;
		$v=~s/^\s*$//;
		if($k eq "-size"){
			$nperm=$v if $v=~m/^\d+/;
		}elsif($k eq "-samples_dir"){
			$v=~s/\/$//;
			$samples_dir=$v;
		}elsif($k eq "-max_samples_indir"){
			$max_samples_indir=$v;
		}elsif($k eq "-stats_ext"){
			$stats_ext=$v;
		}elsif($k eq "-obs_stats_fn"){
			$stats_fn=$v;
		}elsif($k eq "-obs_pairs_fn"){
			$pairs_fn=$v;
		}elsif($k eq "-is_intragene"){
			$f_intragene=$v;
		}elsif($k eq "-pairs_ordering"){
			if(defined($v)&&($v ne "")&&($v==0||$v==1)){
				$f_ordered_pairs=$v;
			}else{
				die "\n\nPValuesPOTApprox::init_(): unvalid value $v was accounted for the parameter '-pairs_ordering': 0 or 1 allowed!";
			}
		}elsif($k eq "-pvalues_fn"){
			$pvalues_fn=$v;
		}elsif($k eq "-pvalue_POThreshold"){
			#Set the upper threshold on empirically estimated pvalues to be approximated by POT
			$pval_POThreshold=$v;
			die "\nPValuesPOTApprox::init_(): The POT approximation is applicable only for high quantiles of a PDF: set -calc_POT_pvalues within the [0,0.5] range!" 
				unless $pval_POThreshold>=0&&$pval_POThreshold<=0.01;
		}else{
			die "\nPValuesPOTApprox::init_(): Unknown parameter: $k!";
		};
	}
	die "\nPValuesPOTApprox::init_(): The POT pvalue approximation could not be applied to all pairs. Set '-pvalue_POThreshold' and '-pvalues_fn'!" unless defined($pval_POThreshold)&&defined($pvalues_fn);
	$f_ordered_pairs=1 unless defined $f_ordered_pairs;
	my @samples_storage;
	if($nperm>$max_samples_indir){
		$samples_heap=IndexedFakeSamplesHeap->new($nperm,$max_samples_indir);
		my @tmp=$samples_heap->get_storage_path_list;
		for(my $i=0;$i<@tmp-1;$i++){
			push @samples_storage,[($max_samples_indir*$i+1,$max_samples_indir*($i+1),"/".$tmp[$i])];
		}
		push @samples_storage,[($max_samples_indir*(@tmp-1)+1,$nperm,"/".$tmp[-1])];
	}else{
		push @samples_storage,[(1,$nperm,"")];
	}
	die "\nParameter -stats_ext has no default value!" unless defined $stats_ext;
	die "\nThe number of samples have to exceed the 50 to apply POT approximation!" unless $nperm>50;
	$self->{SAMPLE_SIZE}=$nperm;
	$nperm-- if defined $self->{SKIP_ID};
	die "\nError PValuesPOTApprox::_init(): No permutations sampled!" unless $nperm>0;
	$self->{OBS_STAT}={};
	$self->{PDF}=[];
	$self->{SITE_PAIR_MTX}=undef;
	#########################
	open INPF, "<$stats_fn" or die "\nUnable to open input file:$stats_fn!";
	my $nlines=0;
	while(<INPF>){
		chomp;
		my @line=split "\t";
		if($line[0] ne ""){
			$self->{OBS_STAT}->{$nlines}=$line[0] if($line[0]>0);
			$nlines++;
		}
	}
	close INPF;
	open INPF, "<$pvalues_fn" or die "\nUnable to open input file:$pvalues_fn!";
	while(<INPF>){
		chomp;
		my @line=split "\t";
		if($line[0] ne ""){
			push @pvalues, $line[0];
		}
	}
	close INPF;
	die "\nPValuesPOTApprox::_init(): Wrong number of records in $pvalues_fn file!" unless @pvalues==$nlines;
	#########################
	if(defined($pairs_fn)&&defined($f_intragene)){
		$self->{SITE_PAIR_MTX}=SitePairMatrix->new();
		$self->{SITE_PAIR_MTX}->init_sites($pairs_fn,$f_intragene);
		$f_square_mtx=$self->{SITE_PAIR_MTX}->is_square;
		for(my $j=0;$j<$nlines;$j++){
			my ($bgr_idx,$fgr_idx)=$self->{SITE_PAIR_MTX}->line2sites_idxs($j);
			if($f_square_mtx){
				my $symj=$self->{SITE_PAIR_MTX}->site_idx_pair2line($fgr_idx,$bgr_idx);
				push @line2compl_line,$symj;
			}
		}
	}
	print STDERR "\nPValuesPOTApprox::init_(): Start reading data files:";
	my $p = Time::Progress->new(min => 1, max => $self->{SAMPLE_SIZE});
	my @pairs_idx;
	my @pairs_data;
	for(my $i=0;$i<$nlines;$i++){
		if($pvalues[$i]<=$pval_POThreshold){
			if($f_ordered_pairs){
				push @pairs_idx,$i;
				push @pairs_data,[];
			}elsif($line2compl_line[$i]>$i){
				push @pairs_idx,$i;
				push @pairs_data,[];
			}
		}
	}
	die "\nNo pairs were found <=$pval_POThreshold" unless scalar(@pairs_data);
	for(my $si=0;$si<@samples_storage;$si++){
		my $from=$samples_storage[$si]->[0];
		my $to=$samples_storage[$si]->[1];
		my $folder=$samples_storage[$si]->[2];
		for(my $i=$from;$i<=$to;$i++){
			print STDERR $p->report("\r%20b  ETA: %E", $i);
			my $str=$samples_dir.$folder."/".$i.$stats_ext;
			my @stat_buff;
			open INPF, "<$str" or die "\nUnable to open input file:$str!";
			while(<INPF>){
				chomp;
				my @line=split "\t";
				push @stat_buff, $line[0] if($line[0] ne "");
			}
			close INPF;
			die "\nPValuesPOTApprox::init_(): Unexpected number of records in the file $str!" unless @stat_buff==$nlines;
			my $ip=0;
			for(my $j=0;$j<$nlines;$j++){
				if($j==$pairs_idx[$ip]){
					my $stat=$stat_buff[$j];
					unless($f_ordered_pairs){
						my $symj=$line2compl_line[$j];
						$stat+=$stat_buff[$symj];
					}
					push @{$pairs_data[$ip]},$stat;
					$ip++;
				}
			}
		}
	}
	
	for(my $i=0;$i<@pairs_data;$i++){
		my $ra=$pairs_data[$i];
		@{$ra}=sort {$a<=>$b} @{$ra};
		my $spi=SitePairStat->new();
		my $spi=SitePairStat->new();
		my $line=$pairs_idx[$i];
		$spi->idx($line);
		$spi->emp_pvalue($pvalues[$line]);
		my $obs=$self->{OBS_STAT}->{$line};
		if(defined $obs){
			$spi->observed($obs);
		}else{
			$spi->observed(0);
		}
		unless($f_ordered_pairs){
			my $symi=$line2compl_line[$line];
			$obs=$self->{OBS_STAT}->{$symi};
			$obs=0 unless defined $obs;
			$obs+=$spi->observed;
			$spi->observed($obs);
		}
		$self->pot($spi,$ra);
		push @{$self->{PDF}},$spi;
	}
	print STDERR "\nPValuesPOTApprox::init_(): done!";
}

sub gen_tempname{
	my $nchar=shift;
	my @chars = ( "A" .. "Z", "a" .. "z", 0 .. 9 );
	return join("", @chars[ map { rand @chars } ( 1 .. $nchar ) ]);
}

sub print_child_termination_status{
	if ($? == -1) {
		print "failed to execute: $!\n";
		exit 1;
	}elsif ($? & 127) {
		printf "\tchild died with signal %d, %s coredump\n",
			($? & 127),  ($? & 128) ? 'with' : 'without';
	}else {
		if($? >> 8!=0){
			printf "\tchild exited with value %d\n", $? >> 8;
		}
	}
}

sub is_square{
	my $self=shift;
	return undef unless defined $self->{SITE_PAIR_MTX};
	return $self->{SITE_PAIR_MTX}->is_square;
}

sub print{
	my $self=shift;
	print "bgr_site\tfgr_site\tobserved\temp_pvalue\tapp_pvalue\tu\tqu\tksi\tbeta";
	foreach my $spi(@{$self->{PDF}}){
		my ($bgs,$fgs)=$self->{SITE_PAIR_MTX}->line2site_pair($spi->idx);
		print "\n$bgs\t$fgs";
		print "\t".$spi->observed;
		print "\t".$spi->emp_pvalue;
		print "\t".$spi->app_pvalue;
		if(defined $spi->gpd_params){
			print "\t".$spi->gpd_params->threshold;
			print "\t".$spi->gpd_params->quantile;
			print "\t".$spi->gpd_params->shape;
			print "\t".$spi->gpd_params->scale;
		}else{
			print "\t\t\t\t";
		}
		print "\n";
	}
}

1;