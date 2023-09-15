#!/usr/bin/env perl
#Select a minimal elements for two vectors
#Usage: <column1_fn>	<column2_fn>
open INPF1, "<$ARGV[0]" or die "\nUnable to open vector file 1: $ARGV[0]!";
open INPF2, "<$ARGV[1]" or die "\nUnable to open vector file 2: $ARGV[1]!";
my ($n1,$n2)=(0,0);
while(<INPF1>){
	if(/\S+/){
		s/^\s+//;
		s/\s+$//;
		if(/^[-+]?([0-9]+(\.[0-9]+)?|\.[0-9]+)$/){
			my $f1=$_;
			$n1++;
			while(<INPF2>){
				if(/\S+/){
					s/^\s+//;
					s/\s+$//;
					if(/^[-+]?([0-9]+(\.[0-9]+)?|\.[0-9]+)$/){
						my $f2=$_;
						my $min=$f1<$f2?$f1:$f2;
						print "$min\n";
						$n2++;
						last;
					}else{
						die "\nIn the vector files nonnumeric values are not allowed!";
					}
				}
			}
		}else{
			die "\nIn the vector files nonnumeric values are not allowed!";
		}
	}
}
die "The numbers of values in the input vectors do not match: $n1<>$n2!" unless $n1==$n2;