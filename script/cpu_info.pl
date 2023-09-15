#!/usr/bin/env perl

open CPU, "/proc/cpuinfo" or die "Can't open cpuinfo: $!\n";
printf "CPUs: %d\n", scalar (map /^processor/, <CPU>) ; 
close CPU;
system("free -m");