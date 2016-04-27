##############################################################################################
## Sasquatch, Sequence based predicting of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##############################################################################################
# Function: Data preprocessing; get first and last DNase cut covered position from wig file
# Usage: perl ./get_footprint_wig_borders.pl input_footprint.wig >out_footprint_borders
# Author: Ron Schwessinger
# Date: 27/04/2016

#!/usr/bin/perl -w
use strict;
use warnings;

#get borders of a wig file (footprints)
my $input=$ARGV[0];

my ($line, @col,$temp,$prev);

open(IN, $input) or die "Cant open input! $! \n";

    while (<IN>) {
        
        if ($_=~/variable/) {
            $prev =~/^(\d+)\s/;
            print "$1\n" if $prev ne "";
            $_=~/chrom=(chr[\d,X,Y,M]+)\s/;
            print "$1\t";
            $temp=<IN>;
            $temp =~/^(\d+)\s/;
            print "$1\t";
        }

        
    $prev =$_;
    
    }
print "\n";
close(IN);
