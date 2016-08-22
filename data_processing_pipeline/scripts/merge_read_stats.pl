##############################################################################################
## Sasquatch, Sequence based predicting of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##############################################################################################
# Function: Data preprocessing; Merge read stats from multiple replicate runs
# Usage: perl ./merge_read_stats.pl ./data1/read_stats.txt ./data2/read_stats.txt | tee ./data_merged/read_stats.txt
# Author: Ron Schwessinger
# Date: 27/04/2016

#!/usr/bin/perl -w
use strict;
use warnings;

#merge multiple replicates read_stats files (read number, peak number, reads in peaks)

my (@col, $line, %hash);

my $number = 0;

$hash{"total"} = 0;
$hash{"peaks"} = 0;
$hash{"inpeaks"} = 0;


while (@ARGV) {

    my $input = shift(@ARGV);
    
    $number++;
    
    #open file
    open(IN, $input) or die "Cant open(input file $input, $!\n";
    
     while (<IN>){
        
            chomp;
            
            $hash{"total"} = $hash{"total"} + $1 if $_ =~ /Total reads: (\d+)/;
            
            $hash{"peaks"} = $hash{"peaks"} + $1 if $_ =~ /Number of Peaks:: (\d+)\s+/;
            
            $hash{"inpeaks"} = $hash{"inpeaks"} + $1 if $_ =~ /Reads in Peaks:: (\d+)/;
            
     }

}     
     
close(IN);

#get average values
$hash{"total"} = $hash{"total"} / $number;
$hash{"peaks"} = $hash{"peaks"} / $number;
$hash{"inpeaks"} = $hash{"inpeaks"} / $number;

#print out
print "Total reads: ".$hash{"total"}."\n";
print "Number of Peaks:: ".$hash{"peaks"}."\n";
print "Reads in Peaks:: ".$hash{"inpeaks"}."\n";
