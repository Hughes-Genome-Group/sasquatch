##############################################################################################
## Sasquatch, Sequence based predicting of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##############################################################################################
# Function: Data preprocessing; Merge kmer profiles from naked chromosome profiles
# Usage: perl ./merge_naked_chrs.pl chr1_profiles.txt chr2_profiles.txt >merged_profiles.txt
# Author: Ron Schwessinger
# Date: 27/04/2016

#!/usr/bin/perl -w
use strict;
use warnings;

#merge multiple naked chromosome cuts together

my (@col, $line, %hash);

while (@ARGV) {
    
    my $input = shift(@ARGV);
    
    #open file
    open(IN, $input) or die "Cant open(input file $input, $!\n";
    
        while (<IN>){
        
            chomp;
            @col = split (/\s+/, $_);
            my $kmer = shift @col;
            
            if (exists $hash{$kmer}) {
                my @temp = split(/\t+/, $hash{$kmer});#split up existsing hash
                my $length = @temp;
                $length--;
                for(my $i=0; $i<=$length; $i++){ #add values
                    $temp[$i] += $col[$i];
                }    
                $hash{$kmer}=join("\t", @temp);    #restore
            }else{
                $hash{$kmer} = join("\t", @col);
            }
            
        }
    
    close(IN);    
    
}

#print combined
my @keys = keys(%hash);
@keys = sort(@keys);

foreach (@keys){
    print "$_\t".join("\t", $hash{$_})."\n";
}

