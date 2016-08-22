##############################################################################################
## Sasquatch, Sequence based predicting of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##############################################################################################
# Function: Data preprocessing; generate all possible kmers after setting kmer length to STDOUT.
# Author: Jim Hughes
# Date: 27/04/2016

#!/usr/bin/perl -w
use Getopt::Long;
use List::Util qw( min max );
use strict;

my $lengths = 6;
my $Name = 'Standard_output';

&GetOptions (
         "lengths=s"=> \ $lengths,
		 "name=s"=> \ $Name,
	     );


my @DesLength = split(/,/, $lengths);

my $length = max @DesLength;

my @bases = qw(G A T C);
my %sequences;
my %evolving;

my $counter = 1;

foreach my $base (@bases){
	$sequences{$counter}{$base}=$base;
}


until ($counter == $length + 1){

	my $CountNext = $counter +1;

	foreach my $StoredSeq (sort keys %{$sequences{$counter}}){
	
		foreach my $base (@bases){
			my $newsequence = "$StoredSeq" . "$base";
			$sequences{$CountNext}{$newsequence}=$newsequence;
		}
	
}

$counter++;
}


foreach my $RequiredLen (@DesLength){

	my $generatedCount = 0;
	open (OUTPUT, ">$Name\_$RequiredLen.txt");

	foreach my $Finalsequences (sort keys %{$sequences{$RequiredLen}}){
		$generatedCount++;
		print OUTPUT "$Finalsequences\n";
	}

	print "\nGenerated $generatedCount sequences of length $RequiredLen\n";
}
