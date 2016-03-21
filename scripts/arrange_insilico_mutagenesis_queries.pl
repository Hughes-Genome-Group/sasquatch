#!/usr/bin/perl -w
use strict;
use warnings;

#for each longer sequence split sequence into kl long kmers and get for every sliding middle base all possible variations (SNPS). Assemble id snpbase refseq varseq (both 6bases surrounding thefor table and pruning to query

my $seq = $ARGV[0];

chomp $seq;

my $stringlength = length($seq);	#get seq length

my $kl = 7; #kmer length to split into

#initialize
my $id = 0;

#move along sequence
for(my $k = 0; $k < $stringlength - 12; $k++){

	$id++;
 	
	#get substring split of the sequence (-+ 6bp surrounding the center)
	my $ref = substr($seq, $k, 13);

	my $refbase = substr($seq,(6+$k),1); #extract reference base (center)

	my @varbases = ();

	if($refbase eq "A"){	#according to refbase select possible SNPs
		@varbases = ("C","G","T");
	}elsif($refbase eq "C"){
		@varbases = ("A","G","T");
	}elsif($refbase eq "G"){
		@varbases = ("A","C","T");
	}elsif($refbase eq "T"){
		@varbases = ("A","C","G");
	}else{
		print "Unknown base in ref sequence $refbase ... exiting\n";
		exit 2;
	}

	#for each possible base make a varstring and print id varbase refseq varseq	
	foreach (@varbases){
	
		my $var = $ref;
		substr($var,6,1) = $_;
	
		#output
		print "$id\t$_\t$ref\t$var\n";
	}

}
