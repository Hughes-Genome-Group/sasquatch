##############################################################################################
## Sasquatch, Sequence based predicting of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##############################################################################################
# Function: Data preprocessing; Create Strandspecific DNase I footprint wig tracks from aligned bam file
# Usage: samtools view ${BAM_FILE} | perl ${SCRIPT_DIR}/Print_It_All_combined.pl --build ${BIGWIG_CHRSIZES} --name ${IDTAG} --type ${SEQ_TYPE} -
# Author: Ron Schwessinger, Jim Hughes
# Date: 27/04/2016

#!/usr/bin/perl -w
use strict;
use Getopt::Long;

&GetOptions (
	     "window=i"=> \my $window_size,
	     "inc=i" => \my $window_incr,
	    "build=s" => \my $chr_lengths_file,
            "name=s" => \my $track_name,
	    "type=s" => \my $seq_type
);

open (OUTPUT1, ">$track_name\_Minus.wig");
open (OUTPUT2, ">$track_name\_Plus.wig");

unless($window_size) {$window_size = 1;} #changed form 3 to zero to get rid of windowing
unless($window_incr) {$window_incr = 1;}

my %chrlength;

open (SIZES, $chr_lengths_file);
	
	while (<SIZES>){
		chomp;
		my ($gchr, $gsize)= split(/\s+/);
		#print "$gchr $gsize\n";
		$gchr =~ s/chr//gi;
		unless ($gchr =~ /M/gi){   
		$chrlength{$gchr}=$gsize;
	}
	
}
	
close SIZES;

my %EndHash;
my %CombinedHashPlus;
my %CombinedHashMinus;
my @col;

if($seq_type eq "pairedend"){	#if singelend input

	while (<STDIN>) {
		#my $samcall=`samtools view -f 0x0042 $bamfile $col[0]:$col[1]-$col[2]`;
		my @splitcall=split(/\n+/,$_);
		foreach (@splitcall) {
			my $Orientation;	
			my ($readsam, $bitwisesam, $chrsam, $startsamm, $qscore, $match, $chr2, $pos2, $insert, @rest)=split(/\t+/,$_);
			unless ($bitwisesam){next;}
			unless  (($bitwisesam & 0x0002) && ($bitwisesam & 0x0040)) {next;}   #the read is mapped in a proper pair and the read is the first read in a pair
				$chrsam =~ s/chr//gi;
			if ($bitwisesam & 0x0002){
			    if ($insert > 0){$Orientation = 1}else{$Orientation = 2};
			    if ($match =~ /(\d+)M/){				
				my $matchSize = $1;
				my $SecondEnd = $pos2;		
				$SecondEnd = $SecondEnd + $matchSize - 1 if $Orientation == 1;	#minus 1 to map dnase cut on minus strand to the first read base (5' shift)
				$startsamm = $startsamm + $matchSize - 1 if $Orientation == 2; 

				$EndHash{$chrsam}{$startsamm}{$SecondEnd}{$Orientation}++;
	
			    }
			}
		}
	}
	
	foreach my $storedChr (sort keys %EndHash){
		foreach my $StoredE1 (sort by_number keys %{$EndHash{$storedChr}}){
				foreach my $StoredE2 (sort by_number keys %{$EndHash{$storedChr}{$StoredE1}}){
						foreach my $StoredORientation (sort by_number keys %{$EndHash{$storedChr}{$StoredE1}{$StoredE2}}){
			    
				if ($StoredORientation == 1) {	#first read plus mate reverse
				    
				    $CombinedHashPlus{$storedChr}{$StoredE1}++;
				    $CombinedHashMinus{$storedChr}{$StoredE2}++;
				    
				 }
				 if ($StoredORientation == 2) {	#first read minus mate plus
				 
				    $CombinedHashMinus{$storedChr}{$StoredE1}++;
				    $CombinedHashPlus{$storedChr}{$StoredE2}++;
				    
				  }
			  
				 }
			}
		}
	}

	
}elsif($seq_type eq "singleend"){	#if singelend input
	
	while(<STDIN>){
		chomp;
		my @splitcall=split(/\n+/,$_);
		foreach (@splitcall) {
			my $Orientation;	
			my ($readsam, $bitwisesam, $chrsam, $startsamm, $qscore, $match, @rest)=split(/\s+/,$_);
			$chrsam =~ s/chr//gi;
			if ($bitwisesam == 0){
				$Orientation = 1;
			}elsif($bitwisesam == 16){
				$Orientation = 2;
			}else{next;}
			if ($match =~ /(\d+)M/){
				my $matchSize = $1;
				$startsamm = $startsamm + $matchSize - 1 if $Orientation == 2;		#minus 1 to map dnase cut on minus strand to the first read base (5' shift)
				$EndHash{$chrsam}{$startsamm}{$Orientation}++;	
			}
		}
	}
	
	foreach my $storedChr (sort keys %EndHash){
    foreach my $StoredE1 (sort by_number keys %{$EndHash{$storedChr}}){
               foreach my $StoredORientation (sort by_number keys %{$EndHash{$storedChr}{$StoredE1}}){
	
	if ($StoredORientation == 1) {	#first read plus mate reverse
		
		$CombinedHashPlus{$storedChr}{$StoredE1}++;
		
	}
	if ($StoredORientation == 2) {	#first read minus mate plus
	
	        $CombinedHashMinus{$storedChr}{$StoredE1}++;
		
	}
   
        }
    }
}

	
}else{
	my $errorstring="no sequencing type defined: singleend or pairedend? exiting ...";
	system("echo $errorstring");
	exit 2;	
}

my %bins;

foreach my $ComBinedChr (sort keys %CombinedHashPlus){
    	print OUTPUT2 "variableStep  chrom=chr$ComBinedChr span=$window_incr\n";
    foreach my $combPosition(sort by_number keys %{$CombinedHashPlus{$ComBinedChr}}){
        
                my $total = $CombinedHashPlus{$ComBinedChr}{$combPosition};
                print OUTPUT2"$combPosition\t$total\n";
    }
}

foreach my $ComBinedChr2 (sort keys %CombinedHashMinus){
    	print OUTPUT1 "variableStep  chrom=chr$ComBinedChr2 span=$window_incr\n";
    foreach my $combPosition(sort by_number keys %{$CombinedHashMinus{$ComBinedChr2}}){
        
                my $total = $CombinedHashMinus{$ComBinedChr2}{$combPosition};
                print OUTPUT1"$combPosition\t$total\n";
    }
}

close(OUTPUT1);
close(OUTPUT2);


sub by_number {
	($a <=> $b);
	}
