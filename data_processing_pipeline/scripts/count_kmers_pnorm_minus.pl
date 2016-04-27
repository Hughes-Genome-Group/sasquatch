##############################################################################################
## Sasquatch, Sequence based predicting of DNase I footprinting potential.                  ##
## Copyright (C) 2016 Genome Biology and Computational Biology Research Group, WIMM, Oxford ##
##############################################################################################
# Function: Data preprocessing; piling up kmer based profiles data minus strand
# Author: Ron Schwessinger
# Date: 27/04/2016

#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::DB::Sam;

#Variables
my (@col,$line,$kmer,%kmers, @kmerstore, $temp, $start, $stop, %wighash,%kmersbarcode);


##INPUTs
my $kmer_in=$ARGV[0]; #input of all possible kmers
my $peak_in=$ARGV[1];  #input of sample DNase peaks called
my $wig_file=$ARGV[2];

#additional parameters
my $length=$ARGV[3]; #kmer length
my $extend_kmersearch=0;   #bp length to extend the peaks up and downstream
my $extend_peak=150;    #bps to extend the peak for storage of the DNase (footprint) signal up and downstream

#reference genome
my $genome_file = $ARGV[4];
#upload
my $fai = Bio::DB::Sam::Fai->load("$genome_file");

#propensity norm factors file
my $propensity_in = $ARGV[5];

#read in propensity file andd store norm factor hash
my %norm_hash;
open(PROP, $propensity_in) or die "Can't open propensity file $propensity_in $!";

	my $disfirstline =<PROP>;	
	while(<PROP>) {
		@col=split(/\t+/,$_);
		$norm_hash{$col[0]} = sprintf("%.5f", $col[5]);	
	}

close(PROP);

###START###
my $barcodelength=(2*$extend_kmersearch)+(2*$extend_peak)+$length;
#read in and store kmers in hash, get length of interest (first kmer);
open(MERS, $kmer_in) or die "Can't open kmers input file $kmer_in $!";

    while (<MERS>) {

        chomp;
        
        if (length($_) != $length) {
            my $systemcall="Not all Kmers match Length: $length ! Check Input !\n";
            system("echo $systemcall");
            exit 2;
        }
      
        $kmers{$_}=0;
        #create second hash store for the barcodes to count up per kmer
        for(my $j=1; $j <= $barcodelength; $j++){   #initialize a barcode hash
            $kmersbarcode{$_}{$j} = 0;
        }
    }
close(MERS);


my $chrstore;
#read in wig file for footprints store in wighash
open(WIG, $wig_file) or die "Can't open wig input file $wig_file $!";
    while (<WIG>) {

        chomp;
        if($_=~/variableStep/) {
            $_=~/chrom=(chr\w+)\s+/;
            $chrstore=$1;
            next;
        }
        @col=split(/\t+/,$_);
        $wighash{"$chrstore\t$col[0]"}=$col[1];
    }    
close(WIG);


##read in peaks and count occurence;
open(IN, $peak_in) or die "Can't open peaks gff input file $peak_in $!";

    while (<IN>) {
        
        chomp;
        @col=split(/\t+/,$_);       #split peak and retrieve genomic coordinates
       
       my ($fai_location,$peakstart,$peakstop);
       my $chromosome = $col[0];

       	next if $chromosome !~ /chr[0,1,2,3,4,5,6,7,8,9,X,Y]+/;

        if( ($peak_in =~ /\.bed$/) || ($peak_in =~ /\.narrowPeak$/) ) {#if bed format
            
            $peakstart = $col[1] - $extend_kmersearch - $extend_peak - 2 + 1;       #extend for kmer search if wanted +1 for bed 0-based
            $peakstop = $col[2] + $extend_kmersearch + $extend_peak + 3;
            
        }elsif( $peak_in =~ /\.gff$/ ) {
            
            $peakstart = $col[3] - $extend_kmersearch - $extend_peak - 2;       #extend for kmer search if wanted
            $peakstop = $col[4] + $extend_kmersearch + $extend_peak + 3;	#extend with 3 bp to alwas complete the 6mer for bias norm
        
        }
         
        $fai_location = $chromosome.':'.$peakstart.'-'.$peakstop; #fetch and store underlying sequence of the reference genome
        my $subseq = $fai->fetch($fai_location);
        #cover small letters / substitute small letters to capital letters
        $subseq=~s/a/A/g;
        $subseq=~s/g/G/g;
        $subseq=~s/t/T/g;
        $subseq=~s/c/C/g;       

	my $startkmerssave = $extend_peak + 2 - 1; #where to start saving the kmers (actual peak is extended by 150 bp +- 3bp for completing the 6mers for duke norm so 150 + 3 and - 1 for perl 0 based indexing

        my $count=0;
        #go through sequence and count up kmers of length found in it
        while (length($subseq) >= ($length + ($extend_peak * 2) + 5) ) {		#until length(6) + 153 flanking
            
            $start=$peakstart+$count+2; 
		#get borders for searching the footprint (kmer length +- $extend_peak up and down stream)
            
            $kmer = substr $subseq, $startkmerssave, $length; #get first $length (5,6) bases
	#prune subseq for sub routine fetching
	my $pruned_subseq = substr $subseq, 0, ($length + ($extend_peak * 2) + 5 );

            if(exists $kmers{$kmer}){   #recent added
                $kmers{$kmer}++;
            
                #get the footprint per region
                &footprint_dukenorm($chromosome, $start ,$kmer, $pruned_subseq);
            }
            $subseq= substr $subseq , 1; #remove first base of the substring
            
            $count++; #set count up as another base is proscced to adjust start and stop in the next loop
        }
        
    }

close(IN);

###sort keys
my @keys = sort { $kmers{$a} <=> $kmers{$b} } keys(%kmers);
@keys = reverse @keys;

foreach my $tempkey (@keys){
    #sort kmersbarcode stored hash entrys
    my $barcodeassemble;
    foreach (sort by_number keys %{$kmersbarcode{$tempkey}}){
 	$barcodeassemble .= "\t".$kmersbarcode{$tempkey}{$_};	##check output       
#$barcodeassemble .= "\t".sprintf("%e" ,$kmersbarcode{$tempkey}{$_});	##check output
    }
    print "$tempkey\t$kmers{$tempkey}".$barcodeassemble."\n";
}


#Subroutine for creating the footprint of a peak region
sub footprint_dukenorm{
    
	$chrstore=shift @_;
	$start=shift @_;
	my $kmer_sub=shift @_;
	my $prunedseq=shift @_;
	my ($pos, $mer, $normf);
    
    #retrieve the footprint values at the specific kmer position in +- window and count kmerbarcoede hash up directly
    for(my $k=0; $k < $barcodelength; $k++){
        
	$mer = substr $prunedseq, $k, 6;	#get underlying 6mer
	
	next if $mer=~/N/;	#test wise skip N containing kmers

	if(exists $norm_hash{$mer}) { 
		$normf = $norm_hash{$mer};
	}else{	
		next;
	}
        
	$pos = $start + $k;    
	#normalize wig entry trimmed to 2 decimal digits
        $kmersbarcode{$kmer_sub}{($k+1)} += sprintf("%.5f", ($wighash{"$chrstore\t$pos"} * $normf) ) if exists $wighash{"$chrstore\t$pos"};
      
    }
    
}


sub by_number {
	($a <=> $b);
}





