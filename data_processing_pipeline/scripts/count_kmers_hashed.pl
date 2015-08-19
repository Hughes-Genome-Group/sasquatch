#!/usr/bin/perl -w
use strict;
use warnings;
#count kmers occurence in DNase peaks
use Bio::DB::Sam;

#Variables
my (@col,$line,$kmer,%kmers, @kmerstore, $temp, $start, $stop, %wighash,%kmersbarcode);

#INPUTs
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


#read in peaks and coutn occurence;
open(IN, $peak_in) or die "Can't open peaks gff input file $peak_in $!";

    while (<IN>) {
        
        chomp;
        @col=split(/\t+/,$_);       #split peak and retrieve genomic coordinates
       
       my ($fai_location,$peakstart,$peakstop);
       my $chromosome = $col[0];

	next if $chromosome !~ /chr[0,1,2,3,4,5,6,7,8,9,X,Y]+/;
                
        if( ($peak_in =~ /\.bed$/) || ($peak_in =~ /\.narrowPeak$/) ) {#if bed format
            
            $peakstart = $col[1] - $extend_kmersearch + 1;       #extend for kmer search if wanted #+1 corrects for bed 0-based
            $peakstop = $col[2] + $extend_kmersearch;
            
        }elsif( $peak_in =~ /\.gff$/ ) {
            
            $peakstart = $col[3] - $extend_kmersearch;       #extend for kmer search if wanted
            $peakstop = $col[4] + $extend_kmersearch;
        
        }
         
        $fai_location = $chromosome.':'.$peakstart.'-'.$peakstop; #fetch and store underlying sequence of the reference genome
        my $subseq = $fai->fetch($fai_location);
        #cover small letters / substitute small letters to capital letters
        $subseq=~s/a/A/g;
        $subseq=~s/g/G/g;
        $subseq=~s/t/T/g;
        $subseq=~s/c/C/g;
        
        
        
        my $count=0;
        #go through sequence and count up kmers of length found in it
        while (length($subseq) >= $length) {
            
            $start=$peakstart+$count-$extend_peak; #get borders for searching the footprint (kmer length +- $extend_peak up and down stream)  #-1 norms for perl array index <. -1 removed 14.04 idx correct (perl array idx 0 correct initially)
            
            $kmer = substr $subseq, 0, $length; #get first $length (5,6) bases
            if(exists $kmers{$kmer}){   #recent added
                $kmers{$kmer}++;
            
                #get the footprint per region
                &footprint($chromosome,$start,$kmer);
            }
            $subseq= substr $subseq , 1; #remove first base of the substring
            
            $count++; #set count up as another base is proscced to adjust start and stop in the next loop
        }
        
    }

close(IN);

#sort keys
my @keys = sort { $kmers{$a} <=> $kmers{$b} } keys(%kmers);
@keys = reverse @keys;

foreach my $tempkey (@keys){
    #sort kmersbarcode stored hash entrys
    my $barcodeassemble;
    foreach (sort by_number keys %{$kmersbarcode{$tempkey}}){
        $barcodeassemble .= "\t".$kmersbarcode{$tempkey}{$_};
    }
    print "$tempkey\t$kmers{$tempkey}".$barcodeassemble."\n";
}


#Subroutine for creating the footprint of a peak region
sub footprint{
    
    $chrstore=shift @_;
    $start=shift @_;
    my $kmer_sub=shift @_;
    my $pos;
    
    #retrieve the footprint values at the specific kmer position in +- window and count kmerbarcoede hash up directly
    for(my $k=0; $k < $barcodelength; $k++){
        
        $pos = $start + $k;            
        $kmersbarcode{$kmer}{($k+1)} += $wighash{"$chrstore\t$pos"} if exists $wighash{"$chrstore\t$pos"};
      
    }
    
}


sub by_number {
	($a <=> $b);
}





