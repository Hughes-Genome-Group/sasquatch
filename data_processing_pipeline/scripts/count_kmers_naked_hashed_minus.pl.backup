#!/usr/bin/perl -w
use strict;
use warnings;
#count kmers occurence in DNase peaks
use Bio::DB::Sam;

#Variables
my (@col,$line,$kmer,%kmers, @kmerstore, $temp, $start, $stop, %wighash,%kmersbarcode);

#INPUTs
my $kmer_in=$ARGV[0]; #input of all possible kmers
my $wig_file=$ARGV[1];

#reference genome
my $genome_file = $ARGV[4]; 
#'/databank/raw/hg18_full/hg18_full.fa';
#upload
my $fai = Bio::DB::Sam::Fai->load("$genome_file"); 

#additional parameters
my $length=$ARGV[2]; #kmer length
my $extend_kmersearch=0;   #bp length to extend the peaks up and downstream
my $extend_peak=150;    #bps to extend the peak for storage of the DNase (footprint) signal up and downstream

my $border_file=$ARGV[3];
my $output=$ARGV[5];

print "starting\n\n";

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

print "finished kmerprocessing...\n";

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

print "WIG store finished\n";

#open output file
open(OUT, ">$output") or die "Cant open output_file $output , $!\n";
#read in borders file 
open(BOR, $border_file) or die "Cant open border_file $border_file , $!\n";

    while (<BOR>) {
    
    chomp;
    my @bor=split(/\t+/, $_);

    print "processing $bor[0] ...\n";

#start at beginning of chromosome #chr21	48129895
    my $chromosome = $bor[0];
    
    my $catchstop = $bor[2]-150-$length;   #48129895 -150 -length
    
#loop through and store
for(my $catchstart=$bor[1]; $catchstart <= $catchstop; $catchstart++){
    
    my $stopcatching = $catchstart + $length - 1; #kmer length to get the sequence from -1 !!
        
    my $fai_location = $chromosome.':'.$catchstart.'-'.$stopcatching; #fetch and store underlying sequence of the reference genome
    my $subseq = $fai->fetch($fai_location);
    #cover small letters / substitute small letters to capital letters
    $subseq=~s/a/A/g;
    $subseq=~s/g/G/g;
    $subseq=~s/t/T/g;
    $subseq=~s/c/C/g;

    next if $subseq=~/[n,N]+/; #avoid uncertainties
    
    $kmer = $subseq; 
    if(exists $kmers{$kmer}){   #recent added
        $kmers{$kmer}++;
        #my $start = $catchstart-$extend_peak-1;		#why -1 ?
	my $start = $catchstart-$extend_peak;		#13.04. changed to fix +1bp shift bug discovered via dukenorm due to bed 0 indexing in peak region counting, -1 was introduced as quick & dirty fix     
        
        #get the footprint per region
        &footprint($chromosome,$start,$kmer);
    }

}#end catch loop

}

#sort keys
my @keys = sort { $kmers{$a} <=> $kmers{$b} } keys(%kmers);
@keys = reverse @keys;

foreach my $tempkey (@keys){
    #sort kmersbarcode stored hash entrys
    my $barcodeassemble;
    foreach (sort by_number keys %{$kmersbarcode{$tempkey}}){
        $barcodeassemble .= "\t".$kmersbarcode{$tempkey}{$_};
    }
    print OUT "$tempkey\t$kmers{$tempkey}".$barcodeassemble."\n";
}

close(BOR);
close(OUT);

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





