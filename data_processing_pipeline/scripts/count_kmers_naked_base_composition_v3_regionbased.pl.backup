##!/usr/bin/perl -w
use strict;
use warnings;
use Bio::DB::Sam;

#Variables
my (@col,$line,$kmer,%kmers, @kmerstore, $temp, $start, $stop, @barcode, @tempsplit, $composition, $dropoffkmerbase);
my @streams = qw(ds us);
my @bases=qw(A C G T);
my @basesN=qw(A C G T N);
#INPUTs
my $kmer_in=$ARGV[0]; #input of all possible kmers
##reference genome
#my $genome_file = '/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa';
my $genome_file = '/databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa';
my $fai = Bio::DB::Sam::Fai->load("$genome_file");


##additional parameters
my $length=$ARGV[1]; #kmer length
my $extend_kmersearch=0;   #bp length to extend the peaks up and downstream
my $extend_peak=150;    #bps to extend the peak for storage of the DNase (footprint) signal up and downstream
my $regions_file=$ARGV[2];
my $output=$ARGV[3];

print "starting\n";

#read in and store kmers in hash, check length of interest (kmers);
open(MERS, $kmer_in) or die "Can't open kmers input file $kmer_in $!";
    while (<MERS>) {
        chomp;      
        if (length($_) != $length) {
            my $systemcall="Not all Kmers match Length: $length ! Check Input !\n";
            system("echo $systemcall"); #exit if one length doesn't fit
            exit 2;
        }
        $kmers{$_}=0;
    }
    $kmers{"NNN"}=0;    #dummy fo ambivalence
close(MERS);
#create second hash store for the combcodes to count up base composition
my %combstore;
my $combcodelength=(2*$extend_kmersearch)+(2*$extend_peak)+$length;

foreach my $k (keys %kmers){    #initialize the base counting uo for all possible kmers
    for(my $i=1; $i <= $extend_peak; $i++){
        $combstore{$k}{"ds"} {$i} {"A"} = 0;
        $combstore{$k}{"ds"} {$i} {"C"} = 0;
        $combstore{$k}{"ds"} {$i} {"G"} = 0;
        $combstore{$k}{"ds"} {$i} {"T"} = 0;
        $combstore{$k}{"ds"} {$i} {"N"} = 0;    #ambivalence dummy
        $combstore{$k}{"us"} {$i} {"A"} = 0;
        $combstore{$k}{"us"} {$i} {"C"} = 0;
        $combstore{$k}{"us"} {$i} {"G"} = 0;
        $combstore{$k}{"us"} {$i} {"T"} = 0;
        $combstore{$k}{"us"} {$i} {"N"} = 0;    #ambivalence dummy
    }
    $combstore {$k} {"count"} = 0;  #count of kmer occurence store
}
print "finished initializing kmer hashes\n";


#read in regions_file
my %regionstore;
open(IN, $regions_file) or die "Can't open regions gff/bed input file $regions_file $!";

    while (<IN>) {

        chomp;
        @col=split(/\t+/,$_);       #split peak and retrieve genomic coordinates
        
        if( $regions_file =~ /\.bed$/ ) {#if bed format 
		
	$col[1]+=1;	#correct for bed 0-based
            $regionstore{"$col[0]\t$col[1]\t$col[2]"}=1;    #store border of region in hash key

        }elsif( $regions_file =~ /\.gff$/ ) {#if bed format 

            $regionstore{"$col[0]\t$col[3]\t$col[3]"}=1;    #store border of region in hash key

        }else{
		print "please use .bed or .gff as regionsfileformat and in the file name\n";
		exit 2;	
	}
    }
close(IN);
#count total number of regions to process
my @regions_array = keys(%regionstore);
my $total_regions = @regions_array;
#$total_regions = $total_regions/10;

my $regionscount = 0;
###actual loop start###
print "\nStarting with Sequence processing ...\n";
system("date");

foreach (@regions_array){

    $regionscount++;
    
    #if ( ($regionscount % $total_regions) == 0 ) {
                
        #        my $tempstepped = int($regionscount / $total_regions);
       #         my $systemstring = "$tempstepped / 10 in progress ...\n"; #print estimated progress
      #          system("echo $systemstring");
                
     #           &print_out; #print temporary output file
    #}
    
    #split region and get chr and m borders
    my ($chromosome, $border1, $border2) = split(/\t+/, $_);
    
    #upload chromosome
    my $fai_location = $chromosome.":".$border1."-".$border2; #fetch and store underlying sequence of the reference genome
    my $chrseq = $fai->fetch($fai_location);
    
    #estimate length to process and set timepoints for temp output
    my $toprocess = $border2 - $border1;
    
    # Search for first ambivalence free substring match found in chrseq
    next if $chrseq !~ m/([A,C,G,T,a,c,g,t]{$combcodelength})/; #break out of region
    #$chrseq =~ m/([A,C,G,T,a,c,g,t]{$combcodelength})/;
    my $seq_to_initialize = $1;
    
    
    
    my $offset = $+[0];     #set offset for running kmer processing
    my $searchstop = $border2 - 150 - $border1; #set stoppoint
    
    #store kmer and up and downstream sequence in temp store object
    $kmer = substr($seq_to_initialize, $extend_peak, $length);  
    $kmer=~s/a/A/g;
    $kmer=~s/g/G/g;
    $kmer=~s/t/T/g;
    $kmer=~s/c/C/g;
    
    my $tmp_ds = substr($seq_to_initialize, ($extend_peak + $length), $extend_peak);
    my $tmp_us = substr($seq_to_initialize, 0, $extend_peak);
    
    my (%usbases, %dsbases); #init hashes to temporary store bases at each us and ds position
    
    #count up base occurene for initial sequence region (also initialize temp base store hash keys and values )
    #upstream
    while($tmp_ds =~ m/A{1}/gi){ $combstore {'ds'} {$+[0]} {'A'} ++; $dsbases {$+[0]} = "A"; }  #end of single letter match $+[1] corresponds 
    while($tmp_ds =~ m/C{1}/gi){ $combstore {'ds'} {$+[0]} {'C'} ++; $dsbases {$+[0]} = "C"; }  #to the 1 start index based position of the letter  
    while($tmp_ds =~ m/G{1}/gi){ $combstore {'ds'} {$+[0]} {'G'} ++; $dsbases {$+[0]} = "G"; }
    while($tmp_ds =~ m/T{1}/gi){ $combstore {'ds'} {$+[0]} {'T'} ++; $dsbases {$+[0]} = "T"; }
    #downstream
    while($tmp_us =~ m/A{1}/gi){ $temp = $extend_peak - $-[0]; $combstore {'us'} {$temp} {'A'} ++; $usbases {$temp} = "A"; }   #convert upstream idxs 1 > 150, 2 > 149, ...  
    while($tmp_us =~ m/C{1}/gi){ $temp = $extend_peak - $-[0]; $combstore {'us'} {$temp} {'C'} ++; $usbases {$temp} = "C"; }   # 150 - start of match = distance to kmer
    while($tmp_us =~ m/G{1}/gi){ $temp = $extend_peak - $-[0]; $combstore {'us'} {$temp} {'G'} ++; $usbases {$temp} = "G"; }
    while($tmp_us =~ m/T{1}/gi){ $temp = $extend_peak - $-[0]; $combstore {'us'} {$temp} {'T'} ++; $usbases {$temp} = "T"; }
    
    $combstore {$kmer} {"count"} ++;    #count occurence up
    
    my $updatekmer; #init parallel kmer sotrage variable for writing ambivalent kmers to dummy hsah store
    
    # move along sequence and change 
    for(my $p = $offset; $p <= $searchstop; $p++) {       #go along accessible sequence starting from offset of first decent match
            
            
            #get next base
            my $nextbase = substr($chrseq, $p, 1);
            $nextbase=~s/a/A/g;
            $nextbase=~s/g/G/g;
            $nextbase=~s/t/T/g;
            $nextbase=~s/c/C/g;
            if ($nextbase =~ /[^A,C,G,T]/i ) {      #treat as ambivalent
                $nextbase = "N";
            }
    
            #update kmer
            $dropoffkmerbase = substr($kmer,0,1);   #store for counting up us base
            $kmer = substr($kmer, 1, ($length - 1)); 
            $kmer .= $dsbases{1};
            
            if( $kmer =~ /N/ ){ $updatekmer = "NNN"; } else{ $updatekmer = $kmer; } #if ambivalent throw into dummy store
            
            $combstore {$kmer} {"count"} ++;    #count up occurence
            
            #counting up upstream position bases and shift the bases from the temporary base store hash( idx new = idx old - 1)
            for(my $i = ($extend_peak - 1); $i >= 2; $i--){
                $usbases{$i} = $usbases{($i-1)};
                $combstore {$updatekmer} {"us"} {($i)} { $usbases {($i)} } ++;            
            }
            #count up upstream position 1
            $combstore {$updatekmer} {"us"} {1} {$dropoffkmerbase} ++;
            $usbases{1} = $dropoffkmerbase; #new upstream position 1 is the first character of the former kmer
            
             #counting up downstream position bases and shift the bases from the temporary base store hash( idx new = idx old + 1)
            for(my $i = 1; $i <= ($extend_peak - 1); $i++){
                $combstore {$updatekmer} {"ds"} {($i)} { $dsbases{($i+1)} } ++;
                $dsbases{$i} = $dsbases{($i+1)};
            }
            #count up downstream position 1
            $combstore {$updatekmer} {"ds"} {1} {$nextbase} ++;
            $dsbases{$extend_peak} = $nextbase; #new downstream 150 is the retrieved nextbase
     
    }

}
print "finished with sequence processing printing output\n";
system("date");

#5##PRINT STORED and SORTED COMPOSITION
&print_out;


#subroutine for output
sub print_out{
    #open output file
    open(OUT, ">$output") or die "Cant open output_file $output , $!\n";
    ###output function##
    foreach my $km (sort keys %kmers){
        
        next if $km eq "NNN";
        
        #assemble composition for identity kmer count in the middle
        my $kmercount = $combstore{$km}{'count'};
        my @identity;
        for(my $i = 0; $i < $length; $i++){
            $temp = substr($km,$i,1);
            if ($temp eq "A") { push(@identity, "$kmercount:0:0:0");
            }elsif($temp eq "C") { push(@identity, "0:$kmercount:0:0");
            }elsif($temp eq "G") { push(@identity, "0:0:$kmercount:0");
            }elsif($temp eq "T") { push(@identity, "0:0:0:$kmercount");
            }        
        }
        my $id_string = join("\t",@identity);
        my @ds_string_store;
        my @us_string_store;
        
        foreach my $po (sort by_number keys $combstore{$km}{"ds"}){
       
            #assemble base counts in a single sting ':' seperated
            $temp = $combstore{$km}{"ds"}{$po}{"A"}.':'.$combstore{$km}{"ds"}{$po}{"C"}.':';
            $temp .= $combstore{$km}{"ds"}{$po}{"G"}.':'.$combstore{$km}{"ds"}{$po}{"T"};
            push(@ds_string_store, $temp);      
        }
        
        foreach my $po (reverse sort by_number keys $combstore{$km}{"us"}){
            #assemble base counts in a single sting ':' seperated
            $temp = $combstore{$km}{"us"}{$po}{"A"}.':'.$combstore{$km}{"us"}{$po}{"C"}.':';
            $temp .= $combstore{$km}{"us"}{$po}{"G"}.':'.$combstore{$km}{"us"}{$po}{"T"};
            push(@us_string_store, $temp);        
        }
    
    #assemble for output: us pos, kmeroccurence, ds pos
        my $output_string = $km."\t".$combstore{$km}{"count"}."\t";
        $output_string .= join("\t", @us_string_store)."\t"."$id_string\t".join("\t", @ds_string_store);
        
        print OUT "$output_string\n";
    
    }
    close(OUT);
}

#mini sub number sorting
sub by_number {
    ($a <=> $b);
}


