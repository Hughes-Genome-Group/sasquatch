#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::DB::Sam;


#reference genome

my $genome = $ARGV[0];

my $genome_file;

if($genome eq "hg18"){
	$genome_file = '/databank/igenomes/Homo_sapiens/UCSC/hg18/Sequence/WholeGenomeFasta/genome.fa';
}elsif($genome eq "hg19"){
	$genome_file = '/databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa';
}elsif($genome eq "mm9"){
	$genome_file = '/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa';
}else{
	print "Please select an implemented reference genome, hg19, hg18, mm9 !\n";
	exit 2;
}

#upload
my $fai = Bio::DB::Sam::Fai->load("$genome_file");

my $chromosome=$ARGV[1];
my $start=$ARGV[2];
my $stop=$ARGV[3];


my $fai_location = $chromosome.':'.$start.'-'.$stop;
my $refseq = $fai->fetch($fai_location);

        $refseq=~s/a/A/g;
        $refseq=~s/g/G/g;
        $refseq=~s/t/T/g;
        $refseq=~s/c/C/g;
        
print "Region:\t$fai_location : \n";
print "$refseq\n";
