#!/usr/bin/perl -w
use strict;
use warnings;

#merge multiple naked chromosome cuts together

my (@col, $line, %hash, %counthash);

while (@ARGV) {
    
    my $input = shift(@ARGV);
    
    #open file
    open(IN, $input) or die "Cant open(input file $input, $!\n";
    
        while (<IN>){
        
            chomp;
            @col = split (/\s+/, $_);
            my $kmer = shift @col;
            my $count = shift @col;
            
            if (exists $counthash{$kmer}) {
                
                $counthash{$kmer} += $count;                
                
            }else{
                $counthash{$kmer} = $count;
            }
            
            
            
            if (exists $hash{$kmer}) {
                my @temp = split(/\t+/, $hash{$kmer});#split up existsing hash
                my $length = @temp;
                $length--;
                for(my $i=0; $i<=$length; $i++){ #add values
                    
                    my ($ca, $cc, $cg, $ct) = split(":",$col[$i]);
                    my ($ha, $hc, $hg, $ht) = split(":",$temp[$i]);
                    
                    $ha += $ca;
                    $hc += $cc;
                    $hg += $cg;
                    $ht += $ct;
                    
                    #join together
                    $temp[$i] = "$ha:$hc:$hg:$ht";
                    
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
    print "$_\t".$counthash{$_}."\t".join("\t", $hash{$_})."\n";
}

