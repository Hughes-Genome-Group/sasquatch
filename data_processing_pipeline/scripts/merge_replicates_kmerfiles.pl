#!/usr/bin/perl -w
use strict;
use warnings;

#merge multiple replicate kmer profiles

my (@col, $line, %hash);
my $number = 0;

while (@ARGV) {
    
	$number++;

    my $input = shift(@ARGV);
    
    #open file
    open(IN, $input) or die "Cant open(input file $input, $!\n";
    
        while (<IN>){
        
            chomp;
            @col = split (/\s+/, $_);
            my $kmer = shift @col;
            
            if (exists $hash{$kmer}) {
                my @temp = split(/\t+/, $hash{$kmer});	#split up existsing hash
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

foreach my $k (@keys){

	#divide by number of input files
	my @temparray = split("\t", $hash{$k});
	
	print "$k";
	
	foreach (@temparray){
	
		my $divided = $_ / $number;
		print "\t$divided";

	}

	print "\n";

}

