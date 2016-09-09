#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::DB::Sam;
use Getopt::Long;
use Pod::Usage;

# 2016/08/11 
# Ron Schwessinger


=head1 USAGE
  
	helper_get_ref_and_var_sequence.pl --bed <bed.like file> --genome <genome.fa> --ref column_of_ref_base --var column_of_var_base --extend sequence_length_to_extract

	--bed  bed like file (0-coordinate based) with at least 5 columns [chr start end ref_base var_base] more columns are allowed and the position of the ref and var base columns can be specified. Expects a 1 bp length input == exact (start)position of the SNP/INDEL
	--genome  genome.fa fasta file of reference genome of interest (can be a short cut like hg19, hg18, mm9 if the appropriate paths have been set within the script)
	--ref  column number (1 based) where to find the reference base(s) 
	--var  column number (1 based) where to find the variant base(s)
	--extend  bp length which sequence length surroudning the variant position should be retrieved [default = 13] for SasQ input. Should be an odd number to ensure centric SNP variant position
	--warnings  warnings.log file to record unmatching sequences [default ./warnings.log]
        
=head1 DESCRIPTION

   Take a bed like file input (chr start stop ...) with variable additional column (content and number). Extract the reference sequence surrounding a region of interest up to an exntending length of interest [default=13]. Expects a 1bp length bed region for the SNP inout. Check and exchange the center position to be the reference base and create a reference and variant sequence from the provided columns of reference and variant base.
        
=cut

# =============================================================================
# ===== GET ARGUMENTS AND SET PARAMETERS ======================================
# =============================================================================

# define default options
my $extend=13;  # target sequence length to extend to [default=13]
my ($genome, $bed, $ref_col, $var_col) = "";
my $warnings = "warnigns.log";

# Get variables to from command line -----------------------------------------
 GetOptions(
	q(help) => \my $help,   # get help
	q(man) => \my $man,     # print description
	"genome=s" => \$genome,       # genome.fa file or short cut if set up in perl script
        "bed=s"   => \$bed,      # bed like file input
        "ref=i"  => \$ref_col,   # column holding the reference base(s)
        "var=i" => \$var_col,      # colum holding the variant basis
        "extend=i" => \$extend,      # target sequence length to extend to [default=13]
	"warnings=s" => \$warnings	# write warings and skipped variants to
        ) or pod2usage(q(-verbose) => 1);

# Help ------------------------------------------------------------------------
pod2usage(1) if ($help);
pod2usage(q(-verbose) => 2) if ($man);
pod2usage(1) if(($genome eq "" || $bed eq "" || $ref_col eq "" || $var_col eq ""));
                                                           

# HELP MESSAGE ----------------------------------------------------------------
if($help){
	print "\nhelper_get_ref_and_var_sequence.pl\nUsage: perl trim_bed_chrom_sizes.pl --bed <bedfile> --sizes <chrom.sizes.file> \n\n";
	exit;
}


# check extension length
if($extend % 2 == 0){
	print "Please specifiy and odd length of sequence window to extract to ensure centric position of the variant!\n";
	exit 2;
}

# select appropriate reference genome -----------------------------------------
my $genome_file = "";

#### DEFINE DEFAULT ROUTES TO GENOME.FA FILES HERE FOR SHORTCUT USAGE <--------
if($genome eq "hg18"){
	$genome_file = '/databank/igenomes/Homo_sapiens/UCSC/hg18/Sequence/WholeGenomeFasta/genome.fa';
}elsif($genome eq "hg19"){
	$genome_file = '/databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa';
}elsif($genome eq "mm9"){
	$genome_file = '/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa';
}

# check if proper genome file selected
if(!-s $genome_file){
	print "Please select an implemented reference genome fasta file!\n or a defined short cut if paths have been specified within the helper script: e.g. hg19, hg18, mm9 !\n";
	exit 2;
}

# load connection to genome file
my $fai = Bio::DB::Sam::Fai->load("$genome_file");

# Adjust ref and var columns for perl 0-based indices
$ref_col--;
$var_col--;

# =============================================================================
# # ===== START ===============================================================
# # ===========================================================================
my $linecounter = 0;

# open warnings file
open(WARN, ">$warnings") or die "Can not open warnings log file $warnings $!";

# Get and CHECK input format ----------------------------------------------------------
open(IN, $bed) or die "Can not open input bed like file $bed ! $!\n";
	
	while(<IN>){

		chomp;
		$linecounter++;

		my @col = split(/\t+/, $_);

		# check column number
		if(scalar(@col)  <= $ref_col || scalar(@col) <= $var_col){
			print "Input file has less columns then required at line $linecounter!\n";
			exit 2;
		}

		my $chromosome = $col[0];
		my $start = $col[1];	
		my $end = $col[2];
	
		# check requested region length
		my $region_length = $end-$start;
                if($region_length > 1){
                        print WARN "#warning: Requested poistion is longer than 1 bp. Use 1 bp position of SNP or start coordinate of INDEL! line $linecounter ... skipping\n";
			print WARN "$_\n";
                        next;
                }

		my $ref_base = $col[$ref_col];
		my $var_base = $col[$var_col];

		# check length and content of ref and variant base
		if($ref_base=~/[^A,C,G,T,N,R,Y,K,M,S,W,B,D,H,V]+/ | $var_base=~/[^A,C,G,T,N,R,Y,K,M,S,W,B,D,H,V]+/){
			print WARN "#warning: Reference or variant have not allowed characters! Only standard fasta characters are allowed! line $linecounter .. skipping\n";
			print WARN "$_\n";
			next;
		}
		# get variant type snp or indel
		my $var_type = "snp";
		$var_type = "indel" if (length($ref_base)>1 || length($var_base)>1);

		
		# Extract sequence surrounding position of interest
		$start++; # adjust for bed 0-based coordinate

		# get extension length and adjust coordinates (take into account length of reference base .. for indels)	
		$start -= ($extend-1)/2;
		$end += ($extend-1)/2 + (length($ref_base)-1);

		# set region to extract
		my $fai_location = $chromosome.':'.$start.'-'.$end;

		# extract
		my $refseq = $fai->fetch($fai_location);

		# clean sequence
	        $refseq=~s/a/A/g;
        	$refseq=~s/g/G/g;
	        $refseq=~s/t/T/g;
        	$refseq=~s/c/C/g;
       

		# Check Ref and make Var sequence
		my $varseq = $refseq;

		# check ref at centric position
		if(substr($refseq, (($extend-1)/2), length($ref_base)) ne $ref_base){
			print WARN "#warning: Ref base does not match center of extracted sequence: $ref_base in $refseq line $linecounter complementary? ... skipping \n";
			print WARN "$_\n";
			next;
		}

		# if matches exchange ref to var base(s)
		substr($varseq, (($extend-1)/2), length($ref_base)) = $var_base;


		#print final line 
		print join("\t", @col)."\t$refseq\t$varseq\n";


	}

close(IN);
close(WARN);
