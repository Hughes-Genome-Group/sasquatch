
#get Arguments
args <- commandArgs(trailingOnly = TRUE)

source(args[3])	#source common functions

kmer=args[1]	#kmer of interest
infile=args[2]	#inout file to count from either tissue specific kmercount file or naked background count (plus or minus)

kmer.list <- decode.kmer(kmer)	#function that decodes ambivalent fasta code  letters like W into a list of all possible kmers ( if W two kmers to acutally search for)

kl=nchar(kmer)	#store kmer length

#initialize variable
count <- 0
profile <- rep(0,(250+kl))

for(km in kmer.list){	#for all entries in the splitted up kmer list

	match <- grep(km, readLines(infile), value=TRUE)	#performs basically a system grep / not sure how fast in relation to systemcommand
	
	split <- strsplit(match, "\t")[[1]]	#split on "\t" and unlist, all flat files are currently saved with a tab as seperator

	count <- count + as.numeric(split[2])	#second entry (first number) is the count how often the kmer appeared in the searched regions, rest is the cutting profile

	profile <- profile + as.numeric(split[c(28:(length(split)-25))])	#currently I trimm the profile for plotting purposes latter (-25 bp in the beginning and end)

}


#profile <- profile / sum(profile)	#normalize the counted cuts to relative cutting frequency in the investigated region

#return
cat(paste0("count=", count, ":::"))
cat(profile, sep=":")




