#get Arguments
args <- commandArgs(trailingOnly = TRUE)

source(args[3]) #source common functions from file

sequence <- args[1]	#kmer of interest

tosplit <- args[2]	#kmer lengths to retrieve

dissect.list <- dissect.sequence(sequence)

if(tosplit == 5){	#diseect into possible kmers
	toreturn <- dissect.list$list5
}else if(tosplit == 6){
	toreturn <- dissect.list$list6
}else if(tosplit == 7){
	toreturn <- dissect.list$list7
}else{
	print("Sequence to split to must be of legnth 5, 6 or 7!")
}

#return list
cat(toreturn, sep="\n")



