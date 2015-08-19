##merge profiles
#get Arguments
args <- commandArgs(trailingOnly = TRUE)

profile.plus=args[1]	

profile.minus=args[2]	

frag.type=args[3]	

kl=as.numeric(args[4])

#split up into profiles from unix stored string
profile.plus <- as.numeric(strsplit(profile.plus, ":")[[1]])
profile.minus <- as.numeric(strsplit(profile.minus, ":")[[1]])

profile.plus <- profile.plus/sum(profile.plus)
profile.minus <- profile.minus/sum(profile.minus)

if(frag.type == "ATAC"){
  
  profile.merge <- (profile.plus + profile.minus)	#sum up profile, average is calculated later anyway
  
}else if(frag.type == "DNase"){
  
  #	profile.merge <- (profile.plus + profile.minus)	#sum up profile, average is calculated later anyway
  
  #strand imbalance merge (half of each strand) and overlap in the kmer center
  ## ONLY FOR DNase !!! ##
  middle <- ( profile.plus[c(126:(125+kl))] + profile.minus[c(126:(125+kl))] ) / 2 #calc average of both strands for the kmer length middle
  
  profile.merge <- c( profile.plus[c(1:125)], middle, profile.minus[c((125+kl+1):length(profile.minus))] )	#combine left flank from plus right flank from minus and kmer middle from the average
  
}else{
  cat("Select according $FRAG_TYPE")
}

cat(profile.merge, sep=":")
