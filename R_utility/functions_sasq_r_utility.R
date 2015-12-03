#####################################
# Functions for Sasquatch R utility #
#####################################

DecodeKmer <- function(kmer){
  # Decode ambiguous FASTA characters of kmer input
  #
  # Args:
  #   kmer: input kmer 
  #
  # Returns:
  #   list of non ambiguous kmers 
  
  mersplit <- unlist(strsplit(kmer,''))
  newlist <- c("")
  for(i in c(1:length(mersplit))){
    #simplebase > just push each newlist entry with base	
    if(mersplit[i] %in% c('A','C','G','T')){ newlist <- lapply(newlist,function(x){ x<-paste0(x,mersplit[i])}); }
    #N every base, repeat newlist 4 times add each base to 1/4 of newlist
    if(mersplit[i] == "N"){
      newlist<-rep(newlist,4)
      l<-length(newlist)
      newlist[c(1:(l/4))] <-paste0(newlist[c(1:(l/4))],"A")
      newlist[c(((l/4)+1):(l/2))] <-paste0(newlist[c(((l/4)+1):(l/2))],"T")
      newlist[c(((l/2)+1):((l/4)*3))] <-paste0(newlist[c(((l/2)+1):((l/4)*3))],"G")
      newlist[c((((l/4)*3)+1):l)] <-paste0(newlist[c((((l/4)*3)+1):l)],"C")
    }
    ###2er##	#K > G or T
    if(mersplit[i] == "K"){
      newlist<-rep(newlist,2)
      l<-length(newlist)
      newlist[c(1:(l/2))] <-paste0(newlist[c(1:(l/2))],"G")
      newlist[c(((l/2)+1):l)] <-paste0(newlist[c(((l/2)+1):l)],"T")
    }
    #M > A or C
    if(mersplit[i] == "M"){
      newlist<-rep(newlist,2)
      l<-length(newlist)
      newlist[c(1:(l/2))] <-paste0(newlist[c(1:(l/2))],"A")
      newlist[c(((l/2)+1):l)] <-paste0(newlist[c(((l/2)+1):l)],"C")
    }
    #R > A or G
    if(mersplit[i] == "R"){
      newlist<-rep(newlist,2)
      l<-length(newlist)
      newlist[c(1:(l/2))] <-paste0(newlist[c(1:(l/2))],"G")
      newlist[c(((l/2)+1):l)] <-paste0(newlist[c(((l/2)+1):l)],"A")
    }
    #Y > C or T
    if(mersplit[i] == "Y"){
      newlist<-rep(newlist,2)
      l<-length(newlist)
      newlist[c(1:(l/2))] <-paste0(newlist[c(1:(l/2))],"C")
      newlist[c(((l/2)+1):l)] <-paste0(newlist[c(((l/2)+1):l)],"T")
    }
    #S > C or G
    if(mersplit[i] == "S"){
      newlist<-rep(newlist,2)
      l<-length(newlist)
      newlist[c(1:(l/2))] <-paste0(newlist[c(1:(l/2))],"C")
      newlist[c(((l/2)+1):l)] <-paste0(newlist[c(((l/2)+1):l)],"G")
    }
    #W > A or T
    if(mersplit[i] == "W"){
      newlist<-rep(newlist,2)
      l<-length(newlist)
      newlist[c(1:(l/2))] <-paste0(newlist[c(1:(l/2))],"A")
      newlist[c(((l/2)+1):l)] <-paste0(newlist[c(((l/2)+1):l)],"T")
    }
    
    ###3er##	#B > C or G or T
    if(mersplit[i] == "B"){
      newlist<-rep(newlist,3)
      l<-length(newlist)
      newlist[c(1:(l/3))] <-paste0(newlist[c(1:(l/3))],"C")
      newlist[c(((l/3)+1):((l/3)*2))] <-paste0(newlist[c(((l/3)+1):((l/3)*2))],"G")
      newlist[c((((l/3)*2)+1):l)] <-paste0(newlist[c((((l/3)*2)+1):l)],"T")
    }
    #V > A or C or G
    if(mersplit[i] == "V"){
      newlist<-rep(newlist,3)
      l<-length(newlist)
      newlist[c(1:(l/3))] <-paste0(newlist[c(1:(l/3))],"A")
      newlist[c(((l/3)+1):((l/3)*2))] <-paste0(newlist[c(((l/3)+1):((l/3)*2))],"C")
      newlist[c((((l/3)*2)+1):l)] <-paste0(newlist[c((((l/3)*2)+1):l)],"G")
    }
    #H > A or C or T
    if(mersplit[i] == "H"){
      newlist<-rep(newlist,3)
      l<-length(newlist)
      newlist[c(1:(l/3))] <-paste0(newlist[c(1:(l/3))],"A")
      newlist[c(((l/3)+1):((l/3)*2))] <-paste0(newlist[c(((l/3)+1):((l/3)*2))],"C")
      newlist[c((((l/3)*2)+1):l)] <-paste0(newlist[c((((l/3)*2)+1):l)],"T")
    }
    #D > A or G or T
    if(mersplit[i] == "D"){
      newlist<-rep(newlist,3)
      l<-length(newlist)
      newlist[c(1:(l/3))] <-paste0(newlist[c(1:(l/3))],"A")
      newlist[c(((l/3)+1):((l/3)*2))] <-paste0(newlist[c(((l/3)+1):((l/3)*2))],"G")
      newlist[c((((l/3)*2)+1):l)] <-paste0(newlist[c((((l/3)*2)+1):l)],"T")
    }
  }
  
  newlist<-unlist(newlist)
  
  return(newlist)
  
}

GrepProfile <- function(kmer, infile){
  # grep 250 bp sourrounding kmer strand specific profile and normalize for 
  # total cuts in 250 bp window
  #
  # Args:
  #   kmer: input kmer 
  #   infile: input processed strandspecific kmer file
  #
  # Returns:
  #   normalized profile (relative cut probability) and kmer occurence count
  
  #decode  ambivalent fasta code char into kmer list
  kmer.list <- DecodeKmer(kmer)  
  kl <- nchar(kmer) #get kmer length  
  count <- 0 #initialise count variable
  profile <- rep(0,(250+kl)) #initialise empty profile
  
  for(km in kmer.list){  #for all entries in the splitted up kmer list
    match <- grep(km, readLines(infile), value=TRUE)  #system grep of kmer flat files
    split <- strsplit(match, "\t")[[1]]	#split on "\t" and unlist, as flat files are a tab seperated
    count <- count + as.numeric(split[2])	#get second column as the count
    split <- split[-c(1,2)] #remove kmer and count from split
    split <- as.numeric(split) #convert to numeric
    profile <- profile + split[c(26:(length(split)-25))]	#get 250 bp surrounding kmer from 300 bp profile
  }
  
  #normalize to all recorded cuts in the 250 surrounding region 
  #--> relative cut probability
  profile <- (profile / sum(profile))
  
  #make list for reporting
  newlist <- list( "profile"=profile, "count"=count ) #assemble return list
  return(newlist)
  
}

GetFootprint <- function(kmer, tissue, data.dir, length, frag.type, smooth){
  # Wrapper to retrieve merged, pruned, smoothed profile of kmer
  #
  # Args:
  #   kmer: input kmer 
  #   tissue: input tissue to query the profile from
  #   data.dir: repository storing processed kmer files per tissue
  #   length: length of window surrounding kmer to retrieve ("even 250 - 2 bps")
  #   frag.type: fragmentation type ("DNase" or "ATAC")
  #   smooth: flag if to smooth the profile (TRUE or FALSE)
  #
  # Returns:
  #   merged, pruned (smoothed) profile and count of the kmer occurence
  
  #get kmer length
  kl=nchar(kmer)
  
  #select flat strand specific files as inputs
  infile.plus=file.path(data.dir, tissue,"counts", paste0("kmers_", kl, "_count_", tissue, "_pnorm_JH60_plus.txt"))
  infile.minus=file.path(data.dir, tissue,"counts", paste0("kmers_", kl, "_count_", tissue, "_pnorm_JH60_minus.txt"))
  
  #grep strand specific profiles & counts
  l.plus <- GrepProfile(kmer, infile.plus)
  l.minus <- GrepProfile(kmer, infile.minus)
  
  #merge profiles accoring to selected fragmentation type
  if(frag.type == "ATAC"){
    #sum up profile as average of both profiles
    profile.merge <- (l.plus$profile + l.minus$profile)/2
  }else if(frag.type == "DNase"){
    #get length of the kmer as middle separately
    middle <- ( l.plus$profile[c(126:(125+kl))] + l.minus$profile[c(126:(125+kl))] ) / 2 #calc average of both strands for the kmer length middle
    #merge stradn specific with average of middle and only strand specific flanks
    profile.merge <- c( l.plus$profile[c(1:125)], middle, l.minus$profile[c((125+kl+1):length(l.minus$profile))] )	#combine left flank from plus right flank from minus and kmer middle from the average
  }else{
    cat("Select according framgentation type \"DNase\" or \"ATAC\"")
  }
  
  #smooth if specified so
  if(smooth){
    profile.merge <- ksmooth(c(1:length(profile.merge)), profile.merge, kernel="normal", bandwidth=5)$y
  }
  
  
  remove.temp <- (250-length)/2 #set positions to subtract from the profile
  
  #prune the profile to the desired length
  profile.merge <- profile.merge[c((1+remove.temp)):(length(profile.merge)-remove.temp)]
  
  return(list(profile=profile.merge, count=l.plus$count)) #return

}



