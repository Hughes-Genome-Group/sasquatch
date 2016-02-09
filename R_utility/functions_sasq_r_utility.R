############################################
# Functions& Theme for Sasquatch R utility #
############################################

library(ggplot2)
# library(gridExtra)
library(grid)
library(RColorBrewer)

### THEMES ####
#define theme for plotting
science_theme <- theme(
  panel.grid.major = element_line(size = 0.5, color = "grey"),
  axis.line = element_line(size = 0.7, color = "black"),
  text = element_text(size = 14)
)

### HELPER FUNCTIONS ####
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

Sobeln <- function(profile){
  # Calculate 1st derivative approximation of profile by 1D sobel filtering
  #
  # Args:
  #   profile: input profile
  #
  # Returns:
  #   1st derivative approximation of profile 
  
  b=rep(0,length(profile))
  for(i in c(2:(length(profile)-1))){
    b[i] = ( profile[i-1] * -1 ) + ( profile[i] * 0 ) + ( profile[i+1] * 1 )
  }
  b[1] = b[2]
  b[length(profile)] <- b[length(profile)-1]
  return(b)
}


### BASIC FUNCTIONS ####
DissectSequence <- function(sequence, kl, list=TRUE){
  # split longer Sequence into list of kmers
  #
  # Args:
  #   kmer: input kmer
  #   kl: length of kmer to split the sequence into 
  #   (currently 5, 6 or 7 is supported for subsequent functions)
  #   list: TRUE/FALSE idnicate if to return a list
  #
  # Returns:
  #   splitted list of k-mers of length kl 
  
  #capture if kl is not 5,6 or 7
  if(!kl %in% c(5,6,7)){
    warning("kl has to be 5, 6 or 7! Returning NA!")
    return("NA")
  }
  
    dissectlist <- list()
    
    for( i in c(1:(nchar(sequence)-kl+1))){	#split in all possible kl-mers
     
       dissectlist <- c(dissectlist, substr(sequence, i, (i+kl-1)))
    
    }

    if(list){    
      return(dissectlist)
    }else if(!list){
      return(unlist(dissectlist))
    }else{
      warning("Please specifiy if to return a list or an array: list TRUE/FLASE")
      return(NA_real_)
    }
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

CalcSFR <- function(profile, 
                    us.mid, 
                    ds.mid, 
                    range.us, 
                    range.ds){
  # Calculate Shoulder to Footprint Ratio of a footprint profile
  #
  # Args:
  #   profile: input profile (should be smoothed)
  #   us.mid: middle point of upstream shoulder
  #   ds.mid: middle point of downstream shoulder
  #   range.us: range of upstream shoulder
  #   range.ds: range of downstream shoulder
  #
  # Returns:
  #   SFR
  
  #window <- 6
  ext.us <- range.us/2
  ext.ds <- range.ds/2
  
  #get the estimated shoulders 
  shoulder <- profile[ c( c((us.mid-ext.us):(us.mid+ext.us)), c((ds.mid-ext.ds):(ds.mid+ext.ds)))]
  shoulder.cuts <- sum(shoulder)	#get all cuts in the shoulders
  footprint <- profile[c((us.mid+ext.us+1):(ds.mid-ext.ds-1))]	#get all cuts in the footprint
  footprint.cuts <- sum(footprint)
  
  #Calculate the SFRatio
  sfr <- (shoulder.cuts/length(shoulder))/((footprint.cuts)/length(footprint))
  
  return(sfr)
  
}

SobelBorders <- function(profile, kl){
  # Estimate the footprint shoulders by based on zero crossings of the 1D
  # 1st derivative approximation of the smoothe profile and  get the optimal 
  # shoulder range by estimating the maximizing the SFRatio
  #
  # Args:
  #   profile: smoothed profile 
  #   kl: length of the centeric kmer
  #
  # Returns:
  #   middle points and bp ranges of upstream and downstream shoulder border 
  #   and flag if borders could be estimated
  
  #Calculate the approximated 1st derivative via 1D discrete sobel operator
  profile.sobel <- Sobeln(profile) 
  
  #get 0 crossings with correct directionality (maxima)
  zero.cross = which(diff(sign(profile.sobel)) == -2)
  
  #separate crossings into upstream and downstream
  us.cross <- subset(zero.cross, zero.cross <= (126))	
  ds.cross <- subset(zero.cross, zero.cross >= (126 + kl - 1))
  
  #only consider crossings in a reasonable window to speed up
  search.range <- 50 #max bp distance to consider zero crossings to speed up search
  
  us.cross <- subset(us.cross, us.cross >= 126 - search.range)
  ds.cross <- subset(ds.cross, ds.cross <= (125 + kl + search.range))
  
  ### === first round with fixed shoulder range to get optimal maximas === ###
  peak.range <- c(6) #select fixed range to speed up first round
  cross.combs <- as.matrix(expand.grid(us.cross, ds.cross, peak.range, peak.range))	#get a matrix with all combinations of the relevant crossings
  #cover cases with no relevant maxima 
  if(dim(cross.combs)[1] == 0) {
    newlist <- list("us"=0, "ds"=0, "range.us"=0, "range.ds"=0, "flag"=FALSE)
    return(newlist)
  }
  
  #calculate SFR for the range of shoulder widths and extract the max one from each
  for(j in c(1:dim(cross.combs)[1])){
    sfr.store <- apply(cross.combs, 1, function(x){ 
      u=x[1]
      d=x[2]
      ru=x[3]
      rd=x[4]
      CalcSFR(profile, u, d, ru, rd)
    })
  }
  
  #get the borders
  us.border <- as.numeric(cross.combs[which(sfr.store == max(sfr.store)),][1])
  ds.border <- as.numeric(cross.combs[which(sfr.store == max(sfr.store)),][2])
  
  #consider valid borders not in kmer range and not overlapping borders
  #as long as the max SFR yielding borders yield invalid borders remove the maximum and get to next
  while( ((us.border + (2)) >= (ds.border - (2))) || ( (us.border ) > 126 ) || ( (ds.border ) < (126 + kl - 1) ) ){
    if(length(sfr.store) <= 1) {	#capture if no suitable borders found
      newlist <- list("us"=0, "ds"=0, "range.us"=0, "range.ds"=0, "flag"=FALSE)
      return(newlist)
    }
    
    sfr.store <- sfr.store[-which(sfr.store == max(sfr.store))] #remove maximum and select next best borders
    us.border <- as.numeric(cross.combs[which(sfr.store == max(sfr.store)), ][1])
    ds.border <- as.numeric(cross.combs[which(sfr.store == max(sfr.store)), ][2])
  }
  
  ### === second round to get optimal range	=== ###
  peak.range <- c(4,6,8,10) #restrict valid border ranges to speed up
  cross.combs <- as.matrix(expand.grid(us.border, ds.border, peak.range, peak.range))
  
  #calculate SFR for the range of shoulder widths and extract the max yielding
  for(j in c(1:dim(cross.combs)[1])){
    sfr.store <- apply(cross.combs, 1, function(x){ 
      u=x[1]
      d=x[2]
      ru=x[3]
      rd=x[4]
      CalcSFR(profile, u, d, ru, rd)
    })
  }
  
  sfr.maxima <- max(sfr.store)
  us.range <- as.numeric(cross.combs[which(sfr.store == max(sfr.store)),][3])
  ds.range <- as.numeric(cross.combs[which(sfr.store == max(sfr.store)),][4])
  
  #consider valid borders not in kmer range and not overlapping borders
  #as long as the max SFR yielding borders yield invalid borders remove the maximum and get to next
  while( ((us.border + (us.range/2)) >= (ds.border - (ds.range/2))) || ( (us.border + (us.range/2)) > (126 + 2 ) ) || ( (ds.border - (ds.range/2)) < (126 + kl - 1 - 2) ) ){
    if(length(sfr.store) <= 1) {	#no suitable borders found
      newlist <- list("us"=0, "ds"=0, "range.us"=0, "range.ds"=0, "flag"=FALSE)
      return(newlist)
    }
    sfr.store <- sfr.store[-which(sfr.store == max(sfr.store))]
    
    us.border <- as.numeric(cross.combs[which(sfr.store == max(sfr.store)),][1])
    ds.border <- as.numeric(cross.combs[which(sfr.store == max(sfr.store)),][2])
    sfr.maxima <- max(sfr.store)
    us.range <- as.numeric(cross.combs[which(sfr.store == max(sfr.store)),][3])
    ds.range <- as.numeric(cross.combs[which(sfr.store == max(sfr.store)),][4])
  }
  
  newlist <- list("us"=us.border, "ds"=ds.border, "range.us"=us.range, "range.ds"=ds.range, "flag"=TRUE)
  
  return(newlist)
  
}

SmoothProfile <- function(profile, bandwidth=5){
  # Smooth a profile given the specified badnwidth and a normal kernel
  #
  # Args:
  #   profile: smoothed profile 
  #   badnwidth: bandwidth for normal kernel smoothing (default=5)
  #
  # Returns:
  #   smoothed profile
  
    profile <- ksmooth(c(1:length(profile)), profile, kernel="normal", bandwidth=bandwidth)$y
    
    return(profile)
  
}

PruneProfile <- function(profile, desired.length){
  # Prune an average profile equally from both directions given the profile and the window size around the kmer which to retrieve
  # requires an even number as length to prune
  # lengthwill always be (desired.length + kmer.length)
  #
  # Args:
  #   profile: input profile
  #   desired.length: number of bps around the kmer to which the profile should be pruned
  #
  # Returns:
  #   pruned profile
  
  #first check if length is equal
  if (desired.length %% 2 != 0) {
    warning("The desired length is not even. Please select an even number to retrieve a profile symetric around the k-mer: returning NA")
    return(NA_real_)
  }
  
  remove.temp <- (250-desired.length)/2 #set positions to subtract from the full profile
  
  #prune
  profile <- profile[c((1+remove.temp)):(length(profile)-remove.temp)]
  
  #return
  return(profile)
  
}

PlotSingle <- function(profile, 
                       kl=6, 
                       plot.shoulders=FALSE, 
                       shoulders=FALSE, 
                       ylim=c(0,0.01), 
                       xlim=c(-125,125), 
                       color="black"){
  # Plot the average profile given an input profile
  # plot the shoulders if plot.shoulders is set to TRUE use shoulders list provided or determine
  #
  # Args:
  #   profile: input profile
  #   kl: k-mer length to plot dashed lines in the center
  #   plot.shoulders: flag if to plot the shoulder regions
  #   shoulders: estiamted shoulder border and ranges (list object as produced from SobelBorders)
  #   ylim: ylim to fix for plot (default c(0, 0.01))
  #   xlim: xlim to fix for plot (default c(-125, 125))
  #   clor: color to plot the profile
  #
  # Returns:
  #   Profile plot
  
  #check if shoulder ranges and border list is provided when shoulders should be plotted
  if((plot.shoulders == TRUE) && (shoulders == FALSE)){
    warning("There is no shoulder details list specified!\n Please run 
            SobelBorders or indicate to not plot the shoulders via plot.shoulders=FALSE!\n
            Will produce the plot without shoulders.")
    plot.shoulders <- FALSE
  }

  #cover if no border could be estimated but should be plotted
  if((plot.shoulders == TRUE) && (shoulders$flag == FALSE)){
    warning("No shoulders could be estimated in the provided shoulder list.\n
            Will plot the profile without shoulders.")
    plot.shoulders <- FALSE
  }

  #get window size of profile arround kmer
  window.size <- length(profile) - kl
  
  #make dataframe for ggplotting
  df <-data.frame(
    x=c(-(window.size/2):((window.size/2)+kl-1)),
    y=profile
  )
    
  #PLOTTING
  if(plot.shoulders){	#if plot.shoulders flag is TRUE
    
    #get borders and range from shoulders list
    us.border <- unlist(shoulders["us"])
    ds.border <- unlist(shoulders["ds"])
    us.range <- shoulders$range.us
    ds.range <- shoulders$range.ds
  
    #convert shoulder positions to relative positions surrounding kmer
    us.border <- us.border - (1+(window.size/2))
    ds.border <- ds.border - (1+(window.size/2))
      
    #MAKE PLOT
    p <- ggplot(df, aes(x=x, y=y)) + geom_line(size=1, color=color) + 
      geom_vline(xintercept = c((0),(kl-1)), linetype = "dashed", size=1) + 
      geom_vline(xintercept = c(
        (us.border+(us.range/2)), ds.border-(ds.range/2)
      ), colour="red", size=1) + 
      geom_vline(xintercept = c(
        (us.border-(us.range/2)), (ds.border+(ds.range/2))
      ), colour="red", linetype = "longdash", size=1) + 
      labs(x="Relative position [bp]", y="Relative cut probability") + 
      coord_cartesian(ylim=ylim, xlim=xlim) + 
      theme_bw() + science_theme + 
      theme(panel.grid = element_blank(), text = element_text(size=20))
    
    }else{		#border could not be estimated
    
    p <- ggplot(df, aes(x=x, y=y)) + geom_line(size=1, color=color) + 
        geom_vline(xintercept = c((0),(kl-1)), linetype = "dashed", size=1) + 
        labs(x="Relative position [bp]", y="Relative cut probability") + 
        coord_cartesian(ylim=ylim, xlim=xlim) + 
        theme_bw() + science_theme + 
        theme(panel.grid = element_blank(), text = element_text(size=20))
      
    
  }
  
  return(p)
    
}

PlotOverlap <- function(profile1, 
                        profile2, 
                        kmer1, 
                        kmer2,
                        count1 = "NA",
                        count2 = "NA",
                        ymode = "separate",
                        ylim = c(0,0.01), 
                        xlim = c(-125,125)){
  # Plot two average profiles overlapping given two input profiles
  # plot the shoulders if plot.shoulders is set to TRUE use shoulders list provided or determine
  #
  # Args:
  #   profile1: input profile1
  #   profile2: input profile2
  #   kmer1: input k-mer 1 (reference)
  #   kmer2: input k-mer 2 (variant)
  #   count1: count of k-mer1 occurence 
  #   count2: count of k-mer2 occurence 
  #   ymode: mode how to plot the overlapping profiles ("merged" or as "separate" profiles above each other)
  #   ylim: ylim to fix for plot (default c(0, 0.01))
  #   xlim: xlim to fix for plot (default c(-125, 125)) 
  #
  # Returns:
  #   Overlapping profile plot
  
  #check if kmer have equal length
  if(nchar(kmer1) != nchar(kmer2)){
    warning("kmers do bot have equal length")
    return("NA")
  }
  #check if ymode is set properly
  if((ymode != "merged") & (ymode != "separate")){
    warning("Select ymode (\"merged\" or \"separate\"")
    return("NA")
  }
  
  #get kmer length
  kl=nchar(kmer1)
  
  #get window size of profile arround kmer
  window.size <- length(profile1) - kl
  
  #make dataframe for ggplotting
  df <-data.frame(
    x=rep(c(-(window.size/2):((window.size/2)+kl-1)), 2),
    y=c(profile1, profile2),
    Source=c(
      rep(kmer1, length(profile1)), 
      rep(kmer2, length(profile2))
      )
  )
 
  #sort source for colors
  df$Source <- factor(df$Source, levels=c(kmer1, kmer2))
  
  #make dataframe for annotation
  anno.df <- data.frame(
    Source = c(kmer1, kmer2),
    x=rep((xlim[2]-(xlim[2]-xlim[1])/9), 2), 
    y=rep((ylim[2]-(ylim[2]-ylim[1])/12), 2),
    label=c(paste0(kmer1," #",count1), paste0(kmer2," #",count2))
  )
  
  #MAKE THE PLOT
  #mode one merged probfiles on same y-axis
  if(ymode == "merged"){
    
    #adjust annotation dataframe y values for nonoverlapping labels
    anno.df$y <- c((ylim[2]-(ylim[2]-ylim[1])/13), (ylim[2]-(ylim[2]-ylim[1])/8))
    
    p <- ggplot( df, aes(x=x, y=y, colour=Source)) + geom_line(size=1) + 
      geom_vline(xintercept = c((0),(kl-1)), linetype = "dashed", size=1) + 
      labs(x="Relative position [bp]", y="Relative cut probability") + 
      coord_cartesian(ylim=ylim, xlim=xlim) + 
      scale_colour_manual(values=rev(brewer.pal(3,"Set1")[c(1,2)])) +
      #add annotation
      geom_text(data=anno.df, aes(x=x, y=y, label=label)) +
      theme_bw() + science_theme + 
      theme(
        legend.position = "none", 
        # legend.title=element_blank(), 
        panel.grid = element_blank(), 
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16)
        )
    
  #mode two separate y-axis above each other  
  }else if(ymode == "separate"){
    
    #part one reference (top) without x-axis
    p <- ggplot(df, aes(x=x, y=y, colour=Source)) + 
      geom_line(size=1) + 
      scale_color_manual(values=rev(brewer.pal(3,"Set1")[c(1,2)])) +
      facet_grid(Source ~ .) +
      labs(x="Relative position [bp]", y="Relative cut probability") +
      geom_vline(xintercept = c((0),(kl-1)), linetype = "dashed", size=1) + 
      #add annotation
      geom_text(data=anno.df, aes(x=x, y=y, label=label)) +
      coord_cartesian(ylim=ylim, xlim=xlim) + 
      theme_bw() + science_theme + 
      theme(
        legend.position = "none",
        panel.grid.major.x = element_blank(), 
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        panel.border = element_blank(),
        panel.margin = unit("1.15", "lines"),
        strip.text.y = element_blank(),  #remove strips from facetting
        strip.background = element_blank()
      )

  }

  return(p)
  
}

QueryJaspar <- function(sequence, 
                        threshold=0.8, 
                        pwm.data){
  # Take a Sequence input and query it against the set of Jaspar2014 PWMs
  #
  # Args:
  #   sequence: input sequence
  #   threshold:percentage relative score thresh over which to report matches (default=0.8)
  #   pwm.data: a stored pwm. RData object as retrieved and save from JASPAR2014 R package
  #
  # Returns:
  #   Character string listing the PWM matches above a certain rel. threshold
  
  require(Biostrings)
  require(TFBSTools)
  
  # 1 pad sequence with N's
  sequence = paste0("NNNNN",sequence,"NNNNN")
  
  #make DNAString
  sequence = DNAString(sequence)
  
  #match
  suppressWarnings(out <- lapply(pwm.data, function(x) searchSeq(x, sequence, strand="*", min.score=threshold )) )
  out <- out[sapply(out, function(x) nrow(writeGFF3(x)) > 0 )]
  
  # combine output
  tmp <- ""
  
  if(length(out) >= 1){
    
    #make dataframes
    out.frame <- lapply(out, function(x) writeGFF3(x) ) #make data frame
    out.frame <- lapply(out.frame, function(x) x <- x[x[,6]==max(x[,6]),]) #pick highest matching strand
    factors <- sapply(out.frame, function(x) strsplit(as.character(unlist(x[9])),";")[[1]][1] )
    factors <- sapply(factors, function(x) gsub("TF=(\\w+)","\\1",x, perl=TRUE) )
    
    #get relative score
    rel.scores <- lapply(out, function(x) relScore(x))
    rel.scores <- lapply(rel.scores, max)
    
    #assemble data frame
    df <- data.frame(factor=factors, rel.scores=unlist(rel.scores) )
    df <- df[order(df$rel.scores, decreasing =TRUE),]
    
    #merge values
    for(i in c(1:nrow(df))){
      tmp <- paste0(tmp,df$factor[i],"=",round(df$rel.scores[[i]],digits=2),";")
    }
    
  }else{	#if no match above threshold report nothing
    tmp <- "no hits"
  }
  
  return(tmp)
  
}


### WRAPPER FUNCTIONS ####
GetFootprint <- function(kmer, 
                         tissue, 
                         data.dir, 
                         frag.type, 
                         smooth=TRUE, 
                         smooth.bandwidth=5){
  # Wrapper to retrieve merged, pruned, smoothed profile of kmer
  #
  # Args:
  #   kmer: input kmer 
  #   tissue: input tissue to query the profile from
  #   data.dir: repository storing processed kmer files per tissue
  #   frag.type: fragmentation type ("DNase" or "ATAC")
  #   smooth: flag if to smooth the profile (TRUE or FALSE)
  #   smooth.bandwidth: bandwidth to use for normal kernel smoothing 
  #   (default 5 which is fixed for most worklfows)
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
    profile.merge <- ksmooth(c(1:length(profile.merge)), profile.merge, kernel="normal", bandwidth=smooth.bandwidth)$y
  }
  
  return(list(profile=profile.merge, count=l.plus$count)) #return
  
}

GetFootprintStrand <- function(kmer, 
                               tissue="", 
                               data.dir="", 
                               frag.type, 
                               smooth, 
                               smooth.bandwidth=5, 
                               background.flag=FALSE){
  # Wrapper to retrieve stradnspecific (smoothed) profiles of kmer from tissue or background
  #
  # Args:
  #   kmer: input kmer 
  #   tissue: input tissue to query the profile from
  #   data.dir: repository storing processed kmer files per tissue
  #   frag.type: fragmentation type ("DNase" or "ATAC")
  #   smooth: flag if to smooth the profile (TRUE or FALSE)
  #   smooth.bandwidth: bandwidth to use for normal kernel smoothing 
  #     (default 5 which is fixed for most worklfows)
  #   background.flag: indicate if to retrieve background cut frequencies 
  #   (other data dir and file structure) [TRUE/FALSE] default FALSE
  # Returns:
  #   (smoothed) strand specific profiles and count of the kmer occurence
  
  #get kmer length
  kl=nchar(kmer)
  
  #select flat strand specific files as inputs (depending on background flag)
  if(!background.flag){
    
    infile.plus=file.path(data.dir, tissue,"counts", paste0("kmers_", kl, "_count_", tissue, "_pnorm_JH60_plus.txt"))
    infile.minus=file.path(data.dir, tissue,"counts", paste0("kmers_", kl, "_count_", tissue, "_pnorm_JH60_minus.txt"))
    
  }else if(background.flag){
    
    infile.plus=file.path(data.dir, tissue, "counts", paste0("kmer_", kl, "_", tissue, "_plus_merged"))
    infile.minus=file.path(data.dir, tissue, "counts", paste0("kmer_", kl, "_", tissue, "_minus_merged"))
    
  }else{
    warning("Specifiy if to retrieve data from the background naked DNA cut frequencies or from tissue\n
            background.flag=FALSE/TRUE !")
    return("NA")
  }
  
  #grep strand specific profiles & counts
  l.plus <- GrepProfile(kmer, infile.plus)
  l.minus <- GrepProfile(kmer, infile.minus)
  
  #smooth if specified so
  if(smooth){
    l.plus$profile <- ksmooth(c(1:length(l.plus$profile)), l.plus$profile, kernel="normal", bandwidth=smooth.bandwidth)$y
    l.minus$profile <- ksmooth(c(1:length(l.minus$profile)), l.minus$profile, kernel="normal", bandwidth=smooth.bandwidth)$y
  }
  
  return(list(profile.plus=l.plus$profile, profile.minus=l.minus$profile, 
              count.plus=l.plus$count, count.minus=l.minus$count)) #return
  
  
}

GetSFR <- function(kmer, 
                   tissue, 
                   data.dir="", 
                   vocab.flag=FALSE, 
                   vocab.file=paste0(data.dir,"/",tissue,"/vocabulary_",tissue,".txt"), 
                   frag.type=""){
  # Wrapper function to get the SFR ratio
  # If indicated and available, use the present vocabulary file to directly grep the SFR
  # Else get the average profile, estimate the borders and calculate the SFR
  # Note that for using the vocabulary file only nonambivlent DNA chars are allowed
  # For the alternative ambivalent chars are decoded
  #
  # Args:
  #   kmer: k-mer to query
  #   tissue: input tissue to query the profile from
  #   data.dir: repository storing processed kmer files per tissue
  #   vocab.flag: indicating if a preprocessed vocabulary file is present and to be used [TRUE/FALSE]
  #   vocab.file= path to preprocessed vocabulary file
  #   frag.type: fragmentation type ("DNase" or "ATAC")
  #
  # Returns:
  #   SFR
  
  #check if vocabulary flag selected and if so a file is present
  if(vocab.flag){
    
    #check if kmer is has ambivalent characters (Only )
    if(grepl("[^(A,C,G,T)]", kmer, perl=T)){
      warning("The Entered kmer hast ambivalent characters, for using the preprocessed vocabulary
              file please specify non-ambivalent chars or run the alternative vocab.flag=FALSE !")
      return(NA_real_)
    }
    #check if indicated vocab file exists
    if(!file.exists(vocab.file)){
      warning(paste0("File ", vocab.file, "does not exists! Please indicate the appropriate 
                     path to the vocabulary file or select vocab.flag=FALSE instead!"))
      return(NA_real_)
    }
    
    #grep kmer line from vocabulary table
    kmer.match <- grep(paste0("^",kmer,"\\s+"), readLines(vocab.file), value=TRUE, perl=TRUE)  #system grep of kmer vocab file
    
    split <- strsplit(kmer.match, "\t")[[1]]	#split on "\t" and unlist, as flat files are tab seperated
    sfr <- as.numeric(split[2])
    
    return(sfr)
    
  }else{
    
    #no vocabfile present --> run retrieving, border estimation and SFR calculation
    
    fp <- GetFootprint(kmer=kmer, tissue=tissue, data.dir=data.dir, frag.type=frag.type, 
                       smooth=TRUE) #get smoothed profile and count
    
    sh <- SobelBorders(profile=fp$profile, kl=nchar(kmer)) #estimate optimal shoulder midpoints and ranges
    
    #calculate SFR
    sfr <- CalcSFR(
      profile=fp$profile, 
      us.mid=sh$us, ds.mid=sh$ds,
      range.us=sh$range.us, range.ds=sh$range.ds
    )
    
    return(sfr)
    
  } 
  
    }

QueryLongSequence <- function(sequence, 
                              kl, 
                              tissue, 
                              data.dir="", 
                              vocab.flag=FALSE, 
                              vocab.file=paste0(data.dir,"/",tissue,"/vocabulary_",tissue,".txt"), 
                              frag.type="", 
                              plots=FALSE, 
                              smooth, 
                              plot.shoulders=TRUE, 
                              ylim=c(0,0.01), 
                              xlim=c(-125,125)){
  # Wrapper function to get split longer sequence into kmer of length kl and return kmer, SFR 
  # plots if specified
  #
  # Args:
  #   sequence: input logner sequence
  #   kl: length of k-mers to qsplit sequence into
  #   tissue: input tissue to query the profile from
  #   data.dir: repository storing processed kmer files per tissue
  #   vocab.flag: indicating if a preprocessed vocabulary file is present and to be used [TRUE/FALSE]
  #   vocab.file= path to preprocessed vocabulary file
  #   frag.type: fragmentation type ("DNase" or "ATAC")
  #   plots: TRUE/FALSE indicate if to make the profile plot for every k-mer
  #   smooth: flag if to smooth the profile (TRUE or FALSE)  
  #   plot.shoulders: flag if to plot the shoulder regions
  #   ylim: ylim to fix for plot (default c(0, 0.01))
  #   xlim: xlim to fix for plot (default c(-125, 125)) 
  #
  # Returns:
  #   data frame listing the splittet kmers with the respective SFR ($df)
  #   if specified list of plots one per splitted k-mer ($plots)
  
  #capture if sequence is to short
  if(nchar(sequence) < kl){
    warnings("Input sequence is shorter than the length of k-mers to split to! Returning NA!")
    return("NA")
  }
  
  #1 dissect sequence 
  dl <- DissectSequence(sequence, kl, list=TRUE)
  
  #2 get SFR per k-mer 
  dl.sfr <- lapply(dl, function(x){
    x <- GetSFR(kmer=x, tissue=tissue, data.dir=data.dir, vocab.flag = vocab.flag, vocab.file = vocab.file, frag.type = frag.type)
  })
  
  #3 compose data frame for output
  dl.df <- data.frame(kmer=unlist(dl), sfr=unlist(dl.sfr))
  
  #4 get plots if specified
  if(plots){
    dl.plots <- lapply(dl, function(x){
      x <- PlotSingleKmer(kmer=x, tissue=tissue, data.dir=data.dir, frag.type=frag.type, 
                          smooth=smooth, plot.shoulders=plot.shoulders, ylim=ylim, xlim=xlim)
    })
  }
  
  #return according to specified desired output
  if(plots){
    newlist <- list(df=dl.df, plots=dl.plots)
    return(newlist)
  }else if(!plots){
    return(dl.df)
  }else{
    warnings("Specifiy if to make plots explicitly if they are desired! plots=TRUE/FALSE  Will return output without plots!")
    return(dl.df)
  }
}

CompareSequences <- function(sequence1, 
                             sequence2, 
                             kl, 
                             damage.mode="exhaustive", 
                             tissue, 
                             data.dir="", 
                             vocab.flag=FALSE, 
                             vocab.file=paste0(data.dir,"/",tissue,"/vocabulary_",tissue,".txt"), 
                             frag.type="", 
                             plots="highest", 
                             smooth=TRUE, 
                             ylim=c(0,0.01), 
                             xlim=c(-125,125)){
  # Wrapper function to analyse multiple Ref-Var-Sequence pairs
  # split each sequence into k-mers of lenfth kl, get their SFRs form vocab file 
  # or calculate new and calculate the damage associated with each kmer pair and 
  # from that the local or exhaustive summed up dmage of entire sequence pair
  #
  # Args:
  #   ref.var.df: three column data frame listing id referecne and variance sequence 
  #   in the following scheme (id reference variance)
  #   kl: length of k-mers to split sequences into
  #   damage.mode: mode for calculating the total damage 
  #   (local: get highest pair, exhaustive: sum up over entire sequence)
  #   tissue: input tissue to query the profile from
  #   data.dir: repository storing processed kmer files per tissue
  #   vocab.flag: indicating if a preprocessed vocabulary file is present and to be used [TRUE/FALSE]
  #   vocab.file: path to preprocessed vocabulary file
  #   frag.type: fragmentation type ("DNase" or "ATAC")
  #   plots: indicate if to retrie overlap plots as well 
  #     with values FALSE for no plots to retrieve, "highest" for only retrieving one plot of 
  #     the highest scoring kmer pair, "all" for retrieving a list of all pairwise overlay plots
  #
  # Returns:
  #   Dataframe listing Ref and Var sequence with highest scoring kmer pair and 
  #   the according SFRs and calculated (exhaustive or local) damage 
  
  #check if indicated vocab file is readable
  if((vocab.flag == TRUE) & (!file.exists(vocab.file))){
    warning(paste0("File ", vocab.file, "does not exists! Please indicate the appropriate 
                     path to the vocabulary file or select vocab.flag=FALSE instead!"))
    return(NA_real_)
  }
  
  #capture if sequences are to short
  if((nchar(sequence1) < kl) | (nchar(sequence2) < kl)){
    warnings("At least one input sequence is shorter than the length of k-mers to split to! Returning NA!")
    return("NA")
  }
  #check if sequences have equal length
  if(nchar(sequence1) != nchar(sequence2)){
    warning("Input sequences have unequal length! Currently only equal length is suported! Returning NA!")
    return("NA")
  }
  
  
  #1 split sequences into two lists
  dl1 <- DissectSequence(sequence1, kl, list=TRUE)
  dl2 <- DissectSequence(sequence2, kl, list=TRUE)
  
  #2 get SFRs for each k-mer 
  dl1.sfr <- lapply(dl1, function(x){
    x <- GetSFR(kmer=x, tissue=tissue, data.dir=data.dir, vocab.flag = vocab.flag, vocab.file = vocab.file, frag.type = frag.type)
  })
  dl2.sfr <- lapply(dl2, function(x){
    x <- GetSFR(kmer=x, tissue=tissue, data.dir=data.dir, vocab.flag = vocab.flag, vocab.file = vocab.file, frag.type = frag.type)
  })
  
  #3 compose data frame
  dl.df <- data.frame(kmer.ref=unlist(dl1), kmer.var=unlist(dl2), sfr.ref=unlist(dl1.sfr), sfr.var=unlist(dl2.sfr))
  
  #4 callc single kmer overlap damag
  dl.df$damage <- dl.df$sfr.ref - dl.df$sfr.var
  
  #5 pick highest scoring kmers
  highest.id <- which(abs(dl.df$damage) == max(abs(dl.df$damage)))[1]
  
  #6 calculate total damage according to selected mode
  if(damage.mode == "exhaustive"){
    total.damage <- sum(dl.df$damage)
  }else if(damage.mode == "local"){
    total.damage <- dl.df[highest.id, ]$damage
  }else{
    warning("Specify a mode to calculate the total damage for comparing the sequences! damage.mode=exhaustive/local")
    return("NA")
  }
  
  #7 make summary line like for batch query
  summary.line <- data.frame(
    sequence.ref=sequence1,
    sequence.var=sequence2,
    kmer.ref=as.character(dl.df[highest.id, ]$kmer.ref), 
    kmer.var=as.character(dl.df[highest.id, ]$kmer.var),
    SFR.ref=dl.df[highest.id, ]$sfr.ref,
    SFR.var=dl.df[highest.id, ]$sfr.var,
    total.damage=total.damage
    )
  
  
  #8 make plots if desired 
  if(plots == "all"){
    
    plot.list <- apply(dl.df[, c("kmer.ref", "kmer.var")], 1, function(x){
      
      pp <- PlotOverlapKmers(
        kmer1=x[1], 
        kmer2=x[2],
        tissue1=tissue, 
        tissue2=tissue, 
        data.dir=data.dir, 
        frag.type=frag.type, 
        smooth=smooth,
        ylim=ylim, 
        xlim=xlim
        )
        
      return(pp)
      
    })

  }else if(plots == "highest"){
    
    plot.list <- PlotOverlapKmers(
      kmer1=as.character(dl.df[highest.id, ]$kmer.ref), 
      kmer2=as.character(dl.df[highest.id, ]$kmer.var),
      tissue1=tissue, 
      tissue2=tissue, 
      data.dir=data.dir, 
      frag.type=frag.type, 
      smooth=smooth, 
      ylim=ylim, 
      xlim=xlim
      )
    
  }else if(!plots){
    plot.list <- "No Plots specified"
  }else{
    warnings("Please select if and what plots to produce: FALSE \"highest\" \"all\"!\n Will return without plots!")
    plot.list <- "No Plots specified"
  }
        
  #8 assemble output
  newlist <- list(df=dl.df, summary=summary.line, plots=plot.list, damage.mode=damage.mode)
  return(newlist)
  
}

RefVarBatch <- function(ref.var.df, 
                        kl, 
                        damage.mode="exhaustive", 
                        tissue, 
                        data.dir, 
                        vocab.flag, 
                        vocab.file=paste0(data.dir,"/",tissue,"/vocabulary_",tissue,".txt"), 
                        frag.type=""){
  # Wrapper function to analyse multiple Ref-Var-Sequence pairs
  # split each sequence into k-mers of lenfth kl, get their SFRs form vocab file 
  # or calculate new and calculate the damage associated with each kmer pair and 
  # from that the local or exhaustive summed up dmage of entire sequence pair
  #
  # Args:
  #   ref.var.df: three column data frame listing id referecne and variance sequence 
  #   in the following scheme (id reference variance)
  #   kl: length of k-mers to split sequences into
  #   damage.mode: mode for calculating the toal damage 
  #   (local: get highest pair, exhaustive: sum up over entire sequence)
  #   tissue: input tissue to query the profile from
  #   data.dir: repository storing processed kmer files per tissue
  #   vocab.flag: indicating if a preprocessed vocabulary file is present and to be used [TRUE/FALSE]
  #   vocab.file= path to preprocessed vocabulary file
  #   frag.type: fragmentation type ("DNase" or "ATAC")
  #
  # Returns:
  #   Dataframe listing Ref and Var sequence with highest scoring kmer pair and 
  #   the according SFRs and calculated (exhaustive or local) damage 
  
  #check format of input.df
  if(ncol(ref.var.df) != 3){
    warnign("The input ref.var.df data.frame requires a three column data frame input format (id ref var)!")
    return("NA")
  }
  
  #1 perform compare sequence for every row
  temp.list <- apply(ref.var.df, 1, function(x){

    comp <- CompareSequences(sequence1=as.character(x[2]), sequence2=as.character(x[3]), kl=kl, damage.mode=damage.mode, 
                             tissue=tissue, data.dir=data.dir, vocab.flag=vocab.flag, 
                             vocab.file=vocab.file, frag.type=frag.type, plots=FALSE)$summary
    return(comp)
    
  })
  
  #2 rbind to output summary data frame
  out.df <- do.call(rbind, temp.list)
  #3 add id columns
  out.df <- cbind(ref.var.df[, 1], out.df)
  names(out.df)[1] <- "id"
   
  return(out.df)
  
}

QueryJasparBatch <- function(df, 
                             damage.threshold=0, 
                             match.threshold=0.8, 
                             pwm.data){
  # Take a data frame from RefVar Query Batch as input and query it against the 
  # set of Jaspar2014 PWMs using a selected match.threshold 
  # Select an absolute sasq damage above which to query jaspar
  #
  # Args:
  #   df: input dataframe from RefVarBatch query 
  #   damage.threshold: absolute damage threshold above which a RefVar pair is queried
  #   match.threshold: percentage relative score thresh over which to report matches (default=0.8)
  #   pwm.data: a stored pwm. RData object as retrieved and save from JASPAR2014 R package
  #
  # Returns:
  #   Dataframe with additional column for jaspar query results
  
  #check input dataframe
  if(ncol(df) != 8){
    warning("Input dataframe df does not have 8 columns! Please make sure RefVarBatch has run properly:\n
            Format: id sequence.ref  sequence.var kmer.ref kmer.var  SFR.ref  SFR.var total.damage")
    return("NA")
  }
  
  
  
  # go throug data frame and query jaspar
  jsp <- apply(df, 1, function(x){

    
    if(abs(as.numeric(x[8])) >= damage.threshold){ #only query if absolute total damage is above threshold
    
      if(x[8] >= 0){ # if total.damage is positive --> quer reference sequence
     
         j <- QueryJaspar(sequence=x[2], threshold=match.threshold, pwm.data=pwm.data)
    
      }else if(x[8] <= 0){# if total.damage is negative query variant sequence
        
          j <- QueryJaspar(sequence=x[3], threshold=match.threshold, pwm.data=pwm.data)
          
      }
    
    }else{#dont query
      
      j <- "."
      
    }
    
    return(j) #return
    
  })
  
  #add new column
  df$jaspar <- jsp
  
  #return data frame
  return(df)
  
}

#PLOT FUNCTION WRAPPER
PlotSingleKmer <- function(kmer, 
                           tissue, 
                           data.dir, 
                           frag.type, 
                           smooth=TRUE, 
                           smooth.bandwidth=5, 
                           plot.shoulders=FALSE, 
                           ylim=c(0,0.01), 
                           xlim=c(-125,125), 
                           color="black"){
  #Wrapper to produce a plot from kmer and tissue input only
  #
  # Args:
  #   kmer: input kmer 
  #   tissue: input tissue to query the profile from
  #   data.dir: repository storing processed kmer files per tissue
  #   frag.type: fragmentation type ("DNase" or "ATAC")
  #   smooth: flag if to smooth the profile (TRUE or FALSE)  
  #   plot.shoulders: flag if to plot the shoulder regions
  #   ylim: ylim to fix for plot (default c(0, 0.01))
  #   xlim: xlim to fix for plot (default c(-125, 125)) 
  #   color: color for profile to plot
  #
  # Returns:
  #   Plot of profile
  
  #1 Get smoothed or unsmoothed profile and shoulders
  fp <- GetFootprint(kmer=kmer, tissue=tissue, data.dir=data.dir, frag.type=frag.type, smooth=smooth, smooth.bandwidth=smooth.bandwidth)
  
  if(plot.shoulders){
    sh <- SobelBorders(fp$profile, kl=nchar(kmer))
  }  
  
  #2 make plot
  p <- PlotSingle(profile=fp$profile, kl=nchar(kmer), plot.shoulders=plot.shoulders, shoulders=sh, ylim=ylim, xlim=xlim, color=color)
  
}

PlotOverlapKmers <- function(kmer1, 
                             kmer2, 
                             tissue1, 
                             tissue2, 
                             data.dir, 
                             frag.type, 
                             smooth=TRUE, 
                             ylim=c(0,0.01), 
                             xlim=c(-125,125)){
  #Wrapper to produce an overly plot from two kmers and tissues input only
  #
  # Args:
  #   kmer1: input k-mer1
  #   kmer2: input k-mer2 
  #   tissue1: input tissue1 to query the profile from
  #   tissue2: input tissue2 to query the profile from
  #   data.dir: repository storing processed kmer files per tissue
  #   frag.type: fragmentation type ("DNase" or "ATAC")
  #   smooth: flag if to smooth the profile (TRUE or FALSE)  
  #   ylim: ylim to fix for plot (default c(0, 0.01))
  #   xlim: xlim to fix for plot (default c(-125, 125)) 
  #
  # Returns:
  #   Overlay plot of profiles
  
  #check if smooth flag set properly
  if((smooth != TRUE) & (smooth != FALSE)){
    warning("Please specify if to smooth the overlay plots: smooth=TRUE/FALSE!")
    return("NA")
  }
  
  #1 Get smoothed profiles  
  fp1 <- GetFootprint(kmer=kmer1, tissue=tissue1, data.dir=data.dir, frag.type=frag.type, smooth=smooth)
  fp2 <- GetFootprint(kmer=kmer2, tissue=tissue2, data.dir=data.dir, frag.type=frag.type, smooth=smooth)

  
  #2 Make plot
  p <- PlotOverlap(profile1 = fp1$profile, profile2 = fp2$profile, 
                   kmer1 = kmer1, kmer2 = kmer2, 
                   count1 = fp1$count, count2 = fp2$count, 
                   ylim=ylim, xlim=xlim)  
  
  return(p)
  
}

PlotSingleStrands <- function(kmer, 
                              tissue, 
                              data.dir, 
                              frag.type, 
                              smooth=TRUE, 
                              smooth.bandwidth=5, 
                              background.flag=FALSE, 
                              ylim=c(0,0.01), 
                              xlim=c(-125,125)){
  #Wrapper to produce a strand specific plots from kmer and tissue input only
  #
  # Args:
  #   kmer: input kmer 
  #   tissue: input tissue to query the profile from
  #   data.dir: repository storing processed kmer files per tissue
  #   frag.type: fragmentation type ("DNase" or "ATAC")
  #   smooth: flag if to smooth the profile (TRUE or FALSE)  
  #   plot.shoulders: flag if to plot the shoulder regions
  #   ylim: ylim to fix for plot (default c(0, 0.01))
  #   xlim: xlim to fix for plot (default c(-125, 125)) 
  #
  # Returns:
  #   Plot of profile for plus ($plot.plus) and minus strand ($plot.minus)
  
  #1 Get smoothed profile and shoulders
  if(smooth){
    
    fp <- GetFootprintStrand(kmer=kmer, tissue=tissue, data.dir=data.dir, frag.type=frag.type, smooth=smooth, smooth.bandwidth=smooth.bandwidth, background.flag=background.flag)
    
  }else if(!smooth){
    
    fp <- GetFootprintStrand(kmer=kmer, tissue=tissue, data.dir=data.dir, frag.type=frag.type, smooth=smooth, smooth.bandwidth=smooth.bandwidth, background.flag=background.flag)
    
  }
  
  #2 make plot
  pplus <- PlotSingle(profile=fp$profile.plus, kl=nchar(kmer), plot.shoulders=FALSE, shoulders=sh, ylim=ylim, xlim=xlim)
  pminus <- PlotSingle(profile=fp$profile.minus, kl=nchar(kmer), plot.shoulders=FALSE, shoulders=sh, ylim=ylim, xlim=xlim)
  
  #return
  newlist <- list(plot.plus=pplus, plot.minus=pminus)
  return(newlist)
  
}

GetPossibleMutations <- function(sequence, 
                                 kl=7, 
                                 chr=".", 
                                 position=1){
  # Take a Sequence input and split into kmers of length kl*2-1 with reference and variance.
  # Take the parsed position as index for the first base to mutate (kl'th sequence to fill window)
  # (e.g. for kl 7 always take the 13 surrounding bases)
  #
  # Args:
  #   sequence: input sequence
  #   kl: k-mer lengh to split the sequence into k-mers
  #   chr: chromosome
  #   position: Start base position where to predict the 
  #
  # Returns:
  #   Six columns dataframe c(chr, position, ref.base, var.base, ref.sequence, var.sequence).

  #check if sequence is long enough
  if(nchar(sequence) < (kl*2-1)){
    warning("Sequence is to short! Has to be a minimum of: kl*2-1 !")
    return(NA_real_)
  }

  window.size <- kl*2-1
  
  #make moving window sequence split
  seq.list <- c()
  for(i in c(1:(nchar(sequence)-window.size+1))){
    
    seq.list <- c(seq.list, substr(sequence, i, i+window.size-1))
  
  }
  
  #get ref base per sequqnce window
  ref.list <- substr(seq.list, kl, kl)
  
  bases <- c("A", "C", "G", "T")
  
  #make dataframe
  df <- data.frame(
    
    chr = rep(chr, length(seq.list)*3),
    
    pos = as.numeric(unlist( lapply( c(position:(length(seq.list)+position-1)), function(x) rep(x, 3) ))),
    
    ref.base = unlist( lapply( ref.list, function(x) rep(x, 3) )) 
    
  )
  
  #add all possible substitutions
  temp <- c()
  for(j in c(position:(length(seq.list)+position-1))){
    
   temp <- c(temp, bases[!(bases %in% c(as.character(df[df$pos == j, ]$ref.base[1])))])
   
  }
  df$var.base <- temp
  
  #add ref sequence windows
  df$ref.seq <- unlist( lapply( seq.list, function(x) rep(x, 3) ))
  
  #add mutated sequence windows
  df$var.seq <- df$ref.seq 
  substr(df$var.seq, kl, kl) <- df$var.base
  
  #return data frame
  return(df)
  
}

InSilicoMutation <- function(  sequence,
                                      kl=7,
                                      chr=".",
                                      position=1,
                                      report="all",
                                      damage.mode="exhaustive",
                                      tissue=tissue,
                                      data.dir=data.dir,
                                      vocab.flag=TRUE,
                                      vocab.file=paste0(data.dir,"/",tissue,"/vocabulary_",tissue,".txt"),
                                      frag.type=frag.type){
# Wrapper for max/abs damage insilico mutation
# Take a sequence, split into dataframe of kmerlength matching window and 
# compare reference an possible mutation sequences
# report according to report mode ("all", "max", "maxabs")
#
# Args:
#   sequence: input sequence
#   kl: k-mer lengh to split the sequence into k-mers
#   chr: chromosome
#   position: Start base position where to predict the base substitution
#   report: which mutations to report; 
#     "all" = report all 3 possible substitutions per position
#     "max" = only report substitution with highest positive damage
#     "maxabs" = only report substitution with highest absolute damage
#   damage.mode: mode for calculating the total damage 
#   tissue: input tissue to query the profile from
#   data.dir: repository storing processed kmer files per tissue
#   vocab.flag: indicating if a preprocessed vocabulary file is present and to be used [TRUE/FALSE]
#   vocab.file: path to preprocessed vocabulary file
#   frag.type: fragmentation type ("DNase" or "ATAC")
#
# Returns:
#   Seven columns dataframe c(chr, position, ref.base, var.base, ref.sequence, var.sequence, damage).  

require(pbapply)  
  
# make split sequence dataframe for query
df <- GetPossibleMutations(sequence=sequence, kl=7, chr=chr, position=position)

print(paste0("Processing ",nrow(df)," sequence windows:"))

#inslico mutation for set mode
df$damage <- pbapply(df, 1, function(x) CompareSequences(
  sequence1=x[5], 
  sequence2=x[6], 
  kl=kl, 
  damage.mode=damage.mode,
  tissue=tissue, 
  data.dir=data.dir,
  vocab.flag=vocab.flag,
  vocab.file=vocab.file,
  frag.type=frag.type,
  plots=FALSE
  )$summary$total.damage
  )

#filter dataframe according to report mode
if(report == "all"){

   return(df)

}else if(report == "max"){
  
  dftemp <- df[0,]
  
  #select max of 3 matching rows
  for(i in df$pos[seq(from=3, to=nrow(df), by=3)]){
    
    temp <- df[df$pos %in% i, ]
    
    temp <- temp[which(temp$damage == max(temp$damage)), ]
    
    if(dim(temp)[1] > 1){  #sample randomly if equal damages
      temp <- temp[sample(c(1:dim(temp)[1]),1), ]
    }
    
    dftemp <- rbind(dftemp, temp)
    
  }
  
  return(dftemp)
  
}else if(report == "maxabs"){
  
  dftemp <- df[0,]
  
  #select max of absolute of 3 matching rows
  for(i in df$pos[seq(from=3, to=nrow(df), by=3)]){
    
    temp <- df[df$pos %in% i, ]
    
    temp <- temp[which(temp$damage == max(abs(temp$damage))), ]
    
    if(dim(temp)[1] > 1){  #sample randomly if equal damages
      temp <- temp[sample(c(1:dim(temp)[1]),1), ]
    }
    
    dftemp <- rbind(dftemp, temp)
    
  }
  
  return(dftemp)
  
}else{
  
  warning("Select a mode for reporting! Will report default (all possible substitutions)!")
  return(df)

}
  
}

#wrapper for rainbow plot
RainbowPlot <- function(df, ylim=c(-2,2)){
  # Wrapper for max/abs damage insilico mutation
  # Take a sequence, split into dataframe of kmerlength matching window and 
  # compare reference an possible mutation sequences
  # report according to report mode ("all", "max", "maxabs")
  #
  # Args:
  #   df: Dataframe from InSilicoMutation (report="all")
  #   ylim: y-limits for plot
  #
  # Returns:
  #   Seven columns dataframe c(chr, position, ref.base, var.base, ref.sequence, var.sequence, damage).  
  
  
  #check data frame
  
  #make plot
  p <- ggplot(df, aes(x=pos, y=damage, col=var.base)) + 
  geom_hline(yintercept=0) +
  geom_segment(aes(x=pos, xend=pos, y=0, yend=damage), col="Grey") + 
  geom_point(size=2) + 
  coord_cartesian(xlim=c(df$pos[1], tail(df$pos, 1)), ylim=ylim) + 
  labs(x=df$chr[1], y="Sum of Damage") +
  scale_color_manual(values=brewer.pal(6, "Set1")[c(1,2,4,3)], name="Variant") + 
  theme_bw() + science_theme +
  theme(
    panel.grid.major.x=element_blank(), 
    panel.border = element_blank()
  )

  return(p)
  
}
