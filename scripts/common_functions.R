### CONTAINS ALL FUNCTION NECESSARY TO PRODUCE THE KMER PROFILE DATABASE ###
### hopefully ###

################################################################
###create list of kmers segregating FaStA code MINIMUM  5mer ###
################################################################
decode.kmer <- function(kmer){

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


###############################################################################
### Calculate the simplified adjusted FSR score using the estimated borders ###
###############################################################################
adj.fsr <- function (profile, kl, usb, dsb, usr, dsr){

	footprint <- profile[c((usb + usr/2):(dsb-dsr/2))]  #get cuts in footprint based on estimated borders
	fl <- length(footprint)
	footprint.cuts <- sum(footprint)	#calc average cuts

	shoulder <- profile[ c( c( (usb - usr/2):(usb + usr/2) ), c( (dsb - dsr/2):(dsb + dsr/2) ) ) ]	#get cuts in shoulders based on est borders
	sl <- length(shoulder)
	shoulder.cuts <- sum(shoulder)	#calc average cuts

	#calc adj fsr
	fsr <- (shoulder.cuts/sl)/(footprint.cuts/fl)
	return(fsr)

}

################################################################################################
### Function to scale Seq Specific Variations, Gaussian Smooth and Substract them afterwards ###
################################################################################################
gaussian.median.diff <- function(profile, profile.naked, kl, skip.flag, extension){

#handle variation outside kmer window to flag as such and skip normalization
	profile.naked.temp <- ksmooth(c(1:length(profile.naked)), profile.naked, kernel="normal", bandwidth=5)$y
	profile.naked.temp <- profile.naked.temp/sum(profile.naked.temp)
#read skip flag to decide if and how to normalize
	if(skip.flag == 2 ){	
		norm.flag <- FALSE
		profile <- ksmooth(c(1:length(profile)), profile, kernel="normal", bandwidth=5)$y
		profile <- profile/sum(profile)
		newlist <- list("profile"=profile, "scale.factor"=0, "norm.flag"=norm.flag)
		return(newlist)
	}else if(skip.flag == 1){
		norm.flag <- TRUE
	#get center regions
		trimm.profile <- profile[c((126-extension):(126+kl-1+extension))] #kmer region only
		trimm.profile.naked <- profile.naked[c((126-extension):(126+kl-1+extension))]
	#derive scale factor
		scale.factor <- sd(trimm.profile)/sd(trimm.profile.naked)
	#apply scale factor to naked
		profile.naked <- profile.naked * scale.factor
	#difference from median
		difftomedian <- profile.naked - median(profile.naked)
		profile <- profile - difftomedian
	#smooth and normalize back to relative cut frequency
		profile <- ksmooth(c(1:length(profile)), profile, kernel="normal", bandwidth=5)$y
		profile <- profile/sum(profile)

		newlist <- list("profile"=profile, "scale.factor"=scale.factor, "norm.flag"=norm.flag)
		return(newlist)
	}else{
		norm.flag <- TRUE
		trimm.profile <- profile[c(125:(126+kl-1+1))] #kmer region only
		trimm.profile.naked <- profile.naked[c(125:(126+kl-1+1))]
	#derive scale factor
		scale.factor <- sd(trimm.profile)/sd(trimm.profile.naked)
	#apply scale factor to naked
		profile.naked <- profile.naked * scale.factor
	#difference from median
		difftomedian <- profile.naked - median(profile.naked)
		profile <- profile - difftomedian
	#smooth and normalize back to relative cut frequency
		profile <- ksmooth(c(1:length(profile)), profile, kernel="normal", bandwidth=5)$y
		profile <- profile/sum(profile)

		newlist <- list("profile"=profile, "scale.factor"=scale.factor, "norm.flag"=norm.flag)
		return(newlist)
	}
}

########################################################################################################################
### FUNCTION wrapping sobeld 1st derivative approximation, zero crossing detection and border estimation baed on that###
########################################################################################################################
sobel.borders <- function(profile, kl){

	profile.sobel <- sobeln(profile) #approx 1st derivative via 1D discrete sobel operator

	zero.cross = which(diff(sign(profile.sobel)) == -2)	#get 0 crossings with correct directionality > maxima

	us.cross <- subset(zero.cross, zero.cross <= (126))	#seperate crossings into upstream and downstream
	ds.cross <- subset(zero.cross, zero.cross >= (126 + kl - 1))

	us.cross <- subset(us.cross, us.cross >= 76)	#only consider crossings in a reasonable window to speed up
	ds.cross <- subset(ds.cross, ds.cross <= (175+kl))

	peak.range <- c(6)

	cross.combs <- as.matrix(expand.grid(us.cross, ds.cross, peak.range, peak.range))	#get a matrix with all combinations of the relevant crossings
	#cover cases with no relevant maxima 
	if(dim(cross.combs)[1] == 0) {
		newlist <- list("us"=0, "ds"=0, "range.us"=0, "range.ds"=0, "flag"=FALSE)
		return(newlist)
	}
#first round get optimal maximas
	for(j in c(1:dim(cross.combs)[1])){	#calculcate fos for the range of shoulder widths and extract the amx one from each
		fos.store <- apply(cross.combs, 1, function(x){ u=x[1]; d=x[2]; ru=x[3]; rd=x[4]; fos.for.max(u, d, profile, ru, rd)})
	}

	us.border <- as.numeric(cross.combs[which(fos.store == max(fos.store)),][1])
	ds.border <- as.numeric(cross.combs[which(fos.store == max(fos.store)),][2])

	while( ((us.border + (2)) >= (ds.border - (2))) || ( (us.border ) > 126 ) || ( (ds.border ) < (126 + kl - 1) ) ){
		if(length(fos.store) <= 1) {	#no suitable borders found
			newlist <- list("us"=0, "ds"=0, "range.us"=0, "range.ds"=0, "flag"=FALSE)
			return(newlist)
		}
		fos.store <- fos.store[-which(fos.store == max(fos.store))]
		us.border <- as.numeric(cross.combs[which(fos.store == max(fos.store)),][1])
		ds.border <- as.numeric(cross.combs[which(fos.store == max(fos.store)),][2])
	}
#second round get optimal range	
	peak.range <- c(4,6,8,10)
	cross.combs <- as.matrix(expand.grid(us.border, ds.border, peak.range, peak.range))
	for(j in c(1:dim(cross.combs)[1])){	#calculcate fos for the range of shoulder widths and extract the amx one from each
		fos.store <- apply(cross.combs, 1, function(x){ u=x[1]; d=x[2]; ru=x[3]; rd=x[4]; fos.for.max(u, d, profile, ru, rd)})
	}
	fos.maxima <- max(fos.store)
	us.range <- as.numeric(cross.combs[which(fos.store == max(fos.store)),][3])
	ds.range <- as.numeric(cross.combs[which(fos.store == max(fos.store)),][4])
	
	#cover cases with overlapping shoulders
	while( ((us.border + (us.range/2)) >= (ds.border - (ds.range/2))) || ( (us.border + (us.range/2)) > (126 + 2 ) ) || ( (ds.border - (ds.range/2)) < (126 + kl - 1 - 2) ) ){
		if(length(fos.store) <= 1) {	#no suitable borders found
			newlist <- list("us"=0, "ds"=0, "range.us"=0, "range.ds"=0, "flag"=FALSE)
			return(newlist)
		}
		fos.store <- fos.store[-which(fos.store == max(fos.store))]

		us.border <- as.numeric(cross.combs[which(fos.store == max(fos.store)),][1])
		ds.border <- as.numeric(cross.combs[which(fos.store == max(fos.store)),][2])
		fos.maxima <- max(fos.store)
		us.range <- as.numeric(cross.combs[which(fos.store == max(fos.store)),][3])
		ds.range <- as.numeric(cross.combs[which(fos.store == max(fos.store)),][4])
	}
	
	newlist <- list("us"=us.border, "ds"=ds.border, "range.us"=us.range, "range.ds"=ds.range, "flag"=TRUE)
	return(newlist)

}

###################################################################################
### CALCULATE THE FOS SCORE USING A >10< b window around each crossing as "peak"###
###################################################################################
fos.for.max <- function(u, d, profile, range.us, range.ds){

	#window <- 6
	ext.us <- range.us/2
	ext.ds <- range.ds/2

	shoulder <- profile[ c( c((u-ext.us):(u+ext.us)), c((d-ext.ds):(d+ext.ds)) ) ]  	#extend the region around the crossing to get an approximation 
	shoulder.cuts <- sum(shoulder)	#of the peak itself, retireve profile values in that window
	footprint <- profile[c((u+ext.us+1):(d-ext.ds-1))]	#region in between = footprint
	footprint.cuts <- sum(footprint)

	fos <- (shoulder.cuts/length(shoulder))/((footprint.cuts)/length(footprint))
	return(fos)
}

######################################################
### SOBEL FUNCTION FOR 1st DERIVATIVE APPROXIMATION###
######################################################
sobeln<-function(x){

	b=rep(0,length(x))
	for(i in c(2:(length(x)-1))){
		b[i] = ( x[i-1] * -1 ) + ( x[i] * 0 ) + ( x[i+1] * 1 )
	}
	b[1] = b[2]
	b[length(x)] = b[length(x)-1]
	return(b)
}

#################################################################
### PLOT FUNCTION BACKGROUND INCLUDING COUNT OF OCCURENCE     ###
#################################################################
plot.profile.count.noborders <- function(kmer,combined.list){

	kl=nchar(kmer)
	kmer.list <- decode.kmer(kmer)
	count.plus <- 0
	count.minus <- 0
	profile.plus <- rep(0,(250+kl))
	profile.minus <- rep(0,(250+kl))

	for(k in kmer.list){
		count.plus <- count.plus + combined.list[[k]]$count.plus
		count.minus <- count.minus + combined.list[[k]]$count.minus
		profile.plus <- profile.plus + combined.list[[k]]$profile.plus 
		profile.minus <- profile.minus + combined.list[[k]]$profile.minus	
	}
	
	#smooth
	profile.plus <- ksmooth(c(1:length(profile.plus)), profile.plus, kernel="normal", bandwidth=5)$y
	profile.minus <- ksmooth(c(1:length(profile.minus)), profile.minus, kernel="normal", bandwidth=5)$y
	#calculate cut probability in that region
	profile.plot.plus <- profile.plus/sum(profile.plus)
	profile.plot.minus <- profile.minus/sum(profile.minus)

#estimatefootprint and plot in output dir if lot flag set TRUE
	list.plus <- sobel.borders(profile=profile.plot.plus,kl=kl)
	list.minus <- sobel.borders(profile=profile.plot.minus,kl=kl)
#cover if no border could be estimated
	border.flag.plus <- list.plus$flag
	border.flag.minus <- list.minus$flag
#get best shoulder ranges
	us.range.plus <- list.plus$range.us
	ds.range.plus <- list.plus$range.ds
	us.range.minus <- list.minus$range.us
	ds.range.minus <- list.minus$range.ds

##PLOTTING
	pp <- qplot(x=c(1:length(profile.plot.plus)),y=profile.plot.plus,geom='line', ylab="cut probability",main=paste0("Cut prob. profile ",kmer, ' #', count.plus)) + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))
  
	pm <- qplot(x=c(1:length(profile.plot.minus)),y=profile.plot.minus,geom='line', ylab="cut probability",main=paste0("Cut prob. profile ",kmer, ' #', count.minus )) + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))

	newlist <- list( "pp"=pp, "pm"=pm)
	return(newlist)

}


#tiny function to create the reverse complement (inverted) of a kmer
reverse.kmer <- function(kmer){
	mersplit <- unlist(strsplit(kmer,''))
	#create reverse complement
	rc <- paste0(rev(unlist(lapply(mersplit,function(x){
		if(x == 'A'){ x='T'; }
		else if(x == 'G'){ x='C'; }
		else if(x == 'C'){ x='G'; }
		else if(x == 'T'){ x='A'; }
	else if(x == 'W'){ x='W'; }
	else if(x == 'S'){ x='S'; }
	else if(x == 'N'){ x='N'; }
	else if(x == 'K'){ x='M'; }
	else if(x == 'M'){ x='K'; }
	else if(x == 'R'){ x='Y'; }
	else if(x == 'Y'){ x='R'; }
	else if(x == 'T'){ x='A'; }
	else if(x == 'D'){ x='H'; }
	else if(x == 'H'){ x='D'; }
	else if(x == 'B'){ x='V'; }
	else if(x == 'V'){ x='B'; }
	}))),collapse="")

	return(rc)
}

#wrapper to disect a region sequence in possible 5,6,7mes containing the disered central bp
dissect.sequence <- function(regionseq){

	list5 <- NULL
	list6 <- NULL
	list7 <- NULL

	for( i in c(1:(nchar(regionseq)-6))){	#split in all possible 7mers
		list7 <- c(list7, substr(regionseq, i, (i+6)))
	}
	for( i in c(1:(nchar(regionseq)-5))){	#split in all possible 6mers
		list6 <- c(list6, substr(regionseq, i, (i+5)))
	}
	for( i in c(1:(nchar(regionseq)-4))){	#split in all possible 5mers
		list5 <- c(list5, substr(regionseq, i, (i+4)))
	}
	
	newlist <- list("list5"=list5, "list6"=list6, "list7"=list7)
	return(newlist)
}


#heavy smooth correlation
hscorr <- function(kmer, profile, profile.naked, bandwidth, plotflag){
	
	pp <- ksmooth(c(1:length(profile)), profile, kernel="normal", bandwidth=bandwidth)$y	#smooth
	pn <- ksmooth(c(1:length(profile.naked)), profile.naked, kernel="normal", bandwidth=bandwidth)$y

	hscorrelation <- cor(pp, pn)	#calc correlation
	
	if(plotflag == 0){
		newlist <- list("hscor"=hscorrelation)
		return(newlist)
	}else if(plotflag == 1){
	#assemble df
	df <- data.frame(x=c(1:length(pp)), value=c(pp,pn), species=c(rep("data profile", length(pp)),rep("background profile", length(pn))))

	theplot <- ggplot(df, aes(x=x, y=value , group=species, colour=species),  xlab="bp idx", ylab="cut probability", main=paste0("Heavy Smoothed Profiles ",kmer)) + geom_line( ) + coord_cartesian(ylim=c(0,0.01))
	newlist <- list("plot"=theplot, "hscor"=hscorrelation)
	return(newlist)
	}

}

##############################
###wrapper for tripple plot###
##############################
tri.plot <- function(kmer, plot.list, tofile=FALSE, plot.dir){
	if(tofile){
		png.file=paste0(plot.dir,'/',kmer,'_plus.png')
		png(png.file, width=900, height=900)
			grid.newpage() # Open a new page on grid device
			pushViewport(viewport(layout = grid.layout(3, 1)))
			print( plot.list[[kmer]]$plots$raw.plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))	#real
			print( plot.list[[kmer]]$plots$naked.plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)) #naked
			print( plot.list[[kmer]]$plots$norm.plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 1)) #window 
		dev.off()
	}else{
		grid.newpage() # Open a new page on grid device
		pushViewport(viewport(layout = grid.layout(3, 1)))
		print( plot.list[[kmer]]$plots$raw.plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))	#real
		print( plot.list[[kmer]]$plots$naked.plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)) #naked
		print( plot.list[[kmer]]$plots$norm.plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 1)) #window normalized
	}
}

#####################################################################
###function to query sor and store a kmer from the jaspar database###
#####################################################################
query.jaspar <- function(kmer, matrixobject, min.score){

	j.res <- searchSeq(matrixobject, kmer, strand="*", min.score=min.score)
	temp <- writeGFF3(j.res, scoreType="relative")
	temp <- temp[order(temp[,"score"], decreasing=TRUE),]
	temp <- temp[c(3,4,5,6,7,9)]
	return(temp)

}

###WRAPPER FUNCTION TO RETRIEVE RATIOS FROM A GIVEN KMER AND THE ACCORDIGN LISTS###
get.comp <- function(kmer, infile.base.comp){
	
	kl=nchar(kmer)
	kmer.list <- decode.kmer(kmer)

	t <- g <- c <- a <- rep(0, (149*2 + kl)) 
	occur.count <- 0
	for(km in kmer.list){

		match <- grep(km, readLines(infile.base.comp), value=TRUE)
		split <- strsplit(match, "\t")[[1]]	#split on \t and unlist
		a <- a + as.numeric( unlist( lapply(split[c(4:(length(split)-1))], function(x){ x <- strsplit(as.character(x), ":")[[1]][1] }) ))
		c <- c + as.numeric( unlist( lapply(split[c(4:(length(split)-1))], function(x){ x <- strsplit(as.character(x), ":")[[1]][2] }) ))
		g <- g + as.numeric( unlist( lapply(split[c(4:(length(split)-1))], function(x){ x <- strsplit(as.character(x), ":")[[1]][3] }) ))
		t <- t + as.numeric( unlist( lapply(split[c(4:(length(split)-1))], function(x){ x <- strsplit(as.character(x), ":")[[1]][4] }) ))
		occur.count <- as.numeric(split[2])

	}
	#calc ratios
	for(i in c(1:length(a))){
		summe <- sum(a[i],c[i],g[i],t[i])
		a[i] <- a[i] / summe
		c[i] <- c[i] / summe
		g[i] <- g[i] / summe
		t[i] <- t[i] / summe
	}
	#name
	names(a) <- c(1:length(a))
	names(c) <- c(1:length(c))
	names(g) <- c(1:length(g))
	names(t) <- c(1:length(t))

	newlist <-list("a"=a,"c"=c,"g"=g,"t"=t, "occur.count"=occur.count)
	return(newlist)
}

#median squared error for base comp	table.count table ratio tables have to be lready loaded
msq.base <- function(kmer, a, c, g, t){ 

	kl=nchar(kmer)

	#take medians across regions
	median.a <- median(a)
	median.c <- median(c)
	median.g <- median(g)
	median.t <- median(t)

	#calculate difference to median
	err.a <- abs(a - median.a)
 	err.c <- abs(c - median.c)
	err.g <- abs(g - median.g)
	err.t <- abs(t - median.t)
	
	err.total <- err.a + err.c + err.g + err.t #sum up
	err.sum <- sum(err.total) #sum up sums

	if(err.sum == 0){ 
		newlist <-list("err.sum"=NA, "err.iqr"=NA) #"err.profile"=NA,
		return(newlist) 
	}

	names(err.total) <- c(c(-149:-1),c(1:149))
	#calc disperion via IQR like measure
	perc.err.sum <- err.sum/4 #calc 25 % / 50 % of total variation
	i=0	
	check.sum <- 0
	while( check.sum <= perc.err.sum ){	

		check.sum <- sum( err.total[c((149-i):(150+i))] )
		i <- i+1		

	}
	err.iqr <- i	#get es range
	#err.total <- ksmooth(c(1:length(err.total)), err.total, kernel="normal", bandwidth=5)$y #smooth
	newlist <-list("err.sum"=err.sum, "err.iqr"=err.iqr) #"err.profile"=err.total,
	return(newlist)

}


##WRAPPER MEDIAN SQERROR AND IQR LIKE FOR DNASE
msq.dnase <- function(kmer, profile, fraction, scale.thresh = 0.9, max.ext=5){	#combinedlist has to be loaded
	
	kl=nchar(kmer)

	profile <- profile[-c(126:(125+kl))]#cut out kmer sequence positions

	names(profile) <- c(c(-125:-1),c(1:125))	#rename as flanking regions

	p.median <- median(profile)	#get median of profile
	err.p <- abs(profile - p.median)	#calc abs difference to median in region
	err.sum <- sum(err.p) #sum up sums

	perc.err.sum <- err.sum/fraction #25%
	i=0	
	check.sum <- 0
	while( check.sum <= perc.err.sum ){	#extend until error is smaller then the chosen scale thresh or max extension has been reached

		check.sum <- sum( err.p[c((125-i):(126+i))] )
		i <- i+1		

	}
	err.iqr <- i	#get es range
	  
	skip.flag <- 0
	extension <- 1
	#check if extending the scaling window would effectivly reduce the MSQ
	if( err.sum >= scale.thresh ){
		err.sum.temp <- err.sum
		ext <- 0
		while( err.sum.temp >= scale.thresh ){
			
			ext <- ext+1
			if(ext >= (max.ext+1)) { skip.flag = 2; break; }#break and give out not normalizable only try it till 5 bo extent to each side

			profile <- profile [-c((126-ext),(126-ext+1))]		#exclude 2 more bases from profile
			p.median <- median(profile)
			err.p <- abs(profile - p.median)
			err.sum.temp <- sum(err.p) #sum up sums	

			if(err.sum.temp < scale.thresh) { skip.flag=1; extension=ext;  break;} #set skip flag 1 means extend scale window 
		
		}

	}	

	#err.p <- ksmooth(c(1:length(err.p)), err.p, kernel="normal", bandwidth=5)$y #smooth
	newlist <-list("err.sum"=err.sum, "err.iqr"=err.iqr, "skip.flag"=skip.flag, "extension"=extension) #"err.profile"=err.p,
	return(newlist)

}

##WRAPPER FUNCTION FOR PRINTING BASE COMPOSITION
print.base.comp <- function(kmer, bandwidth, a, c, g, t){
	
	kl=nchar(kmer)	

	a <- a[-c(150:(149+kl))]
	c <- c[-c(150:(149+kl))]
	g <- g[-c(150:(149+kl))]
	t <- t[-c(150:(149+kl))]

	a <- ksmooth(c(1:length(a)), a, kernel="normal", bandwidth=bandwidth)$y
	c <- ksmooth(c(1:length(c)), c, kernel="normal", bandwidth=bandwidth)$y
	g <- ksmooth(c(1:length(g)), g, kernel="normal", bandwidth=bandwidth)$y
	t <- ksmooth(c(1:length(t)), t, kernel="normal", bandwidth=bandwidth)$y

	#assemble in df
	df <- data.frame(base=c(rep("a",length(a)), rep("c",length(c)), rep("g",length(g)), rep("t",length(t))), idx=c(c(-149:-1),c(1:149)), ratio=c(a,c,g,t) )

	#plot
	pp <- ggplot(df, aes(x=idx, y=ratio, group=base, colour=base),  ylab="base frequency",  xlab="positon relative to kmer") + geom_line(size=0.7) + geom_vline(xintercept = 0, colour="black") + coord_cartesian(ylim=c(0.1,0.6)) + 
	scale_colour_manual(values=c("green", "blue", "yellow","red")) 
	
	return(pp)

}


######################################################
### PLOT FUNCTION INCLUDING COUNT OF OCCURENCE     ###
######################################################
plot.profile.raw.data <- function(kmer, profile.plus, profile.minus, count.plus, count.minus, smooth.flag){

	kl=nchar(kmer)
	
	#smooth (if plot should be visualized without smoothing store unsmmothed profiles in variable
	if(smooth.flag == 0){
		plotting.profile.plus <- profile.plus
		plotting.profile.minus <- profile.minus
		plotting.profile.plus <- plotting.profile.plus/sum(plotting.profile.plus)
		plotting.profile.minus <- plotting.profile.minus/sum(plotting.profile.minus)
	}	

	#smooth
	profile.plus <- ksmooth(c(1:length(profile.plus)), profile.plus, kernel="normal", bandwidth=5)$y
	profile.minus <- ksmooth(c(1:length(profile.minus)), profile.minus, kernel="normal", bandwidth=5)$y
	#calculate cut probability in that region
	profile.plot.plus <- profile.plus/sum(profile.plus)
	profile.plot.minus <- profile.minus/sum(profile.minus)

#estimatefootprint and plot in output dir if plot flag set TRUE
	list.plus <- sobel.borders(profile=profile.plot.plus,kl=kl)
	list.minus <- sobel.borders(profile=profile.plot.minus,kl=kl)
#cover if no border could be estimated
	border.flag.plus <- list.plus$flag
	border.flag.minus <- list.minus$flag
#get best shoulder ranges
	us.range.plus <- list.plus$range.us
	ds.range.plus <- list.plus$range.ds
	us.range.minus <- list.minus$range.us
	ds.range.minus <- list.minus$range.ds

	#if smooth ON set profile for plotting as the smoothed norm profile
	if(smooth.flag == 0){
		profile.plot.plus <- plotting.profile.plus
		profile.plot.minus <- plotting.profile.minus 
	}

##PLOTTING
	if(border.flag.plus){	#border could be estimated

		us.border.plus <- unlist(list.plus["us"])	#get borders
		ds.border.plus <- unlist(list.plus["ds"])
#plot		
		pp <- qplot(x=c(1:length(profile.plot.plus)), y=profile.plot.plus, geom='line', xlab="bp index", ylab="cut probability", main=paste0("Cut prob. profile ",kmer,' #',count.plus)) + geom_vline(xintercept = c((us.border.plus+(us.range.plus/2)),ds.border.plus-(ds.range.plus/2)), colour="red") + geom_vline(xintercept = c((us.border.plus-(us.range.plus/2)),(ds.border.plus+(ds.range.plus/2))), colour="red", linetype = "longdash") + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))
			
	}else{		#border could not be estimated

		pp <- qplot(x=c(1:length(profile.plot.plus)),y=profile.plot.plus,geom='line', xlab="bp index", ylab="cut probability",main=paste0("Cut prob. profile ",kmer)) + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))
  
	}

	if(border.flag.minus){	#border could be estimated

		us.border.minus <- unlist(list.minus["us"])	#get borders
		ds.border.minus <- unlist(list.minus["ds"])
#unlist borders from the function (get.fp) returned list #when using sobel add 5bp to estimate the peaks
#plot		
		pm <- qplot(x=c(1:length(profile.plot.minus)),y=profile.plot.minus,geom='line', xlab="bp index", ylab="cut probability",main=paste0("Cut prob. profile ",kmer)) + geom_vline(xintercept = c((us.border.minus+(us.range.minus/2)),(ds.border.minus-(ds.range.minus/2))), colour="red") + geom_vline(xintercept = c((us.border.minus-(us.range.minus/2)),(kl-1+ds.border.minus+(ds.range.minus/2))), colour="red", linetype = "longdash") + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))

	}else{		#border could not be estimated

		pm <- qplot(x=c(1:length(profile.plot.minus)),y=profile.plot.minus,geom='line', xlab="bp index", ylab="cut probability",main=paste0("Cut prob. profile ",kmer)) + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))

	}

	newlist <- list( "pp"=pp, "pm"=pm)
	return(newlist)

}

######################################################
### PLOT FUNCTION INCLUDING COUNT OF OCCURENCE     ###
######################################################
plot.profile.background.cutting <- function(kmer, profile.plus, profile.minus, count.plus, count.minus, smooth.flag){

	kl=nchar(kmer)
	
	#smooth (if plot should be visualized without smoothing store unsmmothed profiles in variable
	if(smooth.flag == 0){
		plotting.profile.plus <- profile.plus
		plotting.profile.minus <- profile.minus
		plotting.profile.plus <- plotting.profile.plus/sum(plotting.profile.plus)
		plotting.profile.minus <- plotting.profile.minus/sum(plotting.profile.minus)
	}	
	profile.plus <- ksmooth(c(1:length(profile.plus)), profile.plus, kernel="normal", bandwidth=5)$y
	profile.minus <- ksmooth(c(1:length(profile.minus)), profile.minus, kernel="normal", bandwidth=5)$y
	#calculate cut probability in that region
	profile.plot.plus <- profile.plus/sum(profile.plus)
	profile.plot.minus <- profile.minus/sum(profile.minus)

#estimatefootprint and plot in output dir if lot flag set TRUE
	list.plus <- sobel.borders(profile=profile.plot.plus,kl=kl)
	list.minus <- sobel.borders(profile=profile.plot.minus,kl=kl)
#cover if no border could be estimated
	border.flag.plus <- list.plus$flag
	border.flag.minus <- list.minus$flag
#get best shoulder ranges
	us.range.plus <- list.plus$range.us
	ds.range.plus <- list.plus$range.ds
	us.range.minus <- list.minus$range.us
	ds.range.minus <- list.minus$range.ds

	#if to smooth for plot set profile sfor plotting to the smoothed profiles
	if(smooth.flag == 1){
		plotting.profile.plus <- profile.plot.plus
		plotting.profile.minus <- profile.plot.minus
		plotting.profile.plus <- plotting.profile.plus/sum(plotting.profile.plus)
		plotting.profile.minus <- plotting.profile.minus/sum(plotting.profile.minus)
	}

##PLOTTING
	pp <- qplot(x=c(1:length(plotting.profile.plus)),y=plotting.profile.plus, geom='line', xlab="bp index", ylab="cut probability",main=paste0("Cut prob. profile ",kmer, ' #', count.plus)) + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))
  
	pm <- qplot(x=c(1:length(plotting.profile.minus)),y=plotting.profile.minus, geom='line', xlab="bp index", ylab="cut probability",main=paste0("Cut prob. profile ",kmer, ' #', count.minus )) + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))

	newlist <- list( "pp"=pp, "pm"=pm)
	return(newlist)

}



####################################################################################################################################
### Function to create plots of gaussian smoothed average profiles, normalized to cut probability  and sequence specifc variation###
####################################################################################################################################
plot.profile.normalized <- function(kmer, profile.plot.plus, profile.plot.minus, profile.plot.plus.naked, profile.plot.minus.naked, skip.flag, extension){

	kl=nchar(kmer)

#gaussian smooth after scaled diff to median normalization
	list.plus <- gaussian.median.diff(profile.plot.plus, profile.plot.plus.naked, kl, skip.flag, extension)
	list.minus <- gaussian.median.diff(profile.plot.minus, profile.plot.minus.naked, kl, skip.flag, extension)
	profile.norm.plus <- list.plus$profile
	profile.norm.minus <- list.minus$profile
#retrieve norm.flag
	norm.flag.plus <- list.plus$norm.flag
	norm.flag.minus <- list.minus$norm.flag

#estimatefootprint and plot in output dir if lot flag set TRUE
	list.plus <- sobel.borders(profile=profile.norm.plus,kl=kl)
	list.minus <- sobel.borders(profile=profile.norm.minus,kl=kl)
#cover if no border could be estimated
	border.flag.plus <- list.plus$flag
	border.flag.minus <- list.minus$flag
#get best shoulder ranges
	us.range.plus <- list.plus$range.us
	ds.range.plus <- list.plus$range.ds
	us.range.minus <- list.minus$range.us
	ds.range.minus <- list.minus$range.ds


### IMPROVED PLOTTING PROCEDURE, INCLUDING SOBEL BORDER ESTIMATION AND WIDTH OF SHOULDER OPTIMATION ###
	if(border.flag.plus){	#border could be estimated
#unlist borders from the function (get.fp) returned list #when using sobel add 5bp to estimate the peaks

		us.border.plus <- unlist(list.plus["us"])	#get borders
		ds.border.plus <- unlist(list.plus["ds"])
#plot		
		pp <- qplot(x=c(1:length(profile.norm.plus)),y=profile.norm.plus,geom='line', xlab="bp index", ylab="cut probability",main=paste0("Cut prob. profile ",kmer, ' norm.flag = ',norm.flag.plus)) + geom_vline(xintercept = c((us.border.plus+(us.range.plus/2)),ds.border.plus-(ds.range.plus/2)), colour="red") + geom_vline(xintercept = c((us.border.plus-(us.range.plus/2)),(ds.border.plus+(ds.range.plus/2))), colour="red", linetype = "longdash") + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + geom_vline(xintercept = c((126-extension),(126+kl-1+extension)), linetype = "dotted") + coord_cartesian(ylim=c(0,0.01))
			
	}else{		#border could not be estimated

		pp <- qplot(x=c(1:length(profile.norm.plus)),y=profile.norm.plus,geom='line', xlab="bp index", ylab="cut probability",main=paste0("Cut prob. profile ",kmer)) + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + geom_vline(xintercept = c((126-extension),(126+kl-1+extension)), linetype = "dotted") + coord_cartesian(ylim=c(0,0.01)) + coord_cartesian(ylim=c(0,0.01))
  
	}

	if(border.flag.minus){	#border could be estimated

		us.border.minus <- unlist(list.minus["us"])	#get borders
		ds.border.minus <- unlist(list.minus["ds"])
#unlist borders from the function (get.fp) returned list #when using sobel add 5bp to estimate the peaks
#plot		
		pm <- qplot(x=c(1:length(profile.norm.minus)),y=profile.norm.minus,geom='line', xlab="bp index", ylab="cut probability",main=paste0("Cut prob. profile ",kmer, ' norm.flag = ',norm.flag.minus)) + geom_vline(xintercept = c((us.border.minus+(us.range.minus/2)),(ds.border.minus-(ds.range.minus/2))), colour="red") + geom_vline(xintercept = c((us.border.minus-(us.range.minus/2)),(ds.border.minus+(ds.range.minus/2))), colour="red", linetype = "longdash") + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + geom_vline(xintercept = c((126-extension),(126+kl-1+extension)), linetype = "dotted") + coord_cartesian(ylim=c(0,0.01))

	}else{		#border could not be estimated

		pm <- qplot(x=c(1:length(profile.norm.minus)),y=profile.norm.minus,geom='line', ylab="cut probability", xlab="bp index", main=paste0("Cut prob. profile ",kmer)) + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + geom_vline(xintercept = c((126-extension),(126+kl-1+extension)), linetype = "dotted") + coord_cartesian(ylim=c(0,0.01))

	}

	newlist <- list("pp"=pp, "pm"=pm)
	return(newlist)

}

#######################################################################
### Wraper Function to retrieve profile and Counts from combined.list ###
#######################################################################
retrieve.profile.grep <- function(kmer, infile.p, infile.m){

	kmer.list <- decode.kmer(kmer)
	kl <- nchar(kmer)
	count.plus <- 0
	count.minus <- 0
	profile.plus <- rep(0,(250+kl))
	profile.minus <- rep(0,(250+kl))

	for(km in kmer.list){

		#plus
		match <- grep(km, readLines(infile.p), value=TRUE)
		split <- strsplit(match, "\t")[[1]]	#split on \t and unlist
		count.plus <- count.plus + as.numeric(split[2])
		profile.plus <- profile.plus + as.numeric(split[c(28:(length(split)-25))])	#trimm as before with list -25 on front + end

		#minus
		match <- grep(km, readLines(infile.m), value=TRUE)
		split <- strsplit(match, "\t")[[1]]	#split on \t and unlist
		count.minus <- count.minus + as.numeric(split[2])
		profile.minus <- profile.minus + as.numeric(split[c(28:(length(split)-25))])	#trimm as before with list -25 on front + end		

	}
	profile.plot.plus <- profile.plus/sum(profile.plus)	
	profile.plot.minus <- profile.minus/sum(profile.minus)

	newlist <- list("profile.plus"=profile.plot.plus, "profile.minus"=profile.plot.minus, "count.plus"=count.plus, "count.minus"=count.minus)
	return(newlist)
}

##################################################################
### WRAPPER FUNCTION FOR THE WHOLE PROCESS INCLUDING ALL PLOTS ###
##################################################################
full.work.withplot <- function(kmer, infile.plus, infile.minus, infile.naked.plus, infile.naked.minus, infile.base.comp, smooth.flag){

	#get kmer length and retrieve the profiles from stored objects
	kl=nchar(kmer)
	profiles <- retrieve.profile.grep(kmer, infile.plus, infile.minus)
	profiles.naked <- retrieve.profile.grep(kmer, infile.naked.plus, infile.naked.minus)

	#get HeavySmoothCorrelation value and plot
	hsc.list <- hscorr(kmer, profiles$profile.plus, profiles.naked$profile.plus, bandwidth=10, plotflag=1)

	#get MSQ and IQR like measures for DNase Cutting on background
	msq.dnase.list <- msq.dnase(kmer, profiles.naked$profile.plus, fraction=4, scale.thresh=0.09, max.ext=5)	#fraction 4 ~ 25% IQR
	#get base composition stats and plot
	get.base.list <- get.comp(kmer, infile.base.comp)
	msq.base.list <- msq.base(kmer, get.base.list$a, get.base.list$c, get.base.list$g, get.base.list$t)	#get msq msqprofile and iqr
	base.comp.plot <- print.base.comp(kmer=kmer, bandwidth=3, a=get.base.list$a, c=get.base.list$c, g=get.base.list$g, t=get.base.list$t)

	#get smooted data, smoothed naked and normalized and smoothed profiles
	profile.raw.plot.list <- plot.profile.raw.data(kmer, profile.plus=profiles$profile.plus, profile.minus=profiles$profile.minus, count.plus=profiles$count.plus, count.minus=profiles$count.minus, smooth.flag)
	
	profile.naked.plot.list <- plot.profile.background.cutting(kmer, profile.plus=profiles.naked$profile.plus, profile.minus=profiles.naked$profile.minus, count.plus=profiles.naked$count.plus, count.minus=profiles.naked$count.minus, smooth.flag)

	profile.norm.plot.list <- plot.profile.normalized(kmer, profiles$profile.plus, profiles$profile.minus, profiles.naked$profile.plus, profiles.naked$profile.minus, msq.dnase.list$skip.flag, msq.dnase.list$extension)

	#gaussian smooth after scaled diff to median normalization
	gauss.list.plus <- gaussian.median.diff(profiles$profile.plus, profiles.naked$profile.plus, kl, msq.dnase.list$skip.flag, msq.dnase.list$extension)
	gauss.list.minus <- gaussian.median.diff(profiles$profile.minus, profiles.naked$profile.minus, kl, msq.dnase.list$skip.flag, msq.dnase.list$extension)
	#estimate footprint, borders and best shoulder ranges
	border.list.plus <- sobel.borders(profile=gauss.list.plus$profile, kl=kl)
	border.list.minus <- sobel.borders(profile=gauss.list.minus$profile, kl=kl)
	#calculate adjusted FSR
	if(border.list.plus$flag){

		fsr.plus <- adj.fsr(profile=gauss.list.plus$profile, kl=kl, usb=border.list.plus$us, dsb=border.list.plus$ds, usr=border.list.plus$range.us, dsr=border.list.plus$range.ds)

	}else{	
		fsr.plus="NA" 
	}

	if(border.list.minus$flag){

		fsr.minus <- adj.fsr(profile=gauss.list.minus$profile, kl=kl, usb=border.list.minus$us, dsb=border.list.minus$ds, usr=border.list.minus$range.us, dsr=border.list.plus$range.ds)

	}else{	
		fsr.minus="NA" 
	}


	#assemble statslist
	kmer.stat.list <- list(
		"fsr.plus" = fsr.plus, "fsr.minus" = fsr.minus,
		"scale.plus" = gauss.list.plus$scale.factor, "scale.minus" = gauss.list.minus$scale.factor,
		"hscor" = hsc.list$hscor,
		"dnase.msq" = msq.dnase.list$err.sum, "dnase.iqr" = msq.dnase.list$err.iqr,	
		"base.msq" = msq.base.list$err.sum, "base.iqr" = msq.base.list$err.iqr,
		"norm.flag" = msq.dnase.list$skip.flag,
		"extension" = msq.dnase.list$extension
	
	)
	#assemble plot.list
	kmer.plot.list <- list(
		"raw.plot.plus" = profile.raw.plot.list$pp,
		"raw.plot.minus" = profile.raw.plot.list$pm,
		"naked.plot.plus" = profile.naked.plot.list$pp,
		"naked.plot.minus" = profile.naked.plot.list$pm,
		"norm.plot.plus" = profile.norm.plot.list$pp,
		"norm.plot.minus" = profile.norm.plot.list$pm,
		"base.plot" = base.comp.plot,
		"hsc.plot" = hsc.list$plot
	)

	newlist <- list("stats"=kmer.stat.list, "plots"=kmer.plot.list)
	return(newlist)
}


##################################################################
### WRAPPER FUNCTION FOR THE WHOLE PROCESS EXCLUDING ALL PLOTS ###
##################################################################
full.work.noplot <- function(kmer, infile.plus, infile.minus, infile.naked.plus, infile.naked.minus, infile.base.comp){

	#get kmer length and retrieve the profiles from stored objects
	kl=nchar(kmer)
	profiles <- retrieve.profile.grep(kmer, infile.plus, infile.minus)
	profiles.naked <- retrieve.profile.grep(kmer, infile.naked.plus, infile.naked.minus)

	#get HeavySmoothCorrelation value and plot
	hsc.list <- hscorr(kmer, profiles$profile.plus, profiles.naked$profile.plus, bandwidth=10, plotflag=0)

	#get MSQ and IQR like measures for DNase Cutting on background
	msq.dnase.list <- msq.dnase(kmer, profiles.naked$profile.plus, fraction=4, scale.thresh=0.09, max.ext=5)	#fraction 4 ~ 25% IQR
	#get base composition stats and plot
	get.base.list <- get.comp(kmer, infile.base.comp)
	msq.base.list <- msq.base(kmer, get.base.list$a, get.base.list$c, get.base.list$g, get.base.list$t)	#get msq msqprofile and iqr
	
	#gaussian smooth after scaled diff to median normalization
	gauss.list.plus <- gaussian.median.diff(profiles$profile.plus, profiles.naked$profile.plus, kl, msq.dnase.list$skip.flag, msq.dnase.list$extension)
	gauss.list.minus <- gaussian.median.diff(profiles$profile.minus, profiles.naked$profile.minus, kl, msq.dnase.list$skip.flag, msq.dnase.list$extension)
	#estimate footprint, borders and best shoulder ranges
	border.list.plus <- sobel.borders(profile=gauss.list.plus$profile, kl=kl)
	border.list.minus <- sobel.borders(profile=gauss.list.minus$profile, kl=kl)
	#calculate adjusted FSR
	if(border.list.plus$flag){

		fsr.plus <- adj.fsr(profile=gauss.list.plus$profile, kl=kl, usb=border.list.plus$us, dsb=border.list.plus$ds, usr=border.list.plus$range.us, dsr=border.list.plus$range.ds)

	}else{	
		fsr.plus="NA" 
	}

	if(border.list.minus$flag){

		fsr.minus <- adj.fsr(profile=gauss.list.minus$profile, kl=kl, usb=border.list.minus$us, dsb=border.list.minus$ds, usr=border.list.minus$range.us, dsr=border.list.plus$range.ds)

	}else{	
	fsr.minus="NA" 
	}


	#assemble statslist
	kmer.stat.list <- list(
		"fsr.plus" = fsr.plus, "fsr.minus" = fsr.minus,
		"scale.plus" = gauss.list.plus$scale.factor, "scale.minus" = gauss.list.minus$scale.factor,
		"hscor" = hsc.list$hscor,
		"dnase.msq" = msq.dnase.list$err.sum, "dnase.iqr" = msq.dnase.list$err.iqr,	
		"base.msq" = msq.base.list$err.sum, "base.iqr" = msq.base.list$err.iqr,
		"norm.flag" = msq.dnase.list$skip.flag,
		"extension" = msq.dnase.list$extension
	
	)
	newlist <- list("stats"=kmer.stat.list)
	return(newlist)
}


###########################################################################
### Function to create overlap plots of normalized profiles 2kmers	###
###########################################################################
plot.profile.single.norm.overlap <- function(kmer1, kmer2, profile.plus1, profile.plus2, profile.minus1, profile.minus2, profile.naked.plus1, profile.naked.plus2, profile.naked.minus1, profile.naked.minus2, skip.flag1, skip.flag2, extension1, extension2, count.plus1, count.plus2, count.minus1, count.minus2){

	kl=nchar(kmer1)

#gaussian smooth after scaled diff to median normalization
	list.plus1 <- gaussian.median.diff(profile.plus1, profile.naked.plus1, kl, skip.flag1, extension1)
	list.minus1 <- gaussian.median.diff(profile.minus1, profile.naked.minus1, kl, skip.flag1, extension1)
	profile.norm.plus1 <- list.plus1$profile
	profile.norm.minus1 <- list.minus1$profile
#retrieve norm.flag
	norm.flag.plus1 <- list.plus1$norm.flag
	norm.flag.minus1 <- list.minus1$norm.flag
#2
	list.plus2 <- gaussian.median.diff(profile.plus2, profile.naked.plus2, kl, skip.flag2, extension2)
	list.minus2 <- gaussian.median.diff(profile.minus2, profile.naked.minus2, kl, skip.flag2, extension2)
	profile.norm.plus2 <- list.plus2$profile
	profile.norm.minus2 <- list.minus2$profile
	norm.flag.plus2 <- list.plus2$norm.flag
	norm.flag.minus2 <- list.minus2$norm.flag

#make dataframe
df <- data.frame(x=c(1:length(profile.norm.plus1)), value=c(profile.norm.plus1, profile.norm.plus2), kmer=c(rep(paste0("norm profile ",kmer1), length(profile.norm.plus1)),rep(paste0("norm profile ",kmer2), length(profile.norm.plus2)) ))

#OVERLAPPLOT Plus
pp <-ggplot(df, aes(x=x, y=value , group=kmer, colour=kmer)) + geom_line() + labs(title=paste0(kmer1,' #', count.plus1 ,' norm.flag = ',norm.flag.plus1,' vs. ',kmer2, ' #', count.plus2,' norm.flag = ',norm.flag.plus2), x="bp idx", y="cut probability") + geom_vline(xintercept = c(126,(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))

#make dataframe
df <- data.frame(x=c(1:length(profile.norm.minus1)), value=c(profile.norm.minus1,profile.norm.minus2), kmer=c(rep(paste0("norm profile ",kmer1), length(profile.norm.minus1)),rep(paste0("norm profile ",kmer2), length(profile.norm.minus2)) ))

#OVERLAPPLOT Minus
pm <-ggplot(df, aes(x=x, y=value , group=kmer, colour=kmer)) + geom_line( ) + labs(title=paste0(kmer1,' #', count.minus1 ,' norm.flag = ',norm.flag.minus1,' vs. ',kmer2, ' #', count.minus2,' norm.flag = ',norm.flag.minus2), x="bp idx", y="cut probability") + geom_vline(xintercept = c(126,(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))

	newlist <- list("pp"=pp, "pm"=pm)
	return(newlist)

}

##merged variant#
plot.profile.single.norm.overlap.merge <- function(kmer1, kmer2, profile.plus1, profile.plus2, profile.minus1, profile.minus2, profile.naked.plus1, profile.naked.plus2, profile.naked.minus1, profile.naked.minus2, skip.flag1, skip.flag2, extension1, extension2, count.plus1, count.plus2, count.minus1, count.minus2){

	kl=nchar(kmer1)

#gaussian smooth after scaled diff to median normalization
	list.plus1 <- gaussian.median.diff(profile.plus1, profile.naked.plus1, kl, skip.flag1, extension1)
	list.minus1 <- gaussian.median.diff(profile.minus1, profile.naked.minus1, kl, skip.flag1, extension1)
	profile.norm.plus1 <- list.plus1$profile
	profile.norm.minus1 <- list.minus1$profile
#retrieve norm.flag
	norm.flag.plus1 <- list.plus1$norm.flag
	norm.flag.minus1 <- list.minus1$norm.flag
#2
	list.plus2 <- gaussian.median.diff(profile.plus2, profile.naked.plus2, kl, skip.flag2, extension2)
	list.minus2 <- gaussian.median.diff(profile.minus2, profile.naked.minus2, kl, skip.flag2, extension2)
	profile.norm.plus2 <- list.plus2$profile
	profile.norm.minus2 <- list.minus2$profile
	norm.flag.plus2 <- list.plus2$norm.flag
	norm.flag.minus2 <- list.minus2$norm.flag

	#merge
	profile.norm.plus1 <- ( profile.norm.plus1 + profile.norm.minus1)/2
	profile.norm.plus2 <- ( profile.norm.plus2 + profile.norm.minus2)/2


	#make dataframe
	df <- data.frame(x=c(1:length(profile.norm.plus1)), value=c(profile.norm.plus1, profile.norm.plus2), kmer=c(rep(paste0("norm profile ",kmer1), length(profile.norm.plus1)),rep(paste0("norm profile ",kmer2), length(profile.norm.plus2)) ))

	#OVERLAPPLOT Plus
	pp <-ggplot(df, aes(x=x, y=value , group=kmer, colour=kmer)) + geom_line() + labs(title=paste0(kmer1,' #', count.plus1 ,' norm.flag = ',norm.flag.plus1,' vs. ',kmer2, ' #', count.plus2,' norm.flag = ',norm.flag.plus2), x="bp idx", y="cut probability") + geom_vline(xintercept = c(126,(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))

	newlist <- list("pp"=pp)
	return(newlist)

}

###########################################################################
### Function to create overlap plots of normalized profiles 2kmers PNORM! 	###
###########################################################################
plot.profile.single.pnorm.overlap <- function(kmer1, kmer2, profile.plus1, profile.plus2, profile.minus1, profile.minus2, count.plus1, count.plus2, count.minus1, count.minus2){

	kl=nchar(kmer1)

#make dataframe
df <- data.frame(x=c(1:length(profile.plus1)), value=c(profile.plus1, profile.plus2), kmer=c(rep(paste0("norm profile ",kmer1), length(profile.plus1)),rep(paste0("norm profile ",kmer2), length(profile.plus2)) ))

#OVERLAPPLOT Plus
pp <-ggplot(df, aes(x=x, y=value , group=kmer, colour=kmer)) + geom_line() + labs(title=paste0(kmer1,' #', count.plus1,' vs. ',kmer2, ' #', count.plus2), x="bp idx", y="cut probability") + geom_vline(xintercept = c(126,(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))

#make dataframe
df <- data.frame(x=c(1:length(profile.minus1)), value=c(profile.minus1, profile.minus2), kmer=c(rep(paste0("norm profile ",kmer1), length(profile.minus1)),rep(paste0("norm profile ",kmer2), length(profile.minus2)) ))

#OVERLAPPLOT Minus
pm <-ggplot(df, aes(x=x, y=value , group=kmer, colour=kmer)) + geom_line( ) + labs(title=paste0(kmer1,' #', count.minus1 ,' vs. ',kmer2, ' #', count.minus2), x="bp idx", y="cut probability") + geom_vline(xintercept = c(126,(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))

	newlist <- list("pp"=pp, "pm"=pm)
	return(newlist)

}

## For merged ##
plot.profile.single.pnorm.overlap.merge <- function(kmer1, kmer2, profile.plus1, profile.plus2, count.plus1, count.plus2){

	kl=nchar(kmer1)

#make dataframe
df <- data.frame(x=c(1:length(profile.plus1)), value=c(profile.plus1, profile.plus2), kmer=c(rep(paste0("norm profile ",kmer1), length(profile.plus1)),rep(paste0("norm profile ",kmer2), length(profile.plus2)) ))

#OVERLAPPLOT Plus
pp <-ggplot(df, aes(x=x, y=value , group=kmer, colour=kmer)) + geom_line() + labs(title=paste0(kmer1,' #', count.plus1,' vs. ',kmer2, ' #', count.plus2), x="bp idx", y="cut probability") + geom_vline(xintercept = c(126,(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))
	
	newlist <- list("pp"=pp)
	return(newlist)

}

######################################################
### PLOT FUNCTION INCLUDING COUNT OF OCCURENCE     ###
######################################################
plot.profile.pnorm.data <- function(kmer, profile.plus, profile.minus, count.plus, count.minus, smooth.flag){

	kl=nchar(kmer)
	
	#smooth (if plot should be visualized without smoothing store unsmmothed profiles in variable
	if(smooth.flag == 0){
		plotting.profile.plus <- profile.plus
		plotting.profile.minus <- profile.minus
		plotting.profile.plus <- plotting.profile.plus/sum(plotting.profile.plus)
		plotting.profile.minus <- plotting.profile.minus/sum(plotting.profile.minus)
	}	

	#smooth
	profile.plus <- ksmooth(c(1:length(profile.plus)), profile.plus, kernel="normal", bandwidth=5)$y
	profile.minus <- ksmooth(c(1:length(profile.minus)), profile.minus, kernel="normal", bandwidth=5)$y
	#calculate cut probability in that region
	profile.plot.plus <- profile.plus/sum(profile.plus)
	profile.plot.minus <- profile.minus/sum(profile.minus)

#estimatefootprint and plot in output dir if plot flag set TRUE
	list.plus <- sobel.borders(profile=profile.plot.plus,kl=kl)
	list.minus <- sobel.borders(profile=profile.plot.minus,kl=kl)
#cover if no border could be estimated
	border.flag.plus <- list.plus$flag
	border.flag.minus <- list.minus$flag
#get best shoulder ranges
	us.range.plus <- list.plus$range.us
	ds.range.plus <- list.plus$range.ds
	us.range.minus <- list.minus$range.us
	ds.range.minus <- list.minus$range.ds

	#if smooth ON set profile for plotting as the smoothed norm profile
	if(smooth.flag == 0){
		profile.plot.plus <- plotting.profile.plus
		profile.plot.minus <- plotting.profile.minus 
	}

##PLOTTING
	if(border.flag.plus){	#border could be estimated

		us.border.plus <- unlist(list.plus["us"])	#get borders
		ds.border.plus <- unlist(list.plus["ds"])
#plot		
		pp <- qplot(x=c(1:length(profile.plot.plus)), y=profile.plot.plus, geom='line', xlab="bp index", ylab="cut probability", main=paste0("Norm Cut prob. profile ",kmer,' #',count.plus)) + geom_vline(xintercept = c((us.border.plus+(us.range.plus/2)),ds.border.plus-(ds.range.plus/2)), colour="red") + geom_vline(xintercept = c((us.border.plus-(us.range.plus/2)),(ds.border.plus+(ds.range.plus/2))), colour="red", linetype = "longdash") + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))
			
	}else{		#border could not be estimated

		pp <- qplot(x=c(1:length(profile.plot.plus)),y=profile.plot.plus,geom='line', xlab="bp index", ylab="cut probability",main=paste0("Norm Cut prob. profile ",kmer)) + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))
  
	}

	if(border.flag.minus){	#border could be estimated

		us.border.minus <- unlist(list.minus["us"])	#get borders
		ds.border.minus <- unlist(list.minus["ds"])
#unlist borders from the function (get.fp) returned list #when using sobel add 5bp to estimate the peaks
#plot		
		pm <- qplot(x=c(1:length(profile.plot.minus)),y=profile.plot.minus,geom='line', xlab="bp index", ylab="cut probability",main=paste0("Norm Cut prob. profile ",kmer)) + geom_vline(xintercept = c((us.border.minus+(us.range.minus/2)),(ds.border.minus-(ds.range.minus/2))), colour="red") + geom_vline(xintercept = c((us.border.minus-(us.range.minus/2)),(ds.border.minus+(ds.range.minus/2))), colour="red", linetype = "longdash") + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))

	}else{		#border could not be estimated

		pm <- qplot(x=c(1:length(profile.plot.minus)),y=profile.plot.minus,geom='line', xlab="bp index", ylab="cut probability",main=paste0("Norm Cut prob. profile ",kmer)) + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))

	}

	newlist <- list( "pp"=pp, "pm"=pm)
	return(newlist)

}

############################
### PLOT PNORM MERGED    ###
############################
plot.profile.pnorm.data.single <- function(kmer, profile.plus, count.plus, smooth.flag){

	kl=nchar(kmer)
	
	#smooth (if plot should be visualized without smoothing store unsmmothed profiles in variable
	if(smooth.flag == 0){
		plotting.profile.plus <- profile.plus
		plotting.profile.plus <- plotting.profile.plus/sum(plotting.profile.plus)
	}	

	#smooth
	profile.plus <- ksmooth(c(1:length(profile.plus)), profile.plus, kernel="normal", bandwidth=5)$y
	#calculate cut probability in that region
	profile.plot.plus <- profile.plus/sum(profile.plus)
	
#estimatefootprint and plot in output dir if plot flag set TRUE
	list.plus <- sobel.borders(profile=profile.plot.plus,kl=kl)
	#cover if no border could be estimated
	border.flag.plus <- list.plus$flag
	#get best shoulder ranges
	us.range.plus <- list.plus$range.us
	ds.range.plus <- list.plus$range.ds
	
	#if smooth ON set profile for plotting as the smoothed norm profile
	if(smooth.flag == 0){
		profile.plot.plus <- plotting.profile.plus
	}

##PLOTTING
	if(border.flag.plus){	#border could be estimated

		us.border.plus <- unlist(list.plus["us"])	#get borders
		ds.border.plus <- unlist(list.plus["ds"])
#plot		
		pp <- qplot(x=c(1:length(profile.plot.plus)), y=profile.plot.plus, geom='line', xlab="bp index", ylab="cut probability", main=paste0("Norm Cut prob. profile ",kmer,' #',count.plus)) + geom_vline(xintercept = c((us.border.plus+(us.range.plus/2)),ds.border.plus-(ds.range.plus/2)), colour="red") + geom_vline(xintercept = c((us.border.plus-(us.range.plus/2)),(ds.border.plus+(ds.range.plus/2))), colour="red", linetype = "longdash") + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))
			
	}else{		#border could not be estimated

		pp <- qplot(x=c(1:length(profile.plot.plus)),y=profile.plot.plus,geom='line', xlab="bp index", ylab="cut probability",main=paste0("Norm Cut prob. profile ",kmer)) + geom_vline(xintercept = c((126),(126+kl-1)), linetype = "longdash") + coord_cartesian(ylim=c(0,0.01))
  
	}

	newlist <- list( "pp"=pp )
	return(newlist)

}

### ================================================= ###
### Second Batch of wrapper functions for R usage
### ================================================= ###

#function to retrieve kmer profile and count
grep.profile <- function(kmer, infile){
  
  kmer.list <- decode.kmer(kmer)  #function that decodes ambivalent fasta code  letters like W into a list of all possible kmers ( if W two kmers to acutally search for)
  count <- 0 #initialize variables
  kl <- nchar(kmer)
  profile <- rep(0,(250+kl))
  
  for(km in kmer.list){  #for all entries in the splitted up kmer list
    match <- grep(km, readLines(infile), value=TRUE)  #performs basically a system grep / not sure how fast in relation to systemcommand
    split <- strsplit(match, "\t")[[1]]	#split on "\t" and unlist, all flat files are currently saved with a tab as seperator
    count <- count + as.numeric(split[2])	#second entry (first number) is the count how often the kmer appeared in the searched regions, rest is the cutting profile
    #remove kmer and count from split 
    split <- split[-c(1,2)]
    split <- as.numeric(split)
    profile <- profile + split[c(26:(length(split)-25))]	#get length of split
  }
  profile <- (profile / sum(profile))
 
  newlist <- list( "profile"=profile, "count"=count ) #assemble return list
  return(newlist)
  
}

#function to retrieve a normalized and pruned footprint (load completely prior working)
get.footprint <- function(kmer,tissue, data.dir, length, smooth){
  
  infile.plus=file.path(data.dir,tissue,"counts",paste0("kmers_",kl,"_count_",tissue,"_pnorm_JH60_plus.txt"))
  infile.minus=file.path(data.dir,tissue,"counts",paste0("kmers_",kl,"_count_",tissue,"_pnorm_JH60_minus.txt"))
  
  #grep profiles & counts
  l.plus <- grep.profile(kmer, infile.plus)
  l.minus <- grep.profile(kmer, infile.minus)
  
  if(frag.type == "ATAC"){
    profile.merge <- (l.plus$profile + l.minus$profile)/2  #sum up profile, average is calculated later anyway
  }else if(frag.type == "DNase"){
    middle <- ( l.plus$profile[c(126:(125+kl))] + l.minus$profile[c(126:(125+kl))] ) / 2 #calc average of both strands for the kmer length middle
    profile.merge <- c( l.plus$profile[c(1:125)], middle, l.minus$profile[c((125+kl+1):length(l.minus$profile))] )	#combine left flank from plus right flank from minus and kmer middle from the average
  }else{
    cat("Select according $FRAG_TYPE")
  }
  
  #smooth if specified
  if(smooth){
    profile.merge <- ksmooth(c(1:length(profile.merge)), profile.merge, kernel="normal", bandwidth=5)$y
  }
  
  remove.temp <- (250-length)/2
  profile.merge <- profile.merge[c((1+remove.temp)):(length(profile.merge)-remove.temp)]

  return(list(profile=profile.merge, count=l.plus$count))
  
}


