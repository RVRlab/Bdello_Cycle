#######################################################################################################
########### Disclaimer
# - This functions were created by Fares Osam Yañez Cuna while working in the Leach lab in 2020
# - So far, this work has not been published somewhere. 
# - It is therefore not ok to distribute either the full code or some of its parts.
#######################################################################################################
# - This is a list of all functions I use for reading, analyzig and graphing data for the MFA Analysis.
# - This functions were created to make my R scripts easier to read/understand.
# - Anything that had to be done twice in my script was made into a function.
# - Anything that could be encapsulated, was done so I dont have to touch the original data
# - I know some variable decisions might not be the most efficient, but they are working, are quick
# and allowed me to easily test them whenever I had a mistake. 
# - Any question ask: fyanez@ed.ac.uk, faresosam@gmail.com  
#######################################################################################################


#++++++++++++++++++++++++++++++++ Functions ++++++++++++++++++++++++++++++++++++++++++++++++
#----------- Basic functions -----------

# This function is for getting the genome length size of a .fasta file
get_genome_length <- function(reference_genome_path){
	genome <- read.fasta(reference_genome_path)[[1]]
	genome_seq <- paste(genome[1:length(genome)], collapse = "")
	genome_length <- length(genome)

	return(genome_length)
}

# This function is to read the given depth files and merge them into one raw data dataframe
filesformiles <- function(files,genome_length){
  data <- data.frame(position=seq(from=1,to=genome_length))
  for(i in 1:length(files)){
    temp <- read.table(files[i])
    temp[,1] <- NULL
    log_name_split <- strsplit(as.character(files[i]),"_")[[1]]
    colnames(temp) <- c("position",strsplit(log_name_split[length(log_name_split)],".txt")[[1]])
    data <- merge(data,temp,by="position",all.x=TRUE)
  }
  return(data)
}

# Simple function to divide a data point between the sum of all the values of a colum
#**** It assumes the first column is just the position values (so it skips it)
normalizer <- function(data_frame) {
  temporal <- data_frame
  for (j in 2:ncol(data_frame)){
    temporal[,j] <- data_frame[,j]/sum(data_frame[,j],na.rm=T)
  }
  return (temporal)
}

# A function to average multiple MFA data
simpleAverager <- function (data.frame) {
  AverageableData <- data.frame
  average <- data.frame(matrix(0,ncol=1, nrow=dim(AverageableData)[1]))
  for (k in 2:ncol(AverageableData)){
    average <-average + AverageableData[,k]
  }
  average <- average/(dim(AverageableData)[2]-1)
  average <- data.frame(cbind(AverageableData[,1]),average)
  colnames(average) <- c("position","average")
  return(average)
}

# This function is for binning the depth coverage data into bins (windows) of n size (we use n=1000)
beaner <- function(genome_length,data.frame,n){
  position <- data.frame(position=seq(from=1,to=genome_length))
  data <- merge(position,data.frame,by="position",all.x=TRUE)
  BinnedData <- data.frame(matrix(NA,ncol=dim(data)[2],nrow=floor((dim(data)[1])/n)))
  colnames(BinnedData) <- colnames(data)
  for(j in seq(from=1,to=floor(length(data[,1])/n))){
    BinnedData[j,1] <- (j*n)-(n/2)
  }
  for(i in 2:dim(data)[2]){
    for(j in seq(from=1,to=floor(length(data[,i])/n)))
    {
      BinnedData[j,i] <- mean(data[((j*n)-(n-1)):(j*n),i], na.rm = T)
    }
  }
  return(BinnedData)
}

# This function is for making the loess lines for all vectors in a set of data
linemaker <- function(data.frame,s){
  loess_data <- data.frame
  for(i in 2:ncol(loess_data)){
    loess_data[,i] <- NA
  }
  for(i in 2:ncol(loess_data)){
    y.loess <- loess(y~x, span=s, data.frame(x=data.frame[,1],y=data.frame[,i]))
    loess_data[,i] <- predict(y.loess, data.frame(x=data.frame[,1]))
  }
  return(loess_data)
}

# A function for doing the rolling average of a genome to smooth outliers out of genomic data
RollingAverage <- function(data, rollEvery, windowSize) {
  if(dim(data)[1] <= windowSize) {
      stop('The data you provided is not big enough for doing the analysis')
  }
  ## Generating the result array
  windows <- ((dim(data)[1] - windowSize)%/%rollEvery)+1 # we had to add the +1 as we are adding an extra value when considering the final 5000 bp / 2 on each side
  result <- data.frame(matrix(0,ncol=3, nrow=windows))
  colnames(result) <- c("position","average", "SD")
  result$position <- seq((windowSize/2),(((windows-1)*rollEvery)+(windowSize/2)), by = rollEvery) # Eliminated the +1 from before
  ## Generating the results
  for( i in 1:windows){
    a <- ((result$position[i]) - (windowSize/2)) + 1
    b <- (result$position[i]) + (windowSize/2)
    result$average[i] <- mean(data[a:b,2], na.rm = TRUE)
    result$SD[i] <- sd(data[a:b,2], na.rm = TRUE)
  }
  return(result)
}

# A function for eliminanting outliers from a region of the genome
GenomeQuantileCensorer <- function(data, difference, censorSize){
  if((dim(data)[1]) <= censorSize) {
      stop('The genome length is less than the area you want to censor. If this is what you want, use the function QuantileCensorer; else, check your Censor Size')
  }
  # create final data frame
  final <- data.frame()
  last <- data.frame()

  times <- (dim(data)[1])%/%censorSize
  for(i in 1:times){
    if(i == times){
      a <- ((i-1)*censorSize)+1
      b <- (dim(data)[1])
    } else {
      a <- (i*censorSize)-(censorSize -1)
      b <- (i*censorSize)
    }
  last <- QuantileCensorer(data[a:b,], difference)
  final <- rbind(final,last)
  }
  return(final)  
}

#----------- Secondary functions -----------
# This function is for doing the whole MFA! From reaidng files to getting the Loess data
doMFALoess <- function(reading_Path, genome_length, bin, s){
	get_files_list <- list.files(path = reading_Path, pattern="^depth")
	RawData <- filesformiles(paste0(reading_Path,get_files_list),genome_length)
	NormalizedData <- normalizer(RawData)
	AveragedData <- simpleAverager(NormalizedData)
	BinnedData <- beaner(genome_length,AveragedData,bin)
	LoessData <- linemaker(BinnedData,s)

	return(LoessData)  
}

#This function is for doing the whole MFA! From reading files to binning data
doMFABin <- function(reading_Path, genome_length, bin){
	get_files_list <- list.files(path = reading_Path, pattern="^depth")
	RawData <- filesformiles(paste0(reading_Path,get_files_list),genome_length)
	NormalizedData <- normalizer(RawData)
	AveragedData <- simpleAverager(NormalizedData)
	BinnedData <- beaner(genome_length,AveragedData,bin)

	return(BinnedData)  
}

# This function creates a linear model from MFA data. stratrs from reading data and gives the model as output
doMFALinearModel <- function(reading_Path, genome_length, bin, window, difference, genomeWindow) {
	get_files_list <- list.files(path = reading_Path, pattern="^Depth")
	RawData <- filesformiles(paste0(reading_Path,get_files_list),genome_length)
	colnames(RawData) <- c("postion", "average")
	NormalizedData <- normalizer(RawData)
	AveragedData <- simpleAverager(NormalizedData)
	RolledAveragedData <- RollingAverage(AveragedData, bin, window)
	SmoothData <- GenomeQuantileCensorer(RolledAveragedData, difference, genomeWindow)
	SmoothData[2] <- log(SmoothData[2],2)
	# The linear model
	linearSample1 <- SmoothData[1:1891,]
	linearSample2 <- SmoothData[1892:dim(SmoothData)[1],]
	linearMod1 <- lm(average ~ position, data=linearSample1)
	linearMod2 <- lm(average ~ position, data=linearSample2)
	# The prediction
	predictDF <- as.data.frame(RolledAveragedData$position)
	left <- as.data.frame(predictDF[1:1891,])
	right <- as.data.frame(predictDF[1892:dim(SmoothData)[1],])
	colnames(left) <- "position"
	colnames(right) <- "position"
	left$average <- predict(linearMod1,left)
	right$average <- predict(linearMod2, right)
	MFA_model <- rbind(left,right)
	MFA_model[2] <- 2^MFA_model[2]

	return(MFA_model)
}


