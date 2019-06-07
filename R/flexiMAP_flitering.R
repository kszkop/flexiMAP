flexiMAP_filtering <- function(UTRcounts, listBams, normalise=FALSE, listBam, exprFilt=20)
{
  if(normalise==TRUE){
    #Calculate library sizes
    norm_fact <- as.numeric()
    for(lib in 1:length(listBams)){
      tmpBam <- Rsamtools::BamFile(listBams[lib])
      norm_fact[lib] <- as.numeric(Rsamtools::countBam(tmpBam)[6])/1000000
    }
  }
  #Initate output
  UTRcounts_filt <- list()

  #Extract counts for regions
  long_count <- UTRcounts[[1]]
  short_count <- UTRcounts[[2]]

  #Normalise counts for library sizes if normalise true
  if(normalise==TRUE){
    long_count_norm <- long_count/norm_fact
    short_count_norm <- short_count/norm_fact
  } else {
    long_count_norm <- long_count
    short_count_norm <- short_count
  }
  #Calculate number of samples with reads in short region > exprFilt
  nSamples <- apply(short_count_norm > exprFilt, 1, sum)
  filtOut <- names(nSamples)[nSamples<ncol(short_count_norm)]

  #Remove from short
  short_count_filt <- short_count_norm[!rownames(short_count_norm) %in% filtOut,]

  filtOut_all <- setdiff(row.names(long_count_norm),sub('\\..*','',row.names(short_count_filt)))

  #Remove from long
  long_count_filt <- long_count_norm[!rownames(long_count_norm) %in% filtOut_all,]

  #Calculate
  long_long <- data.frame(matrix(0, nrow=nrow(short_count_filt),ncol=ncol(short_count_filt)))
  row.names(long_long) <- row.names(short_count_filt)
  colnames(long_long) <- colnames(short_count_filt)

  for (s in 1:ncol(short_count_filt)){
    for (v in 1:nrow(short_count_filt)){
      id <- sub('\\..*','',row.names(short_count_filt)[v])
      long_long[v,s] <- long_count_filt[id,s] - short_count_filt[v,s]
    }
  }

  #if long_long is below 0 , replace with 0
  long_long[long_long<0] <- 0

  #Output
  UTRcounts_filt[[1]] <- long_long
  UTRcounts_filt[[2]] <- short_count_filt

  return(UTRcounts_filt)
}
