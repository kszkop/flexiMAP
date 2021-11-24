flexiMAP_counting <- function(listBams, UTRannotation, nthreads, pairedEnd=FALSE)
{
  #Initate output
  UTRcounts <- list()

  #Extract annototion for counting
  featureCounts_annot_long <- UTRannotation[[1]]
  featureCounts_annot_short <- UTRannotation[[2]]

  #Counting
  if(!isTRUE(pairedEnd){
    long_count <- Rsubread::featureCounts(listBams,annot.ext=featureCounts_annot_long,allowMultiOverlap = FALSE, ignoreDup=TRUE,useMetaFeatures=FALSE,nthreads=nthreads, isPairedEnd=FALSE)
    short_count <- Rsubread::featureCounts(listBams,annot.ext=featureCounts_annot_short,allowMultiOverlap=TRUE, ignoreDup=TRUE,useMetaFeatures=FALSE,nthreads=nthreadsis, isPairedEnd=FALSE)
  } else {
    long_count <- Rsubread::featureCounts(listBams,annot.ext=featureCounts_annot_long,allowMultiOverlap = FALSE, ignoreDup=TRUE,useMetaFeatures=FALSE,nthreads=nthreads, isPairedEnd=TRUE, countReadPairs = TRUE, requireBothEndsMapped = FALSE)
    short_count <- Rsubread::featureCounts(listBams,annot.ext=featureCounts_annot_short,allowMultiOverlap=TRUE, ignoreDup=TRUE,useMetaFeatures=FALSE,nthreads=nthreads,isPairedEnd=TRUE, countReadPairs = TRUE, requireBothEndsMapped = FALSE)
  }
  #Output
  UTRcounts[[1]] <- as.data.frame(long_count$counts)
  UTRcounts[[2]] <- as.data.frame(short_count$counts)

  return(UTRcounts)
}
