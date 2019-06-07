flexiMAP_counting <- function(listBams, UTRannotation)
{
  #Initate output
  UTRcounts <- list()

  #Extract annototion for counting
  featureCounts_annot_long <- UTRannotation[[1]]
  featureCounts_annot_short <- UTRannotation[[2]]

  #Counting
  long_count <- Rsubread::featureCounts(listBams,annot.ext=featureCounts_annot_long,allowMultiOverlap = FALSE, ignoreDup=TRUE,useMetaFeatures=FALSE)
  short_count <- Rsubread::featureCounts(listBams,annot.ext=featureCounts_annot_short,allowMultiOverlap=TRUE, ignoreDup=TRUE,useMetaFeatures=FALSE)

  #Output
  UTRcounts[[1]] <- as.data.frame(long_count$counts)
  UTRcounts[[2]] <- as.data.frame(short_count$counts)

  return(UTRcounts)
}
