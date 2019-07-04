flexiMAP_annotPrep <- function(polyAsite, reference, TINfilter)
{
  #Initate output
  UTRannotation <- list()

  colnames(reference)[1:8] <- c('geneID','transID','chr','strand','Tstart','Tend','CDSstart','CDSend')
  colnames(polyAsite)[1:4] <- c('geneID','chr_polyAsite','strand_polyAsite','polyAsite')

  #Extract 3'UTR start site
  reference$UTR_start <- ifelse(reference$strand=='+', reference$CDSend, reference$CDSstart)

  #Remove transcripts of the same gene that have the same UTR start site by taking random of the transcripts as represantative and take only these transcript that have the UTR start site that is equal or further downstream of last part of coding region across isoforms
  reference <- reference[!duplicated(reference[c("geneID", "UTR_start")]),]
  reference <- data.table::setDT(reference)[, .(transID = transID, chr=chr, strand=strand, Tstart=Tstart, Tend=Tend, CDSstart=CDSstart, UTR_start=UTR_start, test = ifelse(strand=='+', max(CDSend) <= UTR_start, min(CDSstart) >= UTR_start)), by = 'geneID']
  reference <- reference[reference$test==TRUE,]

  #Merge transcript annotation with polyAsites annotation
  annot_polyA <- merge(polyAsite[,1:4], reference, by='geneID')

  #Select only sites that are sense to gene
  annot_polyA <- subset(annot_polyA, annot_polyA$strand_polyAsite==annot_polyA$strand)

  #Select only sites that strart downstream from end of CDS descripbed above.
  annot_polyA <- annot_polyA[ifelse(annot_polyA$strand=='+',annot_polyA$polyAsite > annot_polyA$UTR_start, annot_polyA$polyAsite < annot_polyA$UTR_start),]

  #Select genes with at least 2 sites in 3'UTR
  annot_polyA_counts <- table(annot_polyA$geneID)
  annot_polyA <- annot_polyA[annot_polyA$geneID %in% names(annot_polyA_counts[annot_polyA_counts>1]),]

  #If privided filters out transcripts not in the character vector
  if(!is.null(TINfilter)){
    annot_polyA <- annot_polyA[annot_polyA$transID %in% TINfilter,]
  }

  #order the sites
  annot_polyA <- annot_polyA[with(annot_polyA,order(annot_polyA$geneID,ifelse(annot_polyA$strand=='+', -annot_polyA$polyAsite, annot_polyA$polyAsite))),]

  #Select the distal site - which is defined as the furthest downstream
  distal <- subset(annot_polyA, !duplicated(geneID))
  #Select the proximal sites - take all other sites that are in 3UTR
  proximal <- annot_polyA[-c(which(!duplicated(annot_polyA$geneID))),]

  #Prepare output format
  distal <- data.frame(chr = distal$chr, start=ifelse(distal$strand=='+',distal$UTR_start,distal$polyAsite), end=ifelse(distal$strand=='+',distal$polyAsite,distal$UTR_start), geneID = distal$geneID, strand=distal$strand, stringsAsFactors = FALSE)
  proximal <- data.frame(chr = proximal$chr, start=ifelse(proximal$strand=='+',proximal$UTR_start,proximal$polyAsite), end=ifelse(proximal$strand=='+',proximal$polyAsite,proximal$UTR_start), geneID = proximal$geneID, strand=proximal$strand, stringsAsFactors = FALSE)

  #Remove genes with overlapping with 3'UTR
  ranges=GenomicRanges::GRanges(seqnames=distal$chr, IRanges::IRanges(start=distal$start,end=distal$end),names=distal$geneID)
  OverTmp <- GenomicRanges::findOverlaps(ranges)
  overlaps <- IRanges::from(OverTmp)[duplicated(IRanges::from(OverTmp))]
  overlapsGene <- as.character(ranges[overlaps]$names)

  #remove overlapping genes
  distal <- distal[!distal$geneID %in% overlapsGene,]
  proximal <- proximal[!proximal$geneID %in% overlapsGene,]

  #Merge distal and proximal
  comb <- merge(proximal,distal,by='geneID')

  #Remove UTRs with too short distance between UTRstart and Proximal (below 100 nb)
  comb <- comb[which(comb$end.x-comb$start.x>=100),]

  #Remove UTRs with too short distance between Proximal and Distal(below 100 nb)
  len <- ifelse(comb$strand.x=='+', comb$end.y-comb$end.x,comb$start.x-comb$start.y)
  comb <- comb[which(len>100),]

  #Extract final output
  long <- unique(comb[,c(1,6:9)])
  colnames(long) <- gsub('\\.y','',colnames(long))
  short <- comb[,c(1:5)]
  colnames(short) <- gsub('\\.x','',colnames(short))

  #add variant ID, order first
  short <- short[with(short, order(short$geneID, ifelse(short$strand=='+',short$end, short$start))),]
  short <- data.table::setDT(short)[, .(chr=chr,start=start,end=end, strand=strand, IDs=paste(as.character(geneID), seq(1:.N),sep='.')), by = geneID]
  short <- short[,c(6,2:5)]

  #Output
  UTRannotation[[1]] <- long
  UTRannotation[[2]] <- short

  return(UTRannotation)
}
