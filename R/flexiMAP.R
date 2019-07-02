#' @title flexiMAP: Flexible Modelling of Alternative Polyadenylation
#'
#' @description A beta-regression-based method for discovering differential alternative polyadenylation events in standard RNA-seq data.
#'
#' @param InputTable  Dataframe containing names of bam files in the first column, condition in format 1,2,3..(and so on) in the second column. Additional columns with optional covariates for modelling with names of covariates as column names.
#' @param path  Path to folder containg indexed bam files
#' @param covariates names of covariates (additional to condition, the same as column names in InputTable)
#' @param polyAsite Dataframe containing polyadenylation sites in format: GeneID, chromosome, strand, genomic coordinate
#' @param reference Dataframe containing reference in format: GeneID, transID, chromosome, strand, transcript start, transcript end, CDS start, CDS end
#' @param TINfilter optional, a vector of gene symbols to be removed from analysis because of degradation bias. The TIN values can be obtained by running RSeQC software. The genes with low TIN values have high 3' bias and should be removed from analysis
#' @param num_samples minimal number of samples with ratio (0,1)
#' @param link default "logit", options: "logit", "probit", "cloglog", "cauchit", "log", "loglog", see ?betareg for explanation
#' @param type default "ML", options: "ML", "BC", "BR", see ?betareg for explanation
#' @param link.phi default "log", options:  "identity", "log", "sqrt",  see ?betareg for explanation
#' @param name optional name of output file if different from default
#' @param exprFilt minimal number counts in short region in all samples, default 20
#' @param normalise Default FALSE. Optionally the counts can be normalised for library sizes
#'
#' @return An object containing p-values, adjusted p-values for each tested polyadenylation proximal site of a gene. In addition the output contain ratios for each sample. THe output files is also written out to working directory.
#'
#' @examples #Format of InputTable dataframe
#' head(InputTable)
#' Sample         Condition Sex
#' sample1.bam    1         1
#' sample2.bam    1         2
#' sample3.bam    2         1
#' sample4.bam    2         2
#'
#' @examples #Format of reference dataframe
#' head(reference)
#' GeneID  TransID   Chromosome  Strand  TransStart  TransEnd  CDSstart  CDSend
#' BHLHE23 NM_080606 chr20       -       61637330    61638387  61637400  61638126
#' CRYBA2  NM_057094 chr2        -       219854911   219858127 219854973 219857898
#' FNDC8   NM_017559 chr17       +       33448630    33457751  33448712  33457453
#'
#' @examples #Format of polyAsite dataframe
#' head(polyAsite)
#' GeneID Chromosome Strand GenomicCoordinate
#' CYP2R1 chr11      +      14910363
#' LSM1   chr8       +      38021091
#' SIK3   chr11      -      116820211
#'
#' @examples #Highly degraded genes with 3'bias shoud be removed from analysis
#' TINfilter = c(CYP2R1, LSM1)
#'
#' #Run Analysis
#' flexiMAP(InputTable, path='/home/bamFolder/', covariates=c('sex'), polyAsite, reference, TINfilter = c(CYP2R1, LSM1), num_samples = 1, link ='logit', type = 'ML', link.phi='log', name = 'flexiMAP_out' , normalise=FALSE, exprFilt = 20)
#'
#'
#'
#'


flexiMAP <- function(InputTable, path=NULL, covariates=NULL, polyAsite , reference , TINfilter = NULL,
                    num_samples = 1, link ="logit", type ="ML", link.phi="log",
                    name = "flexiMAP_out" , normalise=FALSE, exprFilt = 20)
{
  if(is.null(InputTable)){
    stop("Please provide an InputTable.")
  }
  if(!is.data.frame(InputTable)){
    stop("InputTable must be dataframe")
  }
  if(is.null(polyAsite)){
    stop("Please provide performQC parameter. Must be either TRUE or FALSE.\n")
  }
  if(!is.data.frame(polyAsite)){
    stop("polyAsite must be dataframe")
  }
  if(is.null(reference)){
    stop("Please provide performROT parameter. Must be either TRUE or FALSE.\n")
  }
  if(!is.data.frame(reference)){
    stop("reference must be dataframe")
  }
  #if(!length(which(c("logit", "probit", "cloglog", "cauchit", "log", "loglog")==link))==0){
  #  stop("Please provide link parameter. Must be on of: 'logit','probit','cloglog', 'cauchit', 'log', 'loglog'.For explanation see ?betareg \n")
  #}
  #if(!length(which(c("ML", "BC","BR")==type))==0){
  #  stop("Please provide link parameter. Must be on of: 'ML', 'BC', 'BR'.For explanation see ?betareg \n")
  #}
  #if(!length(which(c("identity", "log", "sqrt")==link.phi))==0){
  #  stop("Please provide link parameter. Must be on of: 'identity', 'log', 'sqrt'.For explanation see ?betareg \n")
  #}
  if(!normalise%in%c("TRUE","FALSE")){
    stop("normalise parameter must be either TRUE or FALSE.\n")
  }
  if(!is.numeric(exprFilt)){
    stop("exprFilt must be a number")
  }
  if(!is.numeric(num_samples)){
    stop("num_samples must be a number")
  }
  if(!is.null(path)){
    if(unlist(strsplit(path,''))[length(unlist(strsplit(path,'')))]!='/'){
      path <- paste(path,'/',sep='')
    }
  }
  #Prepare anotation
  UTRannotation <- flexiMAP_annotPrep(polyAsite=polyAsite,reference=reference,TINfilter=TINfilter)
  #Counting
  listBams <- sapply(as.character(InputTable[,1]), function(x) ifelse(is.null(path), gsub(' ','', x), gsub(' ','',paste(path,x,sep=''))),simplify = TRUE, USE.NAMES = F)

  UTRcounts <- flexiMAP_counting(listBams = listBams, UTRannotation = UTRannotation)
  #Filtering
  UTRcounts_filt <- flexiMAP_filtering(UTRcounts = UTRcounts, listBams = listBams, normalise = normalise, exprFilt = exprFilt)
  #beta-regression
  results <- flexiMAP_stat(UTRcounts_filt = UTRcounts_filt, InputTable = InputTable, num_samples = num_samples, covariates = covariates,link = link , type = type, link.phi = link.phi,name=name)

  return(results)
}
