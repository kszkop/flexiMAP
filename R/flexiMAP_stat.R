flexiMAP_stat <- function(UTRcounts_filt, InputTable, num_samples, covariates=NULL,
                         link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                         type = c("ML", "BC", "BR"), link.phi=c("identity", "log", "sqrt"),name)
{
  #extract filtered counts
  long_long <- UTRcounts_filt[[1]]
  short <- UTRcounts_filt[[2]]

  #Create design formula with optional covariates
  if(!is.null(covariates)){
    design <- paste('ratio ~ condition + ',paste(covariates,collapse=' +'))
  } else {
    design <- paste('ratio ~ condition ')
  }

  #Initate output
  res_br  <- matrix(NA,nrow=nrow(short),ncol=ifelse(!is.null(covariates),3+length(covariates),3))
  if(!is.null(covariates)){
    colnames(res_br) <- c('ID','intercept','condition',covariates)
  } else {
    colnames(res_br) <- c('ID','intercept','condition')
  }
  res_br <- as.data.frame(res_br)
  res_br$ID <-  row.names(short)

  #save ratios
  ratioOut <- list()

  for (i in 1:nrow(short)){
    #Progress
    pb1 <- txtProgressBar(min=1, max=nrow(short), style=3)

    ratioOut[[i]] <- as.numeric(long_long[i,])/(as.numeric(long_long[i,])+as.numeric(short[i,]))
    ratio <- ratioOut[[i]]
    if(length(ratio[ratio==0]) < num_samples){
      rm <- which(ratio <= 0 | ratio >= 1)
      if(length(rm)>0){
        ratio <- ratio[-rm]
        condition <- as.factor(InputTable[,2][-rm])
        if(!is.null(covariates)){
          for(i in 1:length(covariates)){
            assign(covariates[i],as.factor(InputTable[,covariates[i]])[-rm])
          }
        }
      } else {
        condition <- as.factor(InputTable[,2])
        if(!is.null(covariates)){
          for(var in 1:length(covariates)){
            assign(covariates[var],as.factor(InputTable[,covariates[var]]))
          }
        }
      }
      br.out <- betareg::betareg(as.formula(design), link = link, type = type, link.phi=link.phi)
      summary_br <- summary(br.out)
      coef_br <- summary_br$coefficients[[1]][,4]
      res_br[i,2:ifelse(!is.null(covariates),3+length(covariates),3)] <- coef_br
    }
    setTxtProgressBar(pb1, i)
  }
  #adjust for multiple testing
  res_br$adj_cond <- p.adjust(res_br[,3])

  #Prepare and write our output
  ratio_df <- data.frame(matrix(unlist(ratioOut), nrow=length(ratioOut), byrow=T),stringsAsFactors=FALSE)
  res_br_all <- cbind(res_br,ratio_df)
  write.table(res_br_all, file=paste(name, 'txt', sep='.'), sep='\t',quote=F,row.names=F)

  return(res_br_all)
}
