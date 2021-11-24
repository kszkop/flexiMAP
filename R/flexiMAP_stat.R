flexiMAP_stat <- function(UTRcounts_filt, InputTable, num_samples, covariates=NULL,
                                  link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                                  type = c("ML", "BC", "BR"), link.phi=c("identity", "log", "sqrt"),name){
  #extract filtered counts
  long_long <- UTRcounts_filt[[1]]
  short <- UTRcounts_filt[[2]]

  #Create design formula with optional covariates
  if(!is.null(covariates)){
    design <- paste('ratio ~ condition + ',paste(covariates,collapse=' +'))
  } else {
    design <- paste('ratio ~ condition | condition')
  }

  #Initate output
  res  <- matrix(NA,nrow=nrow(short),ncol=ifelse(!is.null(covariates),3+length(covariates),3))
  if(!is.null(covariates)){
    colnames(res) <- c('ID','intercept','condition',covariates)
  } else {
    colnames(res) <- c('ID','intercept','condition')
  }
  res <- as.data.frame(res)
  res$ID <-  row.names(short)

  #save ratios
  ratioOut <- list()

  for (i in 1:nrow(short)){
    #Progress
    pb1 <- txtProgressBar(min=1, max=nrow(short), style=3)
    condition <- as.factor(InputTable[,2])
    print(i)
    print(condition)
    ratioOut[[i]] <- as.numeric(long_long[i,])/(as.numeric(long_long[i,])+as.numeric(short[i,]))
    ratio <- ratioOut[[i]]
    ratio_len1 <- length(which(ratio[condition==1]>0 & ratio[condition==1]<1))
    ratio_len2 <- length(which(ratio[condition==2]>0 & ratio[condition==2]<1))
    if(ratio_len1 >= num_samples & ratio_len2 >= num_samples){
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
        if(!is.null(covariates)){
          for(var in 1:length(covariates)){
            assign(covariates[var],as.factor(InputTable[,covariates[var]]))
          }
        }
      }
      br.out <- tryCatch(betareg::betareg(as.formula(design), link = link, type = type, link.phi=link.phi),error=function(e) NA, warning=function(w) NA)
      if(length(which(is.na(br.out)))==0){
        summary_br <- summary(br.out)
        coef_br <- summary_br$coefficients[[1]][,4]
        res[i,2:ifelse(!is.null(covariates),3+length(covariates),3)] <- coef_br
      } else {
        res[i,2:ifelse(!is.null(covariates),3+length(covariates),3)] <- NA
      }
      print('breg')
    } else if(ratio_len1 >= num_samples | ratio_len2 >= num_samples){
      glm.out <- tryCatch(glm(cbind(as.numeric(long_long[i,]),as.numeric(short[i,])) ~  condition, family = quasibinomial, maxit=1000),error=function(e) NA, warning=function(w) NA)
      #if(length(which(is.na(glm.out)))==0){
      summary_glm <- summary(glm.out)
      coef_glm <- summary_glm$coefficients[,4]
      res[i,2:3] <- coef_glm
      #} else {
      #  res[i,2:3] <- NA
      #}
      print('glm')
    }
    setTxtProgressBar(pb1, i)
  }
  #adjust for multiple testing
  res$adj_cond <- p.adjust(res[,3])

  #Prepare and write our output
  ratio_df <- data.frame(matrix(unlist(ratioOut), nrow=length(ratioOut), byrow=T),stringsAsFactors=FALSE)
  res_all <- cbind(res,ratio_df)
  write.table(res_all, file=paste(name, 'txt', sep='.'), sep='\t',quote=F,row.names=F)

  return(res_all)
}
