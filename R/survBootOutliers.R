# source("R/getwhas100dataset.R")
# source('R/bootstraps_p.R')
# source("R/coxfit.R")
# source("R/wod_4_p.R")
# source("R/wod_6.R")
# source("R/wod_coxdevianceresiduals.R")
# source("R/wod_coxld.R")
# source("R/wod_coxmartingaleresiduals.R")
# source("R/dbht_p.R")
# source("R/paod_t_onesided.R")
#' Extract the most outlying observations following a criteria based on the bootstrapped concordance with parallel processing
#' @param surv.object An obect of type survival::Surv containing lifetimes and right-censoring status
#' @param covariate.data A data frame containing the data with covariate values for each individual
#' @param sod.method One of c("osd","bht","dbht","ld","martingale","deviance")
#' @param B The number of bootstrap samples generated only applicable for "bht" and "dbht" methods
#' @param B.N the number of observations in each bootstrap sample
#' @param max.outliers This parameter is only used for the "osd" method   
#' @param mc.cores The number of cores to execute the methods "bht" and "dbht"
#' 
#' @return A data.frame containing the most outlying observations sorted by outlying score
#' 
#' @examples ## One Step Deletion "osd" method
#' whas <- get.whas100.dataset()
#' print(getwd())
#' outliers_osd <- survBootOutliers( surv.object=Surv(time = whas$times,event = whas$status ) , covariate.data = whas[,2:5] , sod.method = "osd" , max.outliers = 10 )
#' 
#' @examples ## Bootstrap Hypothesis Test "bht" with 100 bootstrap samples, each with 100 individuals and running on 4 cores
#' whas <- get.whas100.dataset()
#' outliers_bht <- survBootOutliers( surv.object=Surv(time = whas$times,event = whas$status ) , covariate.data = whas[,2:5] , sod.method = "bht" , B = 100 , B.N = 100 , parallel.param = MulticoreParam() )
#' 
#' @examples ## Dual Bootstrap Hypothesis Test "dbht" with 100 bootstrap samples, each with 100 individuals and running on 4 cores
#' whas <- get.whas100.dataset()
#' outliers_dbht <- survBootOutliers( surv.object=Surv(time = whas$times,event = whas$status ) , covariate.data = whas[,2:5] , sod.method = "dbht" , B = 100 , B.N = 100 , parallel.param = SnowParam() )
#' 
#' @examples ## One Step Deletion "osd" with an amount of 10 for maximum outlier count
#' whas <- get.whas100.dataset()
#' outliers_osd <- survBootOutliers( surv.object=Surv(time = whas$times,event = whas$status ) , covariate.data = whas[,2:5] , sod.method = "osd" , max.outliers = 10)
#' 
#' @examples ## Likelihood displacement criterion for outlier ranking
#' whas <- get.whas100.dataset()
#' outliers_ld <- survBootOutliers( surv.object=Surv(time = whas$times,event = whas$status ) , covariate.data = whas[,2:5] , sod.method = "ld")
#' 
#' @examples ## Cox regression deviance residuals criterion for outlier ranking
#' whas <- get.whas100.dataset()
#' outliers_deviance <- survBootOutliers( surv.object=Surv(time = whas$times,event = whas$status ) , covariate.data = whas[,2:5] , sod.method = "deviance")
#' 
#' @examples ## Cox regression Martingale residuals criterion for outlier ranking
#' whas <- get.whas100.dataset()
#' outliers_martingale <- survBootOutliers( surv.object=Surv(time = whas$times,event = whas$status ) , covariate.data = whas[,2:5] , sod.method = "martingale")
#'
#' @import survival
#' @import stats
#' @export
survBootOutliers <- function(surv.object, covariate.data, sod.method, B, B.N, max.outliers, parallel.param=BiocParallel::SerialParam() ){
  
  ## load library dependencies
  library(stats)
  library(survival)
  
  if( requireNamespace("BiocParallel", quietly = TRUE) ){
    library(BiocParallel);
    HAS_BIOCPARALLEL = TRUE;
  } else {
    HAS_BIOCPARALLEL = FALSE;
  }
  
  # if( mc.cores > 1){
  #   if( Sys.info()['sysname']=='Windows' ){
  #     library(parallelsugar)
  #   }else{
  #     library(parallel)
  #   }
  # }

  if( ! survival::is.Surv(surv.object) ){
    
    stop("Parameter \"surv.object\" must be an object of type \"survival::Surv\" ")
}
  
  if( ! base::is.data.frame(covariate.data) ){
    
    stop("Parameter \"covariate.data\" must be of type \"data.frame\".")
  } 
  
  if( ! is(parallel.param,"BiocParallelParam") ){
    stop("Parameter \"parallel.param\" must be of type \"BiocParallelParam\".")
  }
  
  # obtain the number of records (individuals)
  N <- attr(surv.object,"dim")[1]
  t <- surv.object[1:N,1]
  s <- surv.object[1:N,2]

  if( sod.method=="osd" ){
    outlier_set <- wod_6( surv.object = surv.object , covariate.data = covariate.data,  max_outliers = max.outliers )
    return( outlier_set)
  }
  
  if( sod.method=="bht" ){  
      
    if(HAS_BIOCPARALLEL){
      outlier_set <- BiocParallel::bplapply(X = 1:N,FUN = wod_4_p, s=surv.object, covariate.data=covariate.data, B=B , B.N=B.N, BPPARAM = parallel.param )
    }  else {
      outlier_set <- BiocParallel::bplapply(X = 1:N,FUN = wod_4_p, s=surv.object, covariate.data=covariate.data, B=B , B.N=B.N, BPPARAM = BiocParallel::SerialParam() )
    }
    ## unlist to matrix and sort by pvalue
    outlier_set_flat <-  matrix(unlist(outlier_set), ncol = 4, byrow = TRUE) 
    colnames(outlier_set_flat) <- c("obs_id","avg_delta","max_delta","pvalue")
    set_sorted  <- outlier_set_flat[order(outlier_set_flat[,4]),] 
    
    return(set_sorted)
  }

  if(sod.method=="dbht" ){    
    
    outlier_set <- BiocParallel::bplapply(X = 1:N, FUN=dbht_p, s=surv.object ,covariate.data = covariate.data , B=B, B.N=B.N, BPPARAM = parallel.param )
    outlier_set_flat <-  matrix(unlist(outlier_set), ncol = 2, byrow = TRUE) 
    colnames(outlier_set_flat) <- c("obs_id","pvalue")
    set_sorted  <- outlier_set_flat[order(outlier_set_flat[,2]),]
    return(set_sorted)
  }

  if( sod.method=="ld" ){
    outlier_set <- wod_coxld(surv.object = surv.object , covariate.data = covariate.data)
    return(outlier_set)
    
  }
  
  if( sod.method=="martingale"){
    outlier_set <- wod_coxmartingaleresiduals(surv.object = surv.object , covariate.data = covariate.data)
    return(outlier_set)
  }
  
  if( sod.method=="deviance"){
    outlier_set <- wod_coxdevianceresiduals(surv.object = surv.object , covariate.data = covariate.data)
    return(outlier_set)
  }
  
  stop("Parameter \"sod.method\" must be one of c(\"osd\",\"bht\",\"dbht\",\"ld\",\"martingale\",\"deviance\") ")
  
}