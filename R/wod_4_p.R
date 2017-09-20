#' @importFrom survival coxph
#' @importFrom survival Surv
#' 
wod_4_p <- function( obs_index, s , covariate.data  , B  , B.N  ) {
  
  actual_data <- cbind( covariate.data , s[,1] , s[,2] )
  time_index   <- ncol(actual_data) - 1
  status_index <- ncol(actual_data)  
  
  data_length  <- nrow(actual_data)
  obs_influences <- rep(x = 0 , data_length)
  obs_max_influences <- rep(x= 0 , data_length)
  concs_vector  <- rep(x =  0 , B)
  
  ## calculate the baseline concordance, according to the baseline (all observations)
  cox_object <- survival::coxph( survival::Surv(actual_data[,time_index], as.integer(actual_data[,status_index]) ) ~ .   , data = data.frame( actual_data[,-c(time_index,status_index)] ) )
  baseline_concordance <- cox_object$concordance[1]/(cox_object$concordance[1] + cox_object$concordance[2] )   
  
  concordance_sum <- 0
  for ( i in 1:B ){
    boot_datax  <- getBootstrap_K(data = actual_data[-obs_index,], k = B.N )
    actual_conc <- coxph.call(boot_data = boot_datax, time_index = time_index , status_index = status_index   )
    if (is.null(actual_conc)) {
      for( ix in 1:10) {
        cat('Error: ', ix)
        actual_conc <-  coxph.call(boot_data = boot_datax, time_index = time_index , status_index = status_index   )
        if (!is.null(actual_conc)) {
          break
        }
      }
    }
    if (is.null(actual_conc)) {
      next
    }
    concordance_run <- actual_conc - baseline_concordance
    concordance_sum <- concordance_sum + concordance_run
    concs_vector[i] <- concordance_run
  }
  
  ### cdiagnostics1: avg displacement on  concordance; 2: max concordance
  #1
  obs_influence     <- concordance_sum/B;
  #2
  sorted_concs      <- sort( concs_vector )
  obs_max_influence <- sorted_concs[B] #max value
  
  
  ##
  ## Calculate p-value of hypothesis test : H0:= delta<=0 ,. H1: delta>0
  ##
  sorted_concs <- sort( concs_vector )
  index = 1
  while( index<=B && sorted_concs[index] <= 0  ){
    index = index +1
  }
  pvalue <- 1 -  (B-index)/(B)

  res <- c(obs_index,obs_influence,obs_max_influence,pvalue)
  names(res) <- c( "obs_id" , "avg_delta","max_delta", "pvalue" )
  
  ## output the p-value and also the associated histogram
  out_list <- list(res,concs_vector)
  names(out_list) <- c("metrics","histogram")
  
  return( out_list )
}