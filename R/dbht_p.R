##
## DBHT parallel  
##

##
##
dbht_p <- function( obs_index, s , covariate.data  , B  , B.N  ) {
  
  actual_data <- cbind( covariate.data , s[,1] , s[,2] )
  time_index   <- ncol(actual_data) - 1
  status_index <- ncol(actual_data)  
  
  data_length  <- nrow(actual_data)

  antidote_concs_vector  <- rep(x =  0 , B)
  poison_concs_vector    <- rep(x =  0 , B)
  
  ## calculate the baseline concordance, according to the baseline (all observations)
  cox_object <- survival::coxph( survival::Surv(actual_data[,time_index], as.integer(actual_data[,status_index]) ) ~ .   , data = data.frame(actual_data[,-c(time_index,status_index)]) )
  baseline_concordance <- cox_object$concordance[1]/(cox_object$concordance[1] + cox_object$concordance[2] )    
  
  
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
    antidote_concs_vector[i] <- actual_conc - baseline_concordance
    
   
  }
  

  for ( i in 1:B ){
    boot_datax  <- getBootstrapBiased_k(data = actual_data , obs_index = obs_index, k = B.N)
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
    poison_concs_vector[i] <- actual_conc - baseline_concordance
    
  }
  #print( paste(B," resamples fit for observation " , obs_index , sep="" ))
  
  
  ##
  ## Calculate p-value of hypothesis test : H0:= delta<=0 ,. H1: delta>0
  ##
  pvalue  <- pvalue_poison_antidote_t(hist_lower = poison_concs_vector , hist_higher = antidote_concs_vector )
  
  histograms         <- list( antidote_concs_vector , poison_concs_vector);
  names(histograms)  <- c("antidote","poison");
  dbht_output        <- list(c(obs_index,pvalue) , histograms );
  names(dbht_output) <- c("metrics","histograms");
    
  return( dbht_output )
}