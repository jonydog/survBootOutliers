##
## wod_4hist: Observation influence estimation by bootstrap , returns the matrix with all measures
##
##
## Function wod_4hist
##
##                  performance metric : concordance index

##  Arguments:      data  : SurvivalDataSet object
##                  nruns  : number of runs of the bootstrap
##                  layer5_model  : the layer5 predictor
##  Output   :      vector with length equal to the number of data rows , each position with the estimated influence on concordance
##

#' @importFrom survival coxph
wod_4hist_k_p <-  function( obs_index, surv.object, covariate.data , nruns  , m ) {
  
  
  
  actual_data <- cbind( covariate.data , surv.object[,1],  surv.object[,2] )
  
  time_index   <- ncol(actual_data) - 1
  status_index <- ncol(actual_data)  
  
  
  data_length  <- nrow(actual_data)
  
  obs_influences <- rep(x = 0 , data_length)
  
  obs_max_influences <- rep(x= 0 , data_length)
  
  
  concs_vector  <- rep(x = 0, nruns  )
  
  ## calculate the baseline concordance, according to the baseline (all observations)
  #
  cox_object <- coxph( survival::Surv(actual_data[,time_index], as.integer(actual_data[,status_index]) ) ~ .   , data = data.frame(actual_data[,-c(time_index,status_index)]) )
  baseline_concordance <- cox_object$concordance[1]/(cox_object$concordance[1] + cox_object$concordance[2] )   
  
  cat("baseline:" , baseline_concordance )
  
    
    concordance_sum <- 0
    
    
    for ( i in 1:nruns ){
      boot_datax  <- getBootstrap_K(data = actual_data[-obs_index,], )
      actual_conc <- coxph.call(boot_data = boot_datax , time_index = time_index, status_index = status_index )
      if (is.null(actual_conc)) {
        for( ix in 1:10) {
          cat('Error: ', ix)
          boot_datax  <- getBootstrap_K(data = actual_data[-obs_index,], )
          actual_conc <- coxph.call(boot_data = boot_datax , time_index = time_index, status_index = status_index )
          if (!is.null(actual_conc)) {
            break
          }
        }
      }
      if (is.null(actual_conc)) {
        next
      }
      
      concordance_run <- actual_conc - baseline_concordance
      
      concs_vector[i] <- concordance_run     
    }
    
    
  #  cat(  "Antidote Bootstrap: " , row_index , "out of", data_length , " complete \n")
  
  return(concs_vector)
  
}





