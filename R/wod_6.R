##
##
## WOD_6 : Greedy maximization of concordance by greedy one step ahead , remove one observation at a time, 
##                                                                       removing the one that improves concordancde the most (estimated without bootstrap)

## Function wod_4
##
##                  performance metric : concordance index

##  Arguments:      data  : SurvivalDataSet object
##                  
##                  layer5_model  : the layer5 predictor


##  Output   :      vector with length equal to the number of data rows , each position with the estimated influence on concordance
##


#' @importFrom survival coxph
wod_6 <-  function( surv.object, covariate.data ,  max_outliers ) {
  
  actual_data <- cbind( covariate.data , surv.object[,1] , surv.object[,2] )
  
  time_index   <- ncol(actual_data) - 1
  status_index <- ncol(actual_data)  
  
  
  data_length  <- nrow(actual_data)
  
  
  ## calculate the baseline concordance, according to the baseline (all observations)
  #
  cox_object <- coxph( survival::Surv(actual_data[,time_index], as.integer(actual_data[,status_index]) ) ~ .   , data = data.frame(actual_data[,-c(time_index,status_index)])  )
  baseline_concordance <- cox_object$concordance[1]/(cox_object$concordance[1] + cox_object$concordance[2] )   

  ## initialization of index vectors
  #
  left_indexes      <- c(1:data_length)   # starts full
  removed_indexes   <- c()                # starts empty
  
    
    
  for( outlier_index in 1:max_outliers ){
    
    actual_concordance_delta  <- 0  # by default
  
    for( i in 1:length(left_indexes) ){
      
        ## remove one more observation beyond the already removed ones
        #
        cox_object <- coxph( Surv(actual_data[-c(removed_indexes,left_indexes[i] ),time_index], as.integer(actual_data[-c(removed_indexes,left_indexes[i]),status_index]) ) ~ .   , data = data.frame(actual_data[-c(removed_indexes,left_indexes[i]),-c(time_index,status_index)] ) )
        
        concordance_run_delta <- cox_object$concordance[1]/(cox_object$concordance[1] + cox_object$concordance[2] ) - baseline_concordance
      
        if( concordance_run_delta > actual_concordance_delta ){
          
          actual_index <- i
          actual_concordance_delta <- concordance_run_delta
        }
        
    }
    
    if( actual_concordance_delta == 0){
      return(removed_indexes)
    }
    else{
      
      removed_indexes <- c(removed_indexes,left_indexes[actual_index] )
      
      
      left_indexes    <- left_indexes[-actual_index]
    }
  
  }
  
  removed_indexes <- data.frame(removed_indexes)
  
  output <- removed_indexes

  return( output )
  
}






