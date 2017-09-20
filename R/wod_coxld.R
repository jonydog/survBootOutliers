##
##
##
## Likelihood displacement statistic
##
##
##
##


## Function wod_ld 
##
##              


##  Output   :      vector with length equal to the number of data rows , each position with the ld statistic
##


#' @importFrom survival coxph
wod_coxld <-  function( surv.object , covariate.data ) {
  
  actual_data <- cbind( covariate.data , surv.object[,1] , surv.object[,2] )
  time_index   <- ncol(actual_data) - 1
  status_index <- ncol(actual_data) 
  
  cox_object  <- coxph( survival::Surv(actual_data[,time_index], as.integer(actual_data[,status_index]) ) ~ .   , data = data.frame(actual_data[,-c(time_index,status_index)]) )
  loglik_baseline <- cox_object$loglik[2]
  
  
  lds_vector <- rep(0,nrow(actual_data))

  for(i in 1:nrow(actual_data)){
    
    cox_fit_i  <- coxph( survival::Surv(actual_data[-i,time_index], as.integer(actual_data[-i,status_index]) ) ~ .   , data = data.frame(actual_data[-i,-c(time_index,status_index)]) )
    
    lds_vector[i] <- 2*(cox_fit_i$loglik[2]-loglik_baseline)
    
  }
  
  ## complete the output with the row.names of the data
  #
  
  out <- cbind(1:nrow(actual_data), lds_vector , as.numeric(row.names(actual_data)) )
  
  out_ordered <- out[order(-out[,2] ), ]
  
  row.names(out_ordered) <- as.character(out_ordered[,3])
  
  out_ordered <- out_ordered[,-3]  
  
    
  return(out_ordered)

}