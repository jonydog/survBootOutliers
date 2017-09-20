#
# This file implements the three alternative methods considered in the Bioninformatics paper: 
#
#   1- Martingale Residuals
#
#

#' @importFrom survival coxph
wod_coxmartingaleresiduals <- function(surv.object , covariate.data ){
  
  
  actual_data <- cbind( covariate.data , surv.object[,1] , surv.object[,2] )
  time_index   <- ncol(actual_data) - 1
  status_index <- ncol(actual_data) 
  
  cox_object  <- survival::coxph( survival::Surv(actual_data[,time_index], as.integer(actual_data[,status_index]) ) ~ .   , data = data.frame(actual_data[,-c(time_index,status_index)]) )
  
  
  
  ress <- stats::resid(cox_object,type="martingale")
  ress <- data.frame(cbind(1:length(ress), abs(ress) ) )
  ordered_ress <- ress[ order( -ress[,2]) , ]
  
  return(ordered_ress)
  
}