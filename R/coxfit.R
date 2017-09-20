##
## File that contains the Cox fitting function used across all implemented  methods
##

#' @importFrom survival coxph
coxph.call <- function (boot_data, time_index, status_index) {
  actual_conc  <- tryCatch({
    cox_object <- coxph( Surv(boot_data[,time_index], as.integer(boot_data[,status_index]) ) ~ .,
                         data = data.frame(boot_data[,-c(time_index,status_index)]),
                         control = survival::coxph.control(iter.max = 1000, outer.max = 1000), ties = 'efron')
    cox_object$concordance[1]/(cox_object$concordance[1] + cox_object$concordance[2] )
  }, error = function(e) {
    cat('Error! trying to resample 10 times')
    return(NULL)
  })
  return(actual_conc)
}