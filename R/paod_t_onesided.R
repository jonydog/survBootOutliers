##
## Receives histograms from poison and antidote bootstraps
##
##
##
#' @importFrom stats t.test
pvalue_poison_antidote_t <- function(hist_lower, hist_higher) {
  hist1 <- hist_lower
  hist2 <- hist_higher
  
  #T <- (mean(hist1)-mean(hist2) )/ sqrt( var(hist1)/length(hist1) + var(hist2)/length(hist2)  )
  #
  #print(hist1)
  #print(hist2)
  #
  test_done <- t.test(x = hist2 ,
              y = hist1,
              alternative = "greater")
  
  ## for large sample T follows standard normal distribution
  return(test_done$p.value)
}


## It receives a vector of values (a matrix) with the runs for poison and antidote bootstraps
##
## Arguments : Matrix where for each variable we have a column, each row represents a bootstrap run
##
paod_t_onesided <- function(poison_boot, antidote_boot) {
  if (ncol(poison_boot) != ncol(antidote_boot) ||
      nrow(poison_boot) != nrow(antidote_boot)) {
    stop("Histograms sizes do not match!")
  }
  
  nvars <- ncol(poison_boot)
  nruns <- nrow(poison_boot)
  
  pvalues <- rep(x = NA, nvars)
  
  
  for (i in 1:nvars) {
    pvalues[i]  <-
      pvalue_poison_antidote_t(hist_lower = poison_boot[, i] , hist_higher = antidote_boot[, i])
  }
  
  pre_out <- cbind(1:nvars , pvalues)
  
  out <- pre_out[order(pre_out[, 2]) , ]
  
  return(out)
}
