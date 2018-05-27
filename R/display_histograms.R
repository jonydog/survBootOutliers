#' Auxiliar function that displays the concordance histogram associated with the observation.
#'
#' @param histograms The histograms object returned by the survBootOutliers function when the method selected is "bht" or "dbht".
#' @param type The type of histogram that is given as input, possible choices are again "bht" or "dbht".
#' @param obs.index The original index of the observation of the concordance histograms to be displayed
#' 
#' @return No value is returned
#' 
#' @examples
#' \donttest{
#' whas <- get.whas100.dataset()
#' outliers_bht <- survBootOutliers( 
#'       surv.object=Surv(time = whas$times,event = whas$status ), 
#'      covariate.data = whas[,2:5], 
#'      sod.method = "bht", 
#'      B = 2000, B.N = 100 , 
#'      parallel.param = BiocParallel::MulticoreParam() 
#' )
#' display.obs.histogram(outliers_bht$histograms, "bht", 67)
#' }
#' 
#' @export
display.obs.histogram <- function(histograms, type, obs.index){
  
  if( type=='bht' ){
    graphics::hist(x = as.numeric(histograms[[obs.index]]) , breaks = 1000 )
  }
  
  else if( type=='dbht' ){
    
    poison_histogram   <- histograms[[obs.index]]$poison
    antidote_histogram <- histograms[[obs.index]]$antidote
    
    graphics::hist(poison_histogram, col=grDevices::rgb(1,0,0,0.5), breaks = 100 , main='Antidote(blue) and Poison(red)', xlab='c-index variation');
    graphics::hist(antidote_histogram, col=grDevices::rgb(0,0,1,0.5),  breaks=100, add=TRUE);
    graphics::box()
  }
  
}