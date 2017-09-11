#' This function retrieves the well known Worcester Heart Attack dataset with 100 individuals (WHAS100). This dataset is taken from the book by Hosmer D, Lemeshow S, May S. Applied Survival Analysis: Regression Modeling of Time to Event Data, 2nd edition. John Wiley and Sons Inc., New York, NY, 2008.
#'
#' @param none
#' 
#' @return A data.frame containing the most outlying observations sorted by outlying score
#' 
#' @examples ## One Step Deletion "osd" method
#' whas <- get.whas100.dataset()
#' 
#' @export
get.whas100.dataset <- function(){

  a <-  utils::read.csv2( file=system.file("extdata", "whas100.csv", package = "survBootOutliers") , sep=";" );
  
	return(a)	

} 