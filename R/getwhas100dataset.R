#' @export
get.whas100.dataset <- function(){


	a <- read.csv2(file="../inst/extdata/whas100.csv", sep=";")
	
	return(a)	

} 