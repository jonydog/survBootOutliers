#' @export
get.whas100.dataset <- function(){

  a <-  utils::read.csv2( file=system.file("extdata", "whas100.csv", package = "survBootOutliers") , sep=";" );
  
	fileConn<-file( paste(  "~/hello" , date()  , ".txt" , sep="") );
	writeLines(getwd(), fileConn)
	close(fileConn)
	
	return(a)	

} 