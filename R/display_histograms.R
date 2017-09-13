display_histogram <- function(histograms, type, obs.index){
  
  if( type=='bht' ){
    hist(x = as.numeric(histograms[obs.index,]) , breaks = 100 )
    return;
  }
  
  if( type=='dbht' ){
    
    return;
  }
  
  return;
}