display_histogram <- function(histograms, type, obs.index){
  
  if( type=='bht' ){
    graphics::hist(x = as.numeric(histograms[[obs.index]]) , breaks = 100 )
  }
  
  else if( type=='dbht' ){
    
    poison_histogram   <- histograms[[obs.index]]$poison
    antidote_histogram <- histograms[[obs.index]]$antidote
    
    hist(poison_histogram, col=rgb(1,0,0,0.5), breaks = 100 , main='Antidote(blue) and Poison(red)', xlab='c-index variation');
    hist(antidote_histogram, col=rgb(0,0,1,0.5),  breaks=100, add=T);
    box()
  }
  
}