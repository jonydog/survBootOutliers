##
##
##
## Bootstrap versions for the BCSOD parallel version
##
##
##

##
## Auxiliar function, that creates bootstrap set from the total set
##
## it receives the whole data 

getBootstrapBiased_k <- function( data = "data.frame", obs_index , k){
  
  
  if( ncol(data)> 1 ){
    
    bootstrap_data <- data[1:k,]   #alocar bootstrap data, data mantem-se igual, bootdata tem menosum elemento
    
    boot_indexes   <- sample(1:nrow(data),(k-1),replace=TRUE) 
    
    bootstrap_data[1:(k-1),]  <- data[boot_indexes,]
    
    bootstrap_data[k,]        <- data[obs_index,] 
    
    return(bootstrap_data)  
  }
  else{ ##  for the case when the covariate data  is not  a matrix (just a vector)
     
    bootstrap_data <- data[1:k,]   #alocar bootstrap data, data mantem-se igual, bootdata tem menosum elemento
    
    boot_indexes   <- sample(1:nrow(data),(k-1),replace=TRUE) 
    
    bootstrap_data[1:(k-1)]  <- data[boot_indexes,]
    
    bootstrap_data[k]        <- data[obs_index,]  
    
    return(bootstrap_data)  
  }
  
  return(bootstrap_data)  
}  



##  
## k out of n bootstrap function
## 
getBootstrap_K <- function( data = "data.frame" , k ){
  
  if( ncol(data)>1 ){
    boot_indexes   <- sample(1:nrow(data),k,replace=TRUE)  
    bootstrap_data <- data[boot_indexes,]
  }
  else{
    boot_indexes   <- sample(1:nrow(data),k,replace=TRUE)  
    bootstrap_data <- data[boot_indexes,]
  }
  return(bootstrap_data)
} 
