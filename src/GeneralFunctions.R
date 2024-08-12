library('abind')

#function to concact lists (c)
conCatList = function(L1, L2){
  if(is.null(L1)){
    return(L2)
  }
  
  for(i in 1:length(L1)){
    if(is.null( dim(L2[[i]]))){
      L1[[i]] = c(L1[[i]], L2[[i]])
    } else{
      L1[[i]] = array( abind(L1[[i]], L2[[i]], along=1), dim=c(dim(L1[[i]])[[1]]+1, dim(L1[[i]])[2:length(dim(L1[[i]]))] )    )
    }
  }
  L1
}

#function to concact lists (rbind)
combineList = function(L1, L2){
  if(is.null(L1)){
    return(L2)
  }
  
  for(i in 1:length(L1)){
    if(is.null( dim(L2[[i]]))){
      if(!is.function(L1[[i]])){
        L1[[i]] = c(L1[[i]], L2[[i]])
      }
    } else{
      if(!is.function(L1[[i]])){
        L1[[i]] =  array( abind(L1[[i]], L2[[i]], along=1), dim=c(dim(L1[[i]])[[1]]+1, dim(L1[[i]])[2:length(dim(L1[[i]]))] )    )
      }}
  }
  L1
}

findCohort = function(d, depthCohorts){
  trackI =NULL
  for(i in 1:length(depthCohorts)){
     trackI=c(trackI, sum(depthCohorts[1:i]) )
  }
  which.min(abs(trackI-d))
}


biomassTransferFUN = function(transferArray, zStar, cc, ...){
  transferArray(zStar,cc, ...)
}
