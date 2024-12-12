

biomassTransferFUN = function(transferArray, inundation, cc, maxB, ...){
  transferArray(inundation,cc, maxB, ...)
}


cumSum = function(arr,d1, d2){
  if(dim(arr)[1]==1){
    return(sum(arr[1,1:d2]) )
  } else{
    if(d2>1 & d2-d1>0){
      sub.arr = arr[,d1:d2]
      return(apply(sub.arr, MARGIN=1, FUN=sum) )
    }
    else {
      return(arr[,d2])
    }
  }
}

findCohort = function(d, depthCohorts){
  trackI =NULL
  for(i in 1:length(depthCohorts)){
     trackI=c(trackI, sum(depthCohorts[1:i]) )
  }
  which.min(abs(trackI-d))
}


growRoots = function(top, bottom, rootDepthMax, shape, expDecay){
  if(shape=='lin'){
    bottom[which(bottom>rootDepthMax)]=rootDepthMax
    top[which(top>rootDepthMax)]=rootDepthMax
    slope = -2*1/(rootDepthMax^2)
    intercept = 2*1/rootDepthMax
    rootMass = intercept*(bottom-top) + slope/2*(bottom^2-top^2)
    return(rootMass)
    
  } else if (shape=='exp') {
    bottom[which(bottom>rootDepthMax)]=rootDepthMax
    top[which(top>rootDepthMax)]=rootDepthMax
    b <- expDecay / rootDepthMax
    a <- 1*(1/b*exp(b*rootDepthMax) - exp(b*rootDepthMax) * rootDepthMax-1/b)^-1
    m <- a*exp(b*rootDepthMax)
    rootMass <- a/b*(exp(b*bottom)-exp(b*top)) +  m*(top-bottom)
    return(rootMass)
  }
}

yr.slice = function(tree.state, tt){
  tree.state.yr = list(
    FineRoots = tree.state$FineRoots[,,tt],
    CoarseRoots = tree.state$CoarseRoots[,,tt],
    StrucRoots = tree.state$StrucRoots[,,tt],
    SmallRoots = tree.state$SmallRoots[,,tt],
    LeafLitter = tree.state$LeafLitter[,,tt],
    WoodLitter = tree.state$WoodLitter[,,tt],
    DeadRoot_Litter = tree.state$DeadRoot_Litter[,,tt])
  return(tree.state.yr)
}


