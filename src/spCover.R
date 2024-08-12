
spCover = function(z, perCover, nicheFUN, inputPars, step){ 
  compA = inputPars$pars$compA 
  aliveIndex = inputPars$pars$aliveIndex
  deadIndex= inputPars$pars$deadIndex
  ALPH = inputPars$pars$initCover
  r1 = inputPars$pars$absGrowth
  seeding = inputPars$pars$seeding
  grassI = inputPars$pars$grassI
  treeI = inputPars$pars$treeI
  mixedComm = inputPars$pars$mixedComm
  inundLimReseed = inputPars$pars$inundLimReseed
  
  inundProb =NULL
  for(i in aliveIndex){
    inundProb = cbind(inundProb, sapply( (z/(TR*100/2)), nicheFUN[[i]] ) )  #nicheFUN uses zStar elevation to calculate niche probability [0-1]
  }
  
  trackAlive=NULL
  trackDead = NULL
  fallenDead = array(0, dim=length(deadIndex))
  
  if(length(aliveIndex)>1){
    shuffled = sample(aliveIndex, length(aliveIndex))
  } else {
    shuffled=1
  }
  for(i in shuffled){
    Niche = inundProb[i]
    
    otherSps = perCover[aliveIndex]* compA
    otherSps[i] = 0 
    comp =  sum(otherSps) #amount of carrying capacity taken up by other living species
    
    dC=NULL
    M = max(1-((perCover[i]) / max( (CC-comp)*Niche, 1e-12))  ,-1) 
    #growth                 #death
    dC = ifelse(M>0, perCover[i]*r1[i]*M/step, perCover[i]*(1-Niche)*M/step)
    
    if(is.na(dC)){
      dC=0
    }
    
    if(dC>=0){
      trackDead= c(trackDead, 0)
      trackAlive = c(trackAlive, dC)
    } else {
      trackDead = c(trackDead, dC)
      trackAlive = c(trackAlive, 0)
    }
    
    #RESEEDING
    if(abs(dC)<(ALPH[i]*CC*r1[i])){
      if( perCover[i] < (ALPH[i]*Niche*(CC-comp))) {
        perCover[i] = 0 #max(ALPH**zProb[i]*(CC-comp),1e-6)
        dC = max(ALPH[i]*Niche*(CC-comp),1e-5)
        if(inundFun(z) > inundLimReseed[i]){ #no reseeding above 95% flood time
          dC=0
        }
      } 
    }
    
    #fall trees first
    fallenDead[i] = perCover[i+length(aliveIndex)]*deadFallRate[i]    #fall rate of dead trees
    perCover[i+length(aliveIndex)] =  perCover[i+length(aliveIndex)]- fallenDead[i]
    
    #then update standing dead
    perCover[c(i,i+length(aliveIndex))] = perCover[c(i,i+length(aliveIndex))] + c(dC, ifelse(dC<=0, abs(dC), 0))  #UPDATING COVER
    
  }
  
  return(list(perCover=perCover, decomp=fallenDead)) 
  
}