### Mangrove model functions

# 
# 
# #bell shaped functions --- inundation 
# flood.x = seq(0, 1, 0.01)
# flood.max = 0.89  #optimal flooding level
# #wi = 50 #inundation tolerance
# 
# flooding = function(w, flood.x, flood.max){
#   exp(-w*(flood.x-flood.max)^2 )
# }
# 
# #Bi = exp(-wi*(flood.x-flood.max)^2 )
# 
# #flood.test = sapply(seq(0, 1, 0.01), FUN=flooding, w=20, flood.max=0.89)
# #plot(flood.x,flood.test, type='l', ylab='Survival Prop', xlab='Inundation Rate')
# 
# ##
# competition = function(r, U, U.max){
#   1/(1+exp(r*(U-U.max)))
# }
# 
# 
# #comp.test = sapply(seq(1,1000,1), FUN=competition, r=-0.0092, U=540.625)
# #plot((comp.test)~seq(1,1000,1), type='l')
# 
# 
# sal_response = function(h, S, S.max){
#   1/(1+exp(h*(S-S.max)))
# }
# 

# 
# porosity = function(massAbv){
#   ifelse(massAbv>0, max(min(-0.082*log(massAbv)+0.6788, 0.9), porLim), 0.9)
# }
# 
# porosity = function(massAbv){
#   ifelse(massAbv>0, max(min(-0.055*log(massAbv)+0.7797, 0.95), porLim), 0.95)
# }
# 
# bulk.density = function(minMass){
#  pmin( ( 1.1551*minMass^0.6657 + 0.0905 ), 0.3*2.65)  #
# }
#   
# 
# bulk.density_perOM = function(perOM){
#    # pmin( ( 3.0705*perOM^-.739 ), 0.3*2.65)  # from SF bay tidal marshes (Browns, Petaluma, Rush ranch)
#     pmin( ( 3.0705*perOM^-.739 ), 0.3*2.65)
#  
  
  
 # ifelse(minMass>0, max(   0.1116+3.91*(minMass/2.65), 0.1), 0.1)
 # ifelse(minMass>0, max(  1.1551*minMass^0.6557, 0.1), 0.1)
#}


# bulk.density = function(perOM){  #Soil Physics, SSAJ vol 73, #3 2009, Ruehlmann & Kurshens
#   #B=0.93
#   soc = perOM
#   B = log(0.1)/99  #0.004  #- 3.8e-4 #0.008 error bar3.8e-4
# (.4*2.65)*exp(B*soc)  #Mg/m3 -> g/cm32.684
#  
# }


# bulk.density = function(perOM){
#  # min(5.7017*perOM^-0.78, 0.4*2.65)  #frp, 26 Pohnpei cores
#   min(6.4646*perOM^-0.819, 0.4*2.65)  #36 Pohnpei cores, 5 outliers removed, r2=0.7169
# }
# 
# 
# bulk.densityMass = function(perOM, massAbv){  #from 36 Pohnpei cores
# # min(1.3006570+log(perOM)*-0.2638220 + min(massAbv, 60)*0.0037391, 0.4*2.65)  #40.3
#   min(1.3006570+log(perOM)*-0.2638220, 0.4*2.65)
#   
#  
# }

# #Accounts for Autocompaction of soils based on mass above cohort
 
#all cores
# bulk.densityMass = function(OM, Mass){  #from Drexler cores
#   OM=OM/100
#   k1=0.06059
#   k2= 0.31428
#   k3= -0.03569
#   min(1/( (OM+k3*(Mass))/k1 + (1-OM+k3*(Mass))/k2), 1.25)
#   # min(1/( (OM+k3*log(Mass))/k1 + (1-OM+k3*log(Mass))/k2), 2.65)
# #  min( 1/( (OM/k1)+   (1-OM)/k2 + k3*log(Mass)*OM^2  ), 2.65)
# }

#only basin cores
# bulk.densityMass = function(OM, Mass){  #from Drexler cores
#   #OM=OM/100
#   k1=0.06169
#   k2= 0.16730
#   k3= -0.09615
#   Mass = min(Mass, 3)
#   max(min( 1/( (OM+k3*(Mass))/k1 + (1-OM+k3*(Mass))/k2), 1.5), 0.04)
#   # min(1/( (OM+k3*log(Mass))/k1 + (1-OM+k3*log(Mass))/k2), 2.65)
#   #  min( 1/( (OM/k1)+   (1-OM)/k2 + k3*log(Mass)*OM^2  ), 2.65)
# }


bulk.densityMass = function(OM, Mass){
  k1coef = 0.08229  
  k2coef = 1.99
  min(1/(OM/k1coef + (1-OM)/k2coef),1.9)
}

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


spCover = function(z, perCover, nicheFUN, inputPars){ 
  compA = inputPars$pars$compA
  aliveIndex = inputPars$pars$aliveIndex
  deadIndex= inputPars$pars$deadIndex
  ALPH = inputPars$pars$initCover
  r1 = inputPars$pars$r1
  seeding = inputPars$pars$seeding
  grassI = inputPars$pars$grassI
  treeI = inputPars$pars$treeI
  mixedComm = inputPars$pars$mixedComm
  inundLimReseed = inputPars$pars$inundLimReseed
  
  inundProb = pptProb=NULL
  for(i in 1:length(nicheFUN)){
    inundProb = cbind(inundProb, sapply( (z/(TR*100/2)), nicheFUN[[i]] ) )  #nicheFUN uses zStar elevation to calculate niche
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
    dC = ifelse(M>0, perCover[i]*r1[i]*M, perCover[i]*(1-Niche)*M)
    
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
        dC = max(ALPH*Niche*(CC-comp),1e-5)
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

findCohort = function(d, depthCohorts){
  trackI =NULL
  for(i in 1:length(depthCohorts)){
     trackI=c(trackI, sum(depthCohorts[1:i]) )
  }
  which.min(abs(trackI-d))
}


biomassFun = function(transferArray, zStar, cc, ...){
  transferArray(zStar,cc, ...)
}

prepFlooding = function(floodFun){
  flooding = seq(0,100, 0.1)
  relFlood =pmax(sapply(flooding, floodFun)/max(sapply(flooding, floodFun)),0)
  # relFlood[(which(relFlood==0))[1]:length(relFlood)] = 0
  approxfun(flooding,relFlood, rule=2) 
}


prepSpecies = function(speciesList, years2CC, inundFun){
  
  nicheFUN=list()
  rootDepthFUN=list()
  growth=NA
  cov2AGB=list()
  rootShoot=list()
  pars=list()
  optRSFUN=list()
  
  typeList=NULL
  for(i in 1:length(speciesList)){
    typeList = c(typeList, get(speciesList[i])$Type)
  }

  mixedComm= ifelse('Tree' %in% typeList & 'Grass' %in% typeList, 1, 0)
  
  floodTime = sapply(seq(-TR/2,TR/2,0.01)*100, inundFun) #flooding time based on elevation
  
  for(i in 1:length(speciesList)){
    sp= get(speciesList[i])
    # if(sp$Type=='Grass'){
       nicheFUN[[i]] = (sp$nicheFUN)    #prepFlooding
       #if(is.function(sp$root_shoot)){
        rootShoot[[i]] = sp$root_shoot
       #} else {
      
       z1=seq(0,100,0.1)
       optFlood= inundFun(z1[ which.max(sapply(z1, nicheFUN[[i]]))])
       optFlood = ifelse(optFlood<0, 0, optFlood)
       optRSFUN[[i]] = approxfun(floodTime, (1+sp$RSslp*(optFlood-floodTime)), rule=2, ties=mean)
    
    if(mixedComm==1){
      if(sp$Type=='Grass' | sp$Type=='Seagrass'){
        growth[i] = sp$absGrowth#*sitePars$CC#*sp$compA
      } else {
        growth[i] = sp$absGrowth#/sitePars$CC
      }
    } else {
      if(sp$Type=='Seagrass'){
        growth[i] = sp$absGrowth*sitePars$CC#*sp$compA
      } else {
        growth[i] =sp$absGrowth
      }}
    
    cov2AGB[[i]] = sp$biomassFun
    rootDepthFUN[[i]]=sp$r.depthFUN
    
    if(i==1){
      pars=sp
    } else {
      pars = combineList(pars, sp)
    }
    
  }
  
  nicheFUN[(length(nicheFUN)+1):(length(nicheFUN)*2)] = nicheFUN
  growth2= 1+(growth-mean(growth))/mean(growth)  #normalized growth rate
  
  rate= rep(NA, length(speciesList))
  for(j in 1:length(speciesList)){
    #find optimum rate for Xyears
    CC= 0.02 # any starting carrying capacity
    Cstart=.01*CC
    Abs_rate=0.01
    Xyears= years2CC[j] #years to 99% CC
    x=seq(1,2000)
    k=0
    test=seq(.02,2,.001)
    Y_CC=array(0,length(test))
    for( Abs_rate in test){
      k=k+1
      C=array(Cstart,length(x))
      for( i in 2:length(x)){
        C[i]=C[i-1]+C[i-1]*Abs_rate*(1-C[i-1]/CC)
      }
      Y_CC[k]=approx(C/CC,x,.99)$y
    }
    OUT=rbind(test,Y_CC)
    
    rate[j]=approx(Y_CC,test,Xyears)$y  #~.1 for 90 years to 99% of CC
  }
  
  pars$r1 = growth2*rate  
  pars$aliveIndex = 1:length(speciesList)
  pars$deadIndex=(length(speciesList)+1):(length(speciesList)*2)
  pars$mixedComm = mixedComm
  pars$grassI= rep(0, length(speciesList)*2)
  pars$grassI[which(c(pars$Type,pars$Type)=='Grass')]=1
  pars$treeI = which(c(pars$Type,pars$Type)=='Tree')
  pars$seagrassI = which(c(pars$Type, pars$Type)=='Seagrass')
  
  # pars$CC=1   #for marsh use CC=1. For forest, use fraction basal area (95th percentile) from inventory plots.
  
  return(list(nicheFUN=nicheFUN, coverAgTransfer=cov2AGB, pars=pars, optRSFUN=optRSFUN, rootDepthFUN=rootDepthFUN, rsFUN=rootShoot))
}




# 
# for(i in 1:5){
# print(paste(spTransLIST$baseGrowthAvg[i], numCoverFuns[[i]] ) )
#}

#######################

# 
# #
# #plot z probability functions
# z.test = seq(-75, 130, 1)
# #
# out1  = sapply(z.test, nicheFUN[[1]])
# 
# 
# 
# mcol = c('black','red','blue','brown','darkgreen')
# 
# png('C:\\MERCC\\input_Pohnpei\\Figures\\spOccurrenceProb.png', width = 200, height=200, res=300, units='mm')
# plot(z.test, out1, type='l', ylim=c(0,1), ylab='Occurrence Probability', xlab='Elevation cm MSL')
# 
# for(i in 1:5){
#   out1  = sapply(z.test, nicheFUN[[i]])
#   lines(z.test, out1, col=mcol[i])
# }
# dev.off()

# 
# 
# IF = floodingFun(elev)
# h=-17
# S=0.3
# sal_IF = function(IF) 40/(1+exp(h*(IF-S)))
# 
# o = sapply(IF, sal_IF)
# plot(o~elev)
# 
# 
# h1 = 0.00125
# S1 = 2000
# sal_shoreDist = function(shore_dist)  1/(1+exp(h1*(shore_dist-S1)))
# #o1 = sapply(shore_dist, sal_shoreDist)
# #plot(o1~shore_dist)
# 
# 
# sal = function(IF, shore_dist){
#   if(!is.na(IF)){
#     return(sal_IF(IF)*sal_shoreDist(shore_dist))
#   } else {
#     return(NA)
#   }
# }


# o2 = mapply(FUN=sal, IF=IF, shore_dist=shore_dist)
# o2
# plot(o2~IF)
# plot(o2~shore_dist)

######################

#Competition matrices for saplings and adults
# 
# compAdultMatrix = function(species.i, ntreesJ, ntreesA){
#   if(species.i<=4){ #red-A, black-A, white-A, button-A
#     ntrees_out = rowSums(ntreesJ %*% diag(c(0.1, 0.1, 0.1, 0.1, 2.5, 0, 0, 0)) + ntreesA %*% diag(c(1, 1, 1, 1, 5, 0, 0, 0)) )
#   }
#   if(species.i==5){ #pine-A
#     ntrees_out = rowSums(ntreesJ %*% diag(c(0.2, 0.2, 0.2, 0.2, 0.1, 0, 0, 0)) + ntreesA %*% diag(c(0.5, 0.5, 0.5, 0.5, 1, 0, 0, 0)) )
#   }
#   if(species.i==6){# & species.i<=7){ #FRESHM-A
#     ntrees_out = rowSums(ntreesJ %*% diag(c(800, 800, 800, 800, 800, 1, 300, 0)) + ntreesA %*% diag(c(1000, 1000, 1000, 1000 ,1000, 1, 300, 0)) )
#   }
#   if(species.i==7){ #SaltM-A
#     ntrees_out = rowSums(ntreesJ %*% diag(c(800, 800, 800, 800, 800, 300, 1, 0)) + ntreesA %*% diag(c(1000, 1000, 1000, 1000 ,1000, 300, 1, 0)) )
#   }
#   if(species.i == 8){
#     ntrees_out = rowSums(ntreesJ + ntreesA)
#   }
#   ntrees_out
# }
# 
# compJuvMatrix = function(species.i, ntreesJ, ntreesA){
#   if(species.i<=4){ #red-j, black-j, white-j, button-j
#     ntrees_out = rowSums(ntreesJ %*% diag(c(1, 1, 1, 1, 3, 0, 0, 0)) + ntreesA %*% diag(c(2, 2, 2, 2, 7, 0, 0, 0)) )
#   }
#   if(species.i==5){ #pine-j
#     ntrees_out = rowSums(ntreesJ %*% diag(c(0.5, 0.5, 0.5, 0.5, 1, 0, 0, 0)) + ntreesA %*% diag(c(0.75, 0.75, 0.75, 0.75, 2, 0, 0, 0)) )
#   }
#   if(species.i==6){# & species.i<=7){ #FRESHM-j
#     ntrees_out = rowSums(ntreesJ %*% diag(c(800, 800, 800, 800, 800, 1, 600, 0)) + ntreesA %*% diag(c(1000, 1000, 1000, 1000 ,1000, 1, 600, 0)) )
#   }
#   if(species.i==7){
#     ntrees_out = rowSums(ntreesJ %*% diag(c(800, 800, 800, 800, 800, 600, 1, 0)) + ntreesA %*% diag(c(1000, 1000, 1000, 1000 ,1000, 600, 1, 0)) )
#   }
#   if(species.i == 8){
#     ntrees_out = rowSums(ntreesJ + ntreesA)
#   }
#   ntrees_out
# }


#set water

setWater = function(otherStuff){
  if(!is.na(otherStuff[1])){
    if(otherStuff[8] == 1){
      out = rep(0, length(otherStuff)-1)
    } else{
      out = otherStuff[1:7]
    }
    return(out)
  } else{
    return(rep(NA, 7))
  }
}


#####################################










