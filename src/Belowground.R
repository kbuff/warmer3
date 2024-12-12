
# tree.state.yr = yr.slice(tree.state, tt=tt+1)
# cohorts=cohorts.state
# z= elev.time[,tt]

# Main belowground soil development function
belowground = function(tree.state.yr, cohorts, tt, step, z,  minDep, wood.den, root_por){
  
  oldCohorts=cohorts
  inund = sapply(z, inundFun)
  
  if(is.na(sum(z))){
    stop()
  }
  
  refracRemainSurf =   exp(-kdec/100*365/step*(decompSens-inund^0.5))   #fraction of pool that remains after decomposition (refractory)
  labileRemain = exp(-kdec*365/step)
  
  #decompose all
  cohorts$LOM[,1] = cohorts$LOM[,1] - cohorts$LOM[,1]*(1-labileRemain)
  cohorts$ROM[,1] = cohorts$ROM[,1]- cohorts$ROM[,1]*(1-refracRemainSurf)
  
  for(sp in 1:nSpecies){
    #calculate surface decomposition and input from wood and leaf litter
    cohorts$ROM[,1] = cohorts$ROM[,1]  + tree.state.yr$LeafLitter[,sp]*leafRefrac[sp]/step#
    cohorts$LOM[,1] = cohorts$LOM[,1] + (tree.state.yr$LeafLitter[,sp]*(1-leafRefrac[sp])/step)
    
    #wood litter moves into structural root litter pool at a given rate [wood2soil]
    cohorts$SRP[,sp,1] = cohorts$SRP[,sp,1] + cohorts$WOODLITTER[,sp,1]*wood2soil/step
    cohorts$WOODLITTER[,sp,1] =  cohorts$WOODLITTER[,sp,1] -  cohorts$WOODLITTER[,sp,1]*wood2soil/step +  tree.state.yr$WoodLitter[,sp]/step
  }
  
  #  print('decomp')
  #roots grow starting in the 2nd cohort layer
  cohortDepth = t(apply(cohorts$VOL[,2:ncohorts], MARGIN=1, FUN=cumsum))
  coefsOut <- array(NA, dim = c(nrow(cohortDepth), nSpecies, ncohorts-1))
  coefsOutBR <- array(NA, dim = c(nrow(cohortDepth), nSpecies, ncohorts-1))
  
  for(sp in 1:nSpecies){
    rootDepthT = r.roots[sp]
    rootDepthStruc=r.rootsStruc[sp]
    coefsOut[,sp,]  =  t(sapply(1:length(z), FUN=function(z2) {
      mapply(growRoots, top=c(0,cohortDepth[z2, 1:(ncohorts-2)]), bottom=cohortDepth[z2, 1:(ncohorts-1)], rootDepthMax=rootDepthT, shape=rootingShape[sp], expDecay=rootingShapeCoef[sp]) } ))
    coefsOutBR[,sp,]  = t(sapply(1:length(z), FUN=function(z2) {
      mapply(growRoots, top=c(0,cohortDepth[z2, 1:(ncohorts-2)]), bottom=cohortDepth[z2, 1:(ncohorts-1)], rootDepthMax=rootDepthStruc, shape=rootingShape[sp], expDecay=rootingShapeCoef[sp]) } ))
  }
  
  for(n in 2:ncohorts) {    #run through the cohorts, adding OM, MIN and ROOT growth, less decomp
    inund = sapply(z - cumSum(cohorts$VOL, 1,n), inundFun)
    refracRemain= exp(-kdec/100*365/step*(decompSens-inund^0.5))
    labileRemain = exp(-kdec*365/step)
    
    cohorts$CO2SoilLoss[,tt] = cohorts$CO2SoilLoss[,tt] + ((cohorts$LOM[,n])*(1-labileRemain) +
                                                             cohorts$ROM[,n]*(1-refracRemain) )*(1-inund)*c.density
    
    #lateral C transport g/cm2
    cohorts$CO2LatLoss[,tt] =  cohorts$CO2LatLoss[,tt] + ((cohorts$LOM[,n])*(1-labileRemain)+
                                                            cohorts$ROM[,n]*(1-refracRemain) )*(inund)*c.density
    
    
    #Update ROM pool.  Dead root matter [FRP, CRP] turns into particulate organic matter at a given rate [particulate; kcp, ksp, kfp]
    for(sp in 1:nSpecies){
      cohorts$LOM[,n] =   cohorts$LOM[,n]+cohorts$FRP[,sp,n]*fc1[sp]*kfp[sp]/step + cohorts$CRP[,sp,n]*fc2[sp]*kcp[sp]/step  +
        cohorts$SRP[,sp,n]*fc2[sp]*ksp[sp]/step
      
      cohorts$ROM[,n] =    cohorts$ROM[,n] + cohorts$FRP[,sp,n]*(1-fc1[sp])*kfp[sp]/step +  cohorts$CRP[,sp,n]*(1-fc2[sp])*kcp[sp]/step  +
        cohorts$SRP[,sp,n]*(1-fc2[sp])*ksp[sp]/step
      
      #add in this years turnover
      cohorts$FRP[,sp,n] = cohorts$FRP[,sp,n]+ kf[sp]*(1-coefsOut[,sp,n-1])*cohorts$FINEROOTS[,sp,n]/step + km[sp]*(1-coefsOut[,sp,n-1])*cohorts$SMALLROOTS[,sp,n]/step -
        cohorts$FRP[,sp,n]*(kfp[sp])/step
      
      cohorts$CRP[,sp,n] = cohorts$CRP[,sp,n]+  kc[sp]*(1-coefsOut[,sp,n-1])*cohorts$COARSEROOTS[,sp,n]/step -  cohorts$CRP[,sp,n]*(kcp[sp])/step
      cohorts$SRP[,sp,n]= cohorts$SRP[,sp,n]+ ks[sp]*(1-coefsOutBR[,sp,n-1])*cohorts$STRUCROOTS[,sp,n]/step - cohorts$SRP[,sp,n]*(ksp[sp])/step  + cohorts$WOODLITTER[,sp,n]*wood2soil/step
      
      cohorts$WOODLITTER[,sp,n] = cohorts$WOODLITTER[,sp,n] - cohorts$WOODLITTER[,sp,n]*wood2soil/step
      
      #update live root mass after turnover
      cohorts$FINEROOTS[,sp,n]   = pmax(1e-8, cohorts$FINEROOTS[,sp,n] -  kf[sp]*(1-coefsOut[,sp,n-1])*cohorts$FINEROOTS[,sp,n]/step )
      cohorts$SMALLROOTS[,sp,n]  = pmax(1e-8, cohorts$SMALLROOTS[,sp,n]-  km[sp]*(1-coefsOut[,sp,n-1])*cohorts$SMALLROOTS[,sp,n]/step )
      cohorts$COARSEROOTS[,sp,n] = pmax(1e-8, cohorts$COARSEROOTS[,sp,n]- kc[sp]*(1-coefsOut[,sp,n-1])*cohorts$COARSEROOTS[,sp,n]/step)
      cohorts$STRUCROOTS[,sp,n]  = pmax(1e-8, cohorts$STRUCROOTS[,sp,n] - ks[sp]*(1-coefsOutBR[,sp,n-1])*cohorts$STRUCROOTS[,sp,n]/step )
    }
    
    cohorts$LOM[,n] = cohorts$LOM[,n]- cohorts$LOM[,n]*(1-labileRemain)
    cohorts$ROM[,n] = cohorts$ROM[,n]- cohorts$ROM[,n]*(1-refracRemain)
    cohorts$CARBON[,n] = c.density*( cohorts$LOM[,n] + cohorts$ROM[,n] + rowSums(cohorts$FRP[,,n] + cohorts$CRP[,,n] + cohorts$SRP[,,n] + cohorts$FINEROOTS[,,n] + cohorts$SMALLROOTS[,,n] + cohorts$COARSEROOTS[,,n] + cohorts$STRUCROOTS[,,n] + cohorts$WOODLITTER[,,n]))
  }
  
  for(sp in 1:nSpecies){
    #living roots below the rooting depth are moved to dead root pools. Allows plant to be buried by a sediment pulse
    deepRoots = (coefsOut[,sp,1:(ncohorts-1)]*0)
    VLdeepRoots = (coefsOutBR[,sp,1:(ncohorts-1)]*0)
    
    #distribute live roots downcore. effectively re-grow each timestep
    cohorts$FINEROOTS[,sp,2:ncohorts] =    coefsOut[,sp,] * tree.state.yr$FineRoots[,sp]
    cohorts$COARSEROOTS[,sp,2:ncohorts] =  coefsOut[,sp,] * tree.state.yr$CoarseRoots[,sp]
    cohorts$SMALLROOTS[,sp,2:ncohorts] =   coefsOut[,sp,] * tree.state.yr$SmallRoots[,sp]
    cohorts$STRUCROOTS[,sp,2:ncohorts] =   coefsOutBR[,sp,]*tree.state.yr$StrucRoots[,sp]
    
    #living roots below the live root zone are moved to root litter. Add in dead roots [from aboveground tree death] and woodlitter
    cohorts$SRP[,sp,2:ncohorts] = cohorts$SRP[,sp,2:ncohorts] + tree.state.yr$DeadRoot_Litter[,sp]*strucRoots[sp]*coefsOutBR[,sp,]/step + cohorts$STRUCROOTS[,sp,2:ncohorts]*VLdeepRoots + cohorts$WOODLITTER[,sp,2:ncohorts]*wood2soil/step
    cohorts$CRP[,sp,2:ncohorts]  = cohorts$CRP[,sp,2:ncohorts]  + tree.state.yr$DeadRoot_Litter[,sp]*(1-strucRoots[sp])*(coarseRoots[sp])*coefsOut[,sp,]/step  + cohorts$COARSEROOTS[,sp,2:ncohorts]*deepRoots
    cohorts$FRP[,sp,2:ncohorts]  = cohorts$FRP[,sp,2:ncohorts]  + coefsOut[,sp,]*tree.state.yr$DeadRoot_Litter[,sp]*(1-strucRoots[sp])*(fineRoots[sp])/step  +
      coefsOut[,sp,]*tree.state.yr$DeadRoot_Litter[,sp]*(1-strucRoots[sp])*(smallRoots[sp])/step  + cohorts$SMALLROOTS[,sp,2:ncohorts]*deepRoots + cohorts$FINEROOTS[,sp,2:ncohorts]*deepRoots
    
    cohorts$FINEROOTS[,sp,2:ncohorts]    = cohorts$FINEROOTS[,sp,2:ncohorts]*(1-deepRoots)
    cohorts$COARSEROOTS[,sp,2:ncohorts]  =  cohorts$COARSEROOTS[,sp,2:ncohorts]*(1-deepRoots)
    cohorts$SMALLROOTS[,sp,2:ncohorts]   = cohorts$SMALLROOTS[,sp,2:ncohorts]*(1-deepRoots)
    cohorts$STRUCROOTS[,sp,2:ncohorts]   = cohorts$STRUCROOTS[,sp,2:ncohorts]*(1-VLdeepRoots)
    
    #root mass cant be negative
    cohorts$FINEROOTS[,sp,2:ncohorts]   = pmax(0,  cohorts$FINEROOTS[,sp,2:ncohorts])
    cohorts$COARSEROOTS[,sp,2:ncohorts] = pmax(0,  cohorts$COARSEROOTS[,sp,2:ncohorts])
    cohorts$SMALLROOTS[,sp,2:ncohorts]  = pmax(0,  cohorts$SMALLROOTS[,sp,2:ncohorts])
    cohorts$STRUCROOTS[,sp,2:ncohorts]  = pmax(0,  cohorts$STRUCROOTS[,sp,2:ncohorts])
  }
  
  
  #Surface Deposition (mineral deposition & wood litter:aboveground production & fallen dead)
  #some faction of deposited sediment is organic (e.g., 10%)
  cohorts$MIN[,1] = cohorts$MIN[,1]+ (minDep)/step
  cohorts$ROM[,1] = cohorts$ROM[,1] + (minDep*surfOMDep)/step
  
  cohorts = updateDepthPOR(cohorts,oldCohorts, wood.den, root.den, root_por, porosityRate = porosityRate, minPorosity = minPorosity, maxPorosity=maxPorosity, maxDeltaPor=maxDeltaPor, step,Initialize=F, SPINUP=F)

  return(cohorts)
}
