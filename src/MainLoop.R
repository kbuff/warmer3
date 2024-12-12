## Main Loop

#mass units g/cm2. Cohort assumed to have area of 1 cm2

MainLoop_wSPINUP = function(elev, years, SLR, spinupYrs, inputPars, sitePars){
  
  {
    #make model parameters global [available outside function]
    strucRoots <<- inputPars$strucRoots
    smallRoots <<- inputPars$smallRoots
    coarseRoots <<- inputPars$coarseRoots
    fineRoots <<- inputPars$fineRoots
    initCover <<- inputPars$initCover
    deadFallRate <<- inputPars$deadFallRate
    CC <<- sitePars$CC
    r.roots <<- inputPars$r.roots
    r.rootsStruc <<- inputPars$r.rootsStruc
    rootingShape <<- inputPars$rooting_shape
    rootingShapeCoef <<- inputPars$rooting_shape_coef
    kf <<- inputPars$kf
    km <<- inputPars$km
    kc <<- inputPars$kc
    ks <<- inputPars$ks
    kfp <<- inputPars$kfp
    kcp <<- inputPars$kcp
    ksp <<- inputPars$ksp
    fc1 <<- inputPars$fc1
    fc2 <<- inputPars$fc2
    leafRefrac <<- inputPars$leafRefrac
    leafProp <<- inputPars$leafProp
    minPorosity <<- sitePars$minPorosity
    porosityRate <<- sitePars$porosityRate
    maxPorosity <<- sitePars$maxPorosity
    maxDeltaPor <<- sitePars$maxDeltaPor
    kdec <<- sitePars$kdec
    decompSens <<- sitePars$decompSens
    wood2soil <<- sitePars$wood2soil
    minSedzStar <<- sitePars$minSedzStar
    litterDep <<- inputPars$litterDep
    litterRetain <<- inputPars$litterRetain
    surfOMDep <<- sitePars$surfOMDep
    root.den <<- inputPars$root.den
    wood.den <<- inputPars$wood.den
    root_por <<- inputPars$root_por
    ssc <<- sitePars$ssc
    
    root_shoot <<- inputPars$rsFUN
    biomassFUN <<- inputPars$biomassFUN
    nicheFUN <<- inputPars$nicheFUN
    rootDepthFUN = inputPars$rootDepthFUN
    aliveIndex <<- inputPars$aliveIndex
    deadIndex <<- inputPars$deadIndex
  
    
    #depth cohort intervals. Calibration will be specific to the intervals selected. 
    #The last cohort will accumulate volume.
    z.depth <<- c(5, 10, 10, 10, 10, 10, 10, 20, 20, 40, 50)
    ncohorts <<- length(z.depth)
    
    
    elev.time =acc = array(dim=c(length(elev), years))
    
    elev.time[,1] = elev
    perCover = array(0, dim=c(length(elev), nSpecies*2, years))
    
    initialCover = array( initCover*CC, dim=c(length(elev), nSpecies*2 )  )
    initialCover[,deadIndex] = 0
    seeding = inputPars$pars$seeding
    age = array(1, dim=c(length(elev),nSpecies))
    perCover[,,1] = initialCover
    
    slr_yrs  =array(0, dim=c(years))
    soilDepth = totMIN=totOM=totROOT_VOL = totROOTS =fineROOTS=coarseROOTS=strucROOTS=smallROOTS=totDEADROOTS=totCARBON= totCARBON_Abv = LOMVol = ROMVol = MINVol = LiveROOTVOL =
      allMIN= DeadROOTVOL = totVol = bd_top1 =bd_top2= bd_top3 = bd= allOM =  totDOWNEDWOOD =topDepth=totLOM = totROM= woodlitter = leaflitter=woodDen=
      totDOWNEDWOOD_V=deadRootROM=deadRootLOM = deadRootROM_V=deadRootLOM_V= totDEADROOTSALL= newDeadRoots=totBelowB =array(0, dim=c(length(elev), years))
    belowgroundBiomass =  Root_Shoot= array(0, dim = c(length(elev),nSpecies*2, years))
    BelowgroundB_AdultLIVING=BelowgroundB_AdultDEAD =BelowgroundB_Dead= array(0, dim=c(length(elev),nSpecies, 2))
    ABGBiomass=array(0, dim=c(length(elev),nSpecies, years))
    
  }
  
  # start of MAIN LOOP
  for(tt in 1:(years-1)){
    {
      # print(tt)
      if(tt==1){        #Initialize soil cohort
        
        tree.state = list(                                            #example size distribution definition. Doesnt impact model, just different source for parameters. Model can also be run with 1 root class
          FineRoots = array(0, dim=c(length(elev),nSpecies, years)), #<=2mm
          SmallRoots = array(0, dim=c(length(elev),nSpecies, years)), #<2cm, >=2mm
          CoarseRoots = array(0, dim=c(length(elev),nSpecies, years)), #>2 cm, <5 cm
          StrucRoots = array(0, dim=c(length(elev),nSpecies, years)), #>5 cm
          LeafLitter = array(0, dim=c(length(elev), nSpecies,years)),
          WoodLitter = array(0, dim=c(length(elev), nSpecies,years)),
          DeadRoot_Litter = array(0, dim=c(length(elev),nSpecies, years))
        )
        
        cohorts.state = list(
          LOM = array(rep(0, each=length(elev)), dim=c(length(elev),ncohorts)),         #g cm-2  labile OM
          ROM =  t(array(z.depth*(sitePars$perOMCohorts)*1 ,dim=c(ncohorts, length(elev)))),        #g cm-2  mineral deposition (surface)
          SAND = array(rep(0*z.depth, each=length(elev)), dim=c(length(elev), ncohorts)),
          MIN =  t(array(z.depth*(1-sitePars$perOMCohorts)*1 ,dim=c(ncohorts, length(elev)))),        #g cm-2  mineral deposition (surface)
          POR = array(0.98, dim=c(length(elev), ncohorts)),  # 0-1
          WOODLITTER = array(0.0, dim=c(length(elev),nSpecies,ncohorts)),          #g cm-2  dead woody litter, surface only. slowly leaches to belowground cohorts
          FINEROOTS= array(0, dim=c(length(elev),nSpecies,ncohorts)),    #g cm-2  live roots
          SMALLROOTS= array(0, dim=c(length(elev),nSpecies,ncohorts)),     #g cm-2  live roots
          COARSEROOTS= array(0, dim=c(length(elev),nSpecies,ncohorts)),     #g cm-2  live roots
          STRUCROOTS= array(0, dim=c(length(elev),nSpecies,ncohorts)),     #g cm-2  live roots
          FRP =  array(0, dim=c(length(elev),nSpecies,ncohorts)),        #g cm-2  fine dead roots
          CRP =  array(0, dim=c(length(elev),nSpecies,ncohorts)),        #g cm-2  coarse dead roots
          SRP =  array(0, dim=c(length(elev),nSpecies,ncohorts)),        #g cm-2  structural dead roots
          CARBON = array(0, dim=c(length(elev),ncohorts)),      #g cm-2  Carbon
          VOL = t(array(rep(z.depth, length(elev)), dim=c(ncohorts,length(elev)))),         #cm3     total cohort volume
          BD =  array(0, dim=c(length(elev),ncohorts)),
          CO2SoilLoss = array(0, dim=c(length(elev), years)),
          CO2LatLoss = array(0, dim=c(length(elev), years))
        )
        cohorts.state$Vvoid = cohorts.state$VOL
        cohorts.state = updateDepthPOR(cohorts.state,cohorts.state, wood.den, root.den, root_por, porosityRate, minPorosity, maxPorosity, maxDeltaPor, 1,T,T)
        soilDepth[,1] = rowSums(cohorts.state$VOL)
      }
      
      step=6
      for(s in 1:step){
        
        # LAG SP  PRODUCTIVITY~ELEVATION RESPONSE
        #   if(tt>5){
        #     # print(elev.time[1,t:(t-4)])
        #      z_lag = sapply(1:length(elev), FUN=function(x) mean(elev.time[x,tt:(tt-4)]) )
        #    } else {
        z_lag = t(t(elev.time[,tt]))
        #   }
        
        if(s==1){
          temp1 = t(sapply(1:length(z_lag),  FUN=function(x) spCover(z_lag[x],  perCover[x,,tt], nicheFUN,  inputPars, step ) ))
        } else {
          temp1 = t(sapply(1:length(z_lag),  FUN=function(x) spCover(z_lag[x],  perCover[x,,tt+1], nicheFUN,  inputPars, step ) ))
        }
       
        ##########################################################
        #Belowground Processes ###################################
        ##########################################################
    
        floodDur = sapply(elev.time[,tt], inundFun)
        zStar = elev.time[,tt]/(TR*100/2)
        
        #calculate the aboveground dead biomass using the transfer function. Only incorporate the newly dead trees
        if(s==1){
          perCover[,,tt+1] = t(array(unlist(temp1[,1]), dim=c(nSpecies*2, length(elev)) ))
          fallenDead = t(array(unlist(temp1[,2]), dim=c( nSpecies*2, length(elev))))  #as %cover  #FALLEN DEAD -- LITTER from DEAD
          
          # remove dead trees from t=2 if no trees can survive in t=1
          if(tt==1){
            perCover[,,tt] = perCover[,,tt+1]
            for(i in 1:length(aliveIndex)){
              for(k in 1:length(elev) ){
                if( perCover[k,i,tt] < initCover[i]*sitePars$CC ){
                  # print('reset')
                  perCover[k,i,tt] = 0
                  perCover[k,i,tt+1] = 0
                  perCover[k,i+length(aliveIndex),tt+1] = 0
                  perCover[k,i+length(aliveIndex),tt] = 0
                }
              }}
          }
          
          AbovegroundB_AdultDEAD =  t(sapply(1:length(elev), perCover,  FUN=function(x, perCover) {mapply(biomassTransferFUN, transferArray=biomassFUN, maxB=inputPars$spMaxBiomass, inundation=floodDur[x],cc = perCover[x,,tt+1]-perCover[x,,tt]+rep(fallenDead[x,],2)  ) } ) )[,deadIndex] # g  only include the newly dead trees
        } else {
          newCover = t(array(unlist(temp1[,1]), dim=c(nSpecies*2, length(elev)) ))
          AbovegroundB_AdultDEAD =  t(sapply(1:length(elev), perCover,  FUN=function(x, perCover) {mapply(biomassTransferFUN, transferArray=biomassFUN, maxB=inputPars$spMaxBiomass,inundation=floodDur[x],cc = newCover[x,]-perCover[x,,tt+1]+rep(fallenDead[x,],2)  ) } ) )[,deadIndex] # g  only include the newly dead trees
          perCover[,,tt+1]=newCover
          fallenDead = t(array(unlist(temp1[,2]), dim=c( nSpecies*2, length(elev))))  #as %cover  #FALLEN DEAD -- LITTER from DEAD
        }
        
        #calculate aboveground biomass using a transfer function of elevation [zStar, z/(TR/2)]
        AbovegroundB_Adult= t(sapply(1:length(elev), zStar,perCover, FUN=function(x, zStar, perCover) { mapply(biomassTransferFUN, transferArray=biomassFUN, maxB=inputPars$spMaxBiomass, inundation=floodDur[x], cc=perCover[x,,tt+1] ) } ) ) # g
        
        #use root:shoot ratio to calculate living root biomass (g/cm2)
        rs=NULL
        for(sp in 1:nSpecies){
          if(inputPars$rsFUN_Inund[[sp]]){
            rs= cbind(rs, sapply(floodDur, FUN=root_shoot[[sp]], rs=inputPars$root_shoot[sp]))
          } else {
            rs=cbind(rs, rep(inputPars$root_shoot[sp], length(elev)))
          }
        }
      
        BelowgroundB_AdultLIVING[,,2]= rs*AbovegroundB_Adult[,aliveIndex]
        totBelowB[,tt] = sapply(1:length(elev), FUN=function(x) sum(BelowgroundB_AdultLIVING[x,,2]))
        
        AbovegroundB_Dead  =t(sapply(1:length(elev), fallenDead,  FUN=function(x, fallenDead) { mapply(biomassTransferFUN, transferArray=biomassFUN, maxB=inputPars$spMaxBiomass, inundation=floodDur[x],cc = fallenDead[x,] ) } ) ) # g
        # BelowgroundB_Dead[,,2] = t(sapply(1:length(elev), fallenDead, FUN = function(z, fallenDead)  mapply( FUN=function(x,y){x*y}, y=root_shoot, x= fallenDead[z,] )))
        ########################
        
        #new root growth gets divided into size classes-- root production already accounts for sub-timestep.
        newGrowth = array(BelowgroundB_AdultLIVING[,,2]-BelowgroundB_AdultLIVING[,,1], dim=c(length(elev), nSpecies))
        newDead = newGrowth
        newDead =  rowSums(t(apply(newDead*-1, MARGIN=1, FUN=pmax, 0))) #(0, (newDead*-1))
        
        standingBGBiomass= array(BelowgroundB_AdultLIVING[,,2], dim=c(length(elev), nSpecies))
        tree.state$StrucRoots[,,tt+1] =   (sapply(1:nSpecies, y=standingBGBiomass, FUN=function(x, y) {mapply(FUN=function(x1,y1) {x1*y1}, y1=strucRoots[x], x1=y[,x]  )} ) )  #totBelowB*strucRoots
        tree.state$FineRoots[,,tt+1]  =   (sapply(1:nSpecies, y=standingBGBiomass, FUN=function(x, y) {mapply(FUN=function(x1,y1) {x1*y1}, y1=(1-strucRoots[x])*fineRoots[x], x1=y[,x]  )} ) )  #totBelowB*(1-strucRoots)*fineBigRoots
        tree.state$CoarseRoots[,,tt+1]   = (sapply(1:nSpecies, y=standingBGBiomass, FUN=function(x, y) {mapply(FUN=function(x1,y1) {x1*y1}, y1=(1-strucRoots[x])*(coarseRoots[x]), x1=y[,x]  )} ) )   # totBelowB*(1-strucRoots)*(coarseRoots)
        tree.state$SmallRoots[,,tt+1]   = (sapply(1:nSpecies, y=standingBGBiomass, FUN=function(x, y) {mapply(FUN=function(x1,y1) {x1*y1}, y1=(1-strucRoots[x])*(smallRoots[x]), x1=y[,x]  )} ) )   # totBelowB*(1-strucRoots)*(smallRoots)
        
        #update tracking
        BelowgroundB_AdultLIVING[,,1]=BelowgroundB_AdultLIVING[,,2]
        
        #Leaf litter production (to surface pool)
        for(sp in 1:nSpecies){
          tree.state$LeafLitter[,sp,tt+1] =     (sapply(1:length(elev), FUN=function(x) sum(AbovegroundB_Adult[x,sp]))*litterDep[sp]*(leafProp[sp])*litterRetain[sp]*(1-sapply(elev, FUN=inundFun)))
          #Woody litter production (alive & dead)  Only litter production from live trees. Aboveground dead biomass does NOT get included in soil
          tree.state$WoodLitter[,sp,tt+1] =  ((sapply(1:length(elev), FUN=function(x) (AbovegroundB_Adult[x,sp]))*(litterDep[sp]*(1-leafProp[sp])))* litterRetain[sp]*(1-sapply(elev, inundFun)))  +
            (sapply(1:length(elev), FUN=function(x) (AbovegroundB_Dead[x,sp]))*litterRetain[sp])
        }
        tree.state$DeadRoot_Litter[,,tt+1] = newDead #standing dead + fallen dead  
        
        #track fallen dead tree volume, but not incorporated into soil.
        totDOWNEDWOOD_V[,tt+1] =NA# rowSums(t(sapply(1:length(elev), AbovegroundB_Dead,  FUN=function(x, AbovegroundB_Dead) { mapply( FUN=function(x,y){x/y}, y=wood.den, x= AbovegroundB_Dead[x,]) } )))
        
        #set minimum level of sediment deposition, here equivalent to minSedzStar
        minDep =  pmax(sapply(FUN=sedinundFun, z1=minSedzStar*(TR/2)*100, rowSums(AbovegroundB_Adult)*10000)*ssc,  mapply(sedinundFun, z1=elev.time[,tt], B=rowSums(AbovegroundB_Adult)*10000)*ssc)#, dim=c(length(z)) )#*minScalor  #*10000 #mineral deposition, g cm-2, [based on inundation frequency from tidal harmonics for one yr (15 min resolution)]
        
        totRootMass =  apply(cohorts.state$FINEROOTS + cohorts.state$SMALLROOTS + cohorts.state$COARSEROOTS + cohorts.state$STRUCROOTS , FUN=sum, MARGIN=c(1,2))
        
        #update soil cohorts -- root growth, decomposition, bulk density & volume
        cohorts.state1 <-  belowground( yr.slice(tree.state, tt+1), cohorts.state, tt+1,  step, elev.time[,tt], minDep, wood.den, root_por)
        
        accStep =  rowSums(cohorts.state1$VOL)- soilDepth[,tt]
        if(is.na(sum(accStep))){
          stop()
        }
        cohorts.state=cohorts.state1
        soilDepth[,tt] = rowSums(cohorts.state$VOL)
        
        if(tt<=spinupYrs){
          elev.time[,tt] =  elev.time[,1]
        } else {
          elev.time[,tt] = elev.time[,tt] + accStep -  SLR[tt,]/step
        }
        
      }  #end of substep loop
      
      totRootMass1 = apply(cohorts.state$FINEROOTS + cohorts.state$SMALLROOTS + cohorts.state$COARSEROOTS + cohorts.state$STRUCROOTS, FUN=sum, MARGIN=c(1,2)) #sum of roots at each elevation for each species
      
      #root mass loss due to deep root death
      BGBMassLoss= totRootMass1-totRootMass
      covLost = (perCover[,aliveIndex,tt+1] * BGBMassLoss/totRootMass)*-1
      
      #if roots have been lost due to deep root death, update the percent cover
      covLost = pmax(0, covLost)
      covLost = ifelse(is.nan(covLost), 0, covLost)
      perCover[,aliveIndex,tt+1]= perCover[,aliveIndex,tt+1] - covLost
      perCover[,deadIndex,tt+1] = perCover[,deadIndex,tt+1] + covLost
   
      elev.time[,tt+1] = elev.time[,tt]
      soilDepth[,tt+1] = rowSums(cohorts.state$VOL)
      
      targetDepth=  apply(cohorts.state$VOL, MARGIN=1, FUN=findCohort, d=100) 
      slr_yrs[tt]  = SLR[tt,1]
      ABGBiomass[,,tt] = AbovegroundB_Adult[,1:nSpecies]
      totMIN[,tt]=  sapply(1:length(elev), FUN= function(x) sum(cohorts.state$MIN[x,1:targetDepth[x]]  )  )
      totOM[,tt] =  sapply(1:length(elev), FUN= function(x) sum(cohorts.state$LOM[x,1:targetDepth[x]]) + sum(cohorts.state$ROM[x,1:targetDepth[x]])  )
      totDEADROOTS[,tt] = sapply(1:length(elev), FUN= function(x) {sum(cohorts.state$FRP[x,,1:targetDepth[x]]) + sum(cohorts.state$CRP[x,,1:targetDepth[x]]) +
          sum(cohorts.state$SRP[x,,1:targetDepth[x]]) })
      
      totROOTS[,tt] = sapply(1:length(elev), FUN= function(x) sum(cohorts.state$STRUCROOTS[x,,1:targetDepth[x]]) ) +
        sapply(1:length(elev), FUN= function(x) sum(cohorts.state$COARSEROOTS[x,,1:targetDepth[x]]) )+
        sapply(1:length(elev), FUN= function(x) sum(cohorts.state$SMALLROOTS[x,,1:targetDepth[x]]) )+
        sapply(1:length(elev), FUN= function(x) sum(cohorts.state$FINEROOTS[x,,1:targetDepth[x]]) )
      fineROOTS[,tt] =    sapply(1:length(elev), FUN= function(x) sum(cohorts.state$FINEROOTS[x,,]) )
      coarseROOTS[,tt] =    sapply(1:length(elev), FUN= function(x) sum(cohorts.state$COARSEROOTS[x,,]) )
      strucROOTS[,tt] =    sapply(1:length(elev), FUN= function(x) sum(cohorts.state$STRUCROOTS[x,,]) )
      smallROOTS[,tt] =    sapply(1:length(elev), FUN= function(x) sum(cohorts.state$SMALLROOTS[x,,]) )
      
      totCARBON[,tt] = sapply(1:length(elev), FUN= function(x) sum(cohorts.state$CARBON[x,]) )  #entire core
      totCARBON_Abv[,tt] = sapply(1:length(elev), FUN= function(x) sum(AbovegroundB_Adult[x]) )*c.density
      totDOWNEDWOOD[,tt] = sapply(1:length(elev), FUN= function(x)  sum(tree.state$WoodLitter[x,,tt]) ) 
      leaflitter[,tt] = rowSums(tree.state$LeafLitter[,,tt])
      woodlitter[,tt] = rowSums(tree.state$WoodLitter[,,tt])
      woodDen[,tt] = mean(wood.den)
      totVol[,tt] =  sapply(1:length(elev), FUN= function(x) sum(cohorts.state$VOL[x,1:targetDepth[x]]) )
      
      allOM[,tt] =  sapply(1:length(elev), FUN= function(x) sum(cohorts.state$LOM[x,] +cohorts.state$ROM[x,] + sum(cohorts.state$FINEROOTS[x,,]) + sum(cohorts.state$FRP[x,,]) + sum(cohorts.state$CRP[x,,]) ))
      allMIN[,tt] =  sapply(1:length(elev), FUN= function(x) sum(cohorts.state$MIN[x,]) )
      
      LOMVol[,tt] = sapply(1:length(elev),    FUN = function(x) sum(cohorts.state$LOM[x,1:targetDepth[x]]/cohorts.state$BD[x,1:targetDepth[x]] ))
      ROMVol[,tt] = sapply(1:length(elev),    FUN = function(x) sum(cohorts.state$ROM[x,1:targetDepth[x]]/cohorts.state$BD[x,1:targetDepth[x]] ))
      MINVol[,tt] = sapply(1:length(elev),    FUN = function(x) sum(cohorts.state$MIN[x,1:targetDepth[x]]/cohorts.state$BD[x,1:targetDepth[x]] ))
      LiveROOTVOL[,tt] =  sapply(1:length(elev), FUN= function(x) sum(  cohorts.state$SMALLROOTS[x,,1:targetDepth[x]]/root.den + cohorts.state$COARSEROOTS[x,,1:targetDepth[x]]/root.den +  cohorts.state$FINEROOTS[x,,1:targetDepth[x]]/root.den +  cohorts.state$STRUCROOTS[x,,1:targetDepth[x]]/wood.den  ))
      DeadROOTVOL[,tt] =  sapply(1:length(elev), FUN= function(x) sum( cohorts.state$FRP[x,,1:targetDepth[x]]/(root.den/(1-root_por))  +  cohorts.state$CRP[x,,1:targetDepth[x]]/(root.den/(1-root_por)) + cohorts.state$SRP[x,,1:targetDepth[x]]/(wood.den/(1-root_por))))
      
      totLOM[,tt] =  sapply(1:length(elev), FUN= function(x) sum(cohorts.state$LOM[x,1:targetDepth[x]]) )
      totROM[,tt] =  sapply(1:length(elev), FUN= function(x) sum(cohorts.state$ROM[x,1:targetDepth[x]]) )
      topDepth[,tt] =  sapply(1:length(elev), FUN= function(x) sum(cohorts.state$VOL[x,1:targetDepth[x]]) )
      bd[,tt] =  sapply(1:length(elev), FUN= function(x) mean(cohorts.state$BD[x,1:targetDepth[x]]) )
      bd_top1[,tt] = sapply(1:length(elev), FUN= function(x) (cohorts.state$BD[x,1]) )
      bd_top2[,tt] = sapply(1:length(elev), FUN= function(x) (cohorts.state$BD[x,2]) )
      bd_top3[,tt] = sapply(1:length(elev), FUN= function(x) (cohorts.state$BD[x,3]) )
      
      if(tt==1){
        acc[,tt] = elev.time[,tt+1]-elev # rep(0,length(elev))
      } else {
        acc[,tt] = (  soilDepth[,tt] - soilDepth[,tt-1])
      }
    }
    
  }
  
  
  #output model results
  mercc.out = list(
    inputPars = inputPars,
    SLR = SLR,
    slr_yrs = slr_yrs,
    elev_msl = elev.time,
    soilDepth = soilDepth,
    acc = acc,
    tree.state=tree.state,
    cohorts.state=cohorts.state,
    min = totMIN,
    om = totOM,
    totABG = ABGBiomass,
    totBGB = totBelowB,
    totVol = totVol,
    roots = totROOTS,
    rootShoots =Root_Shoot,
    fineRoots = fineROOTS,
    smallRoots=smallROOTS,
    coarseRoots=coarseROOTS,
    strucRoots=strucROOTS,
    root_vol = totROOT_VOL,
    deadRoots = totDEADROOTS,
    deadRootVol = DeadROOTVOL,
    liveRootVol = LiveROOTVOL,
    romVol = ROMVol,
    lomVol = LOMVol,
    minVol = MINVol,
    allOM = allOM,
    allMIN = allMIN,
    
    carbon = totCARBON,
    carbon_abv=totCARBON_Abv,
    perCover = perCover,
    downedWood = totDOWNEDWOOD,
    downedWood_V = totDOWNEDWOOD_V,
    leaflitter =leaflitter,
    woodlitter = woodlitter,
    woodDen = woodDen,
    lom = totLOM,
    rom = totROM,
    topDepth = topDepth,
    bd=bd,
    bd1 = bd_top1,
    bd2 = bd_top2,
    bd3 = bd_top3
    
  )
  mercc.out
}


