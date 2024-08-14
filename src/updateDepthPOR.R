
#Calculate bulk density-- only used for the surface cohort
calc.BD = function(cohorts.vol, n){
  if(n>1){
    g = cumSum(cohorts.vol$OM,1,n-1) + cumSum(cohorts.vol$MIN,1, n-1) #sum of the mass of everything above cohort n
    mass = pmax(cohorts.vol$OM[,n]+ cohorts.vol$MIN[,n], 1e-8)
    pOM_n =  sapply(1:length(mass), function(x)  ifelse(mass[x]>0, (cohorts.vol$OM[x,n])/mass[x], 0) )
    pOM_n = pmin(pOM_n,1)
    BD =  mapply( bulk.densityMass, pOM_n, (g))
  } else {
    g = 0
    mass = pmax(cohorts.vol$OM[,n]+cohorts.vol$MIN[,n] , 1e-8)
    pOM_n =  sapply(1:length(mass), function(x)  ifelse(mass[x]>0, (cohorts.vol$OM[x,n])/mass[x], 0) )
    pOM_n = pmin(pOM_n,1)
    BD =  mapply( bulk.densityMass, pOM_n, (g))
    
  }
  return(list(BD=BD, pOM=pOM_n))
}

# porosity = function(massAbv, rate, minPor, maxPor, pOM, mid){
#  # sig=  1/(1+(exp((massAbv*(1-pOM^0.5)-mid)*rate)))
#  sig=  1/(1+(exp((massAbv-mid)*rate)))
#   sigMax = 1/(1+(exp((0-mid)*rate)))
#   sig/sigMax*(maxPor-minPor)+minPor
# }

porosity = function(density, rate, minPor, maxPor, pOM, mid){
  # sig=  1/(1+(exp((massAbv*(1-pOM^0.5)-mid)*rate)))
  sig=  1/(1+(exp((density-mid)*rate)))
  sigMax = 1/(1+(exp((0-mid)*rate)))
  sig/sigMax*(maxPor-minPor)+minPor
}

#den=seq(0,1,0.05)
#plot(den, sapply(den, porosity, rate=5, minPor=0.6, maxPor=0.8, pOM=0.5, mid=0))
#porosityRate=5

#uses rowsums across species
addVOLs=function(Volumesn, doRowSum){
  if(doRowSum){
    Volumesn$vROM+ Volumesn$vMIN + Volumesn$vLOM  + Volumesn$vSAND +
      (Volumesn$vCRP)+ rowSums(Volumesn$vCOARSEROOTS) + rowSums(Volumesn$vSMALLROOTS) + rowSums(Volumesn$vFINEROOTS)+ rowSums(Volumesn$vSTRUCROOTS)+
      (Volumesn$vFRP)+(Volumesn$vSRP)+(Volumesn$vWOODLITTER) + Volumesn$Vvoid
  } else {
    Volumesn$vROM+ Volumesn$vMIN + Volumesn$vLOM  + Volumesn$vSAND +
      (Volumesn$vCRP)+ (Volumesn$vCOARSEROOTS) + (Volumesn$vSMALLROOTS) + (Volumesn$vFINEROOTS)+ (Volumesn$vSTRUCROOTS)+
      (Volumesn$vFRP)+(Volumesn$vSRP)+(Volumesn$vWOODLITTER) + Volumesn$Vvoid
  }
}

#add volume of solid material, not voids
addVOLSTUFF=function(Volumesn, doRowSum){
  if(doRowSum){
    Volumesn$vROM+ Volumesn$vMIN + Volumesn$vLOM  + Volumesn$vSAND +
      (Volumesn$vCRP)+ rowSums(Volumesn$vCOARSEROOTS) + rowSums(Volumesn$vSMALLROOTS) +
      rowSums(Volumesn$vFINEROOTS)+ rowSums(Volumesn$vSTRUCROOTS)+
      (Volumesn$vFRP)+(Volumesn$vSRP)+(Volumesn$vWOODLITTER)
  } else {
    Volumesn$vROM+ Volumesn$vMIN + Volumesn$vLOM  + Volumesn$vSAND +
      (Volumesn$vCRP)+ sum(Volumesn$vCOARSEROOTS) + sum(Volumesn$vSMALLROOTS) +
      sum(Volumesn$vFINEROOTS)+ sum(Volumesn$vSTRUCROOTS)+
      (Volumesn$vFRP)+(Volumesn$vSRP)+(Volumesn$vWOODLITTER)
  }
}


#add root vol - EXCEPT structural roots -- they are always added to total volume and will displace pore space
addVOLROOTS = function(Volumesn, doRowSums){
  if(doRowSums){
    (Volumesn$vCRP)+ rowSums(Volumesn$vCOARSEROOTS) + rowSums(Volumesn$vSMALLROOTS) + rowSums(Volumesn$vFINEROOTS)+ #rowSums(Volumesn$vSTRUCROOTS)+rowSums(Volumesn$vSRP)+rowSums(Volumesn$vWOODLITTER)
      (Volumesn$vFRP)
  } else {
    (Volumesn$vCRP)+ (Volumesn$vCOARSEROOTS) + (Volumesn$vSMALLROOTS) + (Volumesn$vFINEROOTS)+ #(Volumesn$vSTRUCROOTS)+(Volumesn$vSRP)+(Volumesn$vWOODLITTER)
      (Volumesn$vFRP)
  }
}


#during spinup, compact bottom layer too

updateVOLs = function(cohorts, oldCohorts,porosityRate,maxPorosity, minPorosity, maxDeltaPor,midPt,DT, Initialize, SPINUP){
  #cohorts=roundCohorts(cohorts)  #round off small masses
  Act=1:nrow(cohorts$LOM)
  
  cohorts.vol = list(OM=cohorts$ROM[Act,]+cohorts$LOM[Act,], MIN=cohorts$MIN[Act,] + cohorts$SAND[Act,])  #for calculate of %OM
  
  #  print(dim(cohorts.vol$OM))
  # massAbv= cbind(0,t(apply(FUN=cumsum, cohorts$ROM[Act,]  +cohorts$LOM[Act,]+cohorts$MIN[Act,] + cohorts$SAND[Act,], MARGIN=1  )))[,1:ncohorts]
  #  massAbv=round(massAbv,8)
  depth= cbind(0,t(apply(FUN=cumsum, cohorts$VOL[Act,], MARGIN=1  )))[,1:ncohorts]
  mass = cbind(0,t(apply(FUN=cumsum, cohorts.vol$OM+cohorts.vol$MIN, MARGIN=1)))[,1:ncohorts]
  density=mass/depth
  density = ifelse(is.nan(density) | is.infinite(density),0, density)  
  #  massAbv=round(massAbv,8)
  
  #  massAbv1=cbind(0,t(apply(FUN=cumsum, oldCohorts$ROM[Act,]  +oldCohorts$LOM[Act,]+oldCohorts$MIN[Act,] + oldCohorts$SAND[Act,], MARGIN=1  )))[,1:ncohorts]
  
  #  massAbv= cbind(0,t(apply(FUN=cumsum, cohorts$ROM[Act,]  +cohorts$LOM[Act,]+cohorts$MIN[Act,] + cohorts$SAND[Act,] + ( apply(cohorts$FINEROOTS[Act,,]+cohorts$SMALLROOTS[Act,,]+cohorts$COARSEROOTS[Act,,]+cohorts$CRP[Act,,]+cohorts$FRP[Act,,]+ cohorts$SRP[Act,,]+cohorts$STRUCROOTS[Act,,], FUN=sum, MARGIN=c(1,3)) ) ,  MARGIN=1  )))[,1:ncohorts]
  #  massAbv1= cbind(0,t(apply(FUN=cumsum, oldCohorts$ROM[Act,]  +oldCohorts$LOM[Act,]+oldCohorts$MIN[Act,] + oldCohorts$SAND[Act,] + ( apply(oldCohorts$FINEROOTS[Act,,]+oldCohorts$SMALLROOTS[Act,,]+oldCohorts$COARSEROOTS[Act,,]+oldCohorts$CRP[Act,,]+oldCohorts$FRP[Act,,]+ oldCohorts$SRP[Act,,]+oldCohorts$STRUCROOTS[Act,,], FUN=sum, MARGIN=c(1,3)) ) ,  MARGIN=1  )))[,1:ncohorts]
  
  # recalPor= ifelse(rowSums(massAbv-massAbv1)==0, 10,1) #update porosity only if the massAbv has changed
  
  pOM = (cohorts.vol$OM/(cohorts.vol$OM+cohorts.vol$MIN))
  porold = cohorts$POR[Act,]
  
  pornew = array(mapply(FUN=porosity, density=density, MoreArgs=list(rate=porosityRate, minPor=minPorosity, maxPor=maxPorosity, mid=midPt), pOM=pOM), dim=c(length(Act), ncohorts))
  # print(pornew)
  pornew=round(pornew,8)
  porold=round(porold, 8)
  
  # print(pornew-porold)
  if(Initialize){  #for the first run, allow porosity to increase; assumes that cohorts are initially overfilled
    # cohorts$POR = array(pmin(porold,pornew), dim=c(length(elev), ncohorts)) #ratchet; can only compact/decrease porosity
    cohorts$POR[Act,] = pornew #array(pmin(porold, pornew), dim=c(length(elev), ncohorts)) #ratchet; can only compact/decrease porosity
    
  } else {
    if(SPINUP){
      deltaPor = porold[,1:(ncohorts)] - pmin(porold[,1:(ncohorts)],pornew[,1:(ncohorts)])
      #    cohorts$POR[,1:(ncohorts-1)] =  array(pmin(porold[,1:(ncohorts-1)],pornew[,1:(ncohorts-1)]), dim=c(length(elev), ncohorts-1)) #ratchet; can only compact/decrease porosity
      deltaPorSign=sign(deltaPor)
      deltaPor=abs(deltaPor)
      #rateArray= array(rep(maxDeltaPor*DT/cohorts$VOL[Act,1:(ncohorts-1)], length(Act)), dim=c(length(Act),ncohorts-1 ))
      rateArray= array(maxDeltaPor*DT/cohorts$VOL[Act,1:(ncohorts)], dim=c(length(Act),ncohorts ))
      porDif= pmin(deltaPor,rateArray)*deltaPorSign
      porDif= pmax(0, porDif)  #can only compact
      cohorts$POR[Act,1:(ncohorts)] = porold[,1:(ncohorts)] - porDif  #pmin(deltaPor,array(rep(maxDeltaPor*dt/cohorts$VOL[Act,1:(ncohorts-1)], length(Act)), dim=c(length(Act),ncohorts-1 )))*deltaPorSign
      
    }else {
      deltaPor = porold[,1:(ncohorts-1)] - pmin(porold[,1:(ncohorts-1)],pornew[,1:(ncohorts-1)])
      #    cohorts$POR[,1:(ncohorts-1)] =  array(pmin(porold[,1:(ncohorts-1)],pornew[,1:(ncohorts-1)]), dim=c(length(elev), ncohorts-1)) #ratchet; can only compact/decrease porosity
      deltaPorSign=sign(deltaPor)
      deltaPor=abs(deltaPor)
      #rateArray= array(rep(maxDeltaPor*DT/cohorts$VOL[Act,1:(ncohorts-1)], length(Act)), dim=c(length(Act),ncohorts-1 ))
      rateArray= array(maxDeltaPor*DT/cohorts$VOL[Act,1:(ncohorts-1)], dim=c(length(Act),ncohorts-1 ))
      porDif= pmin(deltaPor,rateArray)*deltaPorSign
      porDif= pmax(0, porDif)  #can only compact
      cohorts$POR[Act,1:(ncohorts-1)] = porold[,1:(ncohorts-1)] - porDif  #pmin(deltaPor,array(rep(maxDeltaPor*dt/cohorts$VOL[Act,1:(ncohorts-1)], length(Act)), dim=c(length(Act),ncohorts-1 )))*deltaPorSign
    }
  }
  
  
  #  if(Initialize){
  vstuff= addVOLs(list(vROM = cohorts$ROM/OMden, vLOM=cohorts$LOM/OMden,  vMIN=cohorts$MIN/MINden, vSAND=cohorts$SAND/MINden,
                       vFRP=          array(rowSums(sapply(1:nSpecies, function(x) (cohorts$FRP[,x,]/(root.den[x]/(1-root_por[x]))))), dim=c(length(elev), ncohorts)),
                       vCRP=          array(rowSums(sapply(1:nSpecies, function(x) (cohorts$CRP[,x,]/(root.den[x]/(1-root_por[x]))))), dim=c(length(elev), ncohorts)),
                       vSRP =         array(rowSums(sapply(1:nSpecies, function(x) (cohorts$SRP[,x,])/(wood.den[x]/(1-root_por[x])))), dim=c(length(elev), ncohorts)),
                       vFINEROOTS =    array(rowSums(sapply(1:nSpecies, function(x) (cohorts$FINEROOTS[,x,]/root.den[x])) ),  dim=c(length(elev), ncohorts)), #add together the root volumes for each species. They will get moved down together 
                       vSMALLROOTS =   array(rowSums(sapply(1:nSpecies, function(x) (cohorts$SMALLROOTS[,x,]/root.den[x])) ), dim=c(length(elev), ncohorts)),
                       vCOARSEROOTS =  array(rowSums(sapply(1:nSpecies, function(x) (cohorts$COARSEROOTS[,x,]/root.den[x])) ), dim=c(length(elev), ncohorts)),
                       vSTRUCROOTS =   array(rowSums(sapply(1:nSpecies, function(x) (cohorts$STRUCROOTS[,x,]/wood.den[x])) ), dim=c(length(elev), ncohorts)),
                       vWOODLITTER=    array(rowSums(sapply(1:nSpecies, function(x) (cohorts$WOODLITTER[,x,]/wood.den[x])) ), dim=c(length(elev), ncohorts)),
                       Vvoid = array(0, dim=c(length(elev), ncohorts))),
                  F)
  
  vRoots = addVOLROOTS(list(vFRP=          array(rowSums(sapply(1:nSpecies, function(x) (cohorts$FRP[,x,]/(root.den[x]/(1-root_por[x]))))), dim=c(length(elev), ncohorts)),
                            vCRP=          array( rowSums(sapply(1:nSpecies, function(x) (cohorts$CRP[,x,]/(root.den[x]/(1-root_por[x]))))), dim=c(length(elev), ncohorts)),
                            # vSRP =         array(rowSums(sapply(1:nSpecies, function(x) (cohorts$SRP[,x,])/(wood.den[x]/(1-root_por[x])))), dim=c(length(elev), ncohorts)),
                            vFINEROOTS =    array(rowSums(sapply(1:nSpecies, function(x) (cohorts$FINEROOTS[,x,]/root.den[x])) ),  dim=c(length(elev), ncohorts)), #add together the root volumes for each species. They will get moved down together 
                            vSMALLROOTS =   array(rowSums(sapply(1:nSpecies, function(x) (cohorts$SMALLROOTS[,x,]/root.den[x])) ), dim=c(length(elev), ncohorts)),
                            vCOARSEROOTS =  array(rowSums(sapply(1:nSpecies, function(x) (cohorts$COARSEROOTS[,x,]/root.den[x])) ), dim=c(length(elev), ncohorts)) ),
                       #    vSTRUCROOTS =   array(rowSums(sapply(1:nSpecies, function(x) (cohorts$STRUCROOTS[,x,]/wood.den[x])) ), dim=c(length(elev), ncohorts)),
                       #    vWOODLITTER=    array(rowSums(sapply(1:nSpecies, function(x) (cohorts$WOODLITTER[,x,]/wood.den[x])) ), dim=c(length(elev), ncohorts)) )
                       F)
  vRootsOLD = addVOLROOTS(list(vFRP=          array(rowSums(sapply(1:nSpecies, function(x) (oldCohorts$FRP[,x,]/(root.den[x]/(1-root_por[x]))))), dim=c(length(elev), ncohorts)),
                               vCRP=          array( rowSums(sapply(1:nSpecies, function(x) (oldCohorts$CRP[,x,]/(root.den[x]/(1-root_por[x]))))), dim=c(length(elev), ncohorts)),
                               # vSRP =         array(rowSums(sapply(1:nSpecies, function(x) (cohorts$SRP[,x,])/(wood.den[x]/(1-root_por[x])))), dim=c(length(elev), ncohorts)),
                               vFINEROOTS =    array(rowSums(sapply(1:nSpecies, function(x) (oldCohorts$FINEROOTS[,x,]/root.den[x])) ),  dim=c(length(elev), ncohorts)), #add together the root volumes for each species. They will get moved down together 
                               vSMALLROOTS =   array(rowSums(sapply(1:nSpecies, function(x) (oldCohorts$SMALLROOTS[,x,]/root.den[x])) ), dim=c(length(elev), ncohorts)),
                               vCOARSEROOTS =  array(rowSums(sapply(1:nSpecies, function(x) (oldCohorts$COARSEROOTS[,x,]/root.den[x])) ), dim=c(length(elev), ncohorts)) ),
                          #    vSTRUCROOTS =   array(rowSums(sapply(1:nSpecies, function(x) (cohorts$STRUCROOTS[,x,]/wood.den[x])) ), dim=c(length(elev), ncohorts)),
                          #    vWOODLITTER=    array(rowSums(sapply(1:nSpecies, function(x) (cohorts$WOODLITTER[,x,]/wood.den[x])) ), dim=c(length(elev), ncohorts)) )
                          F)
  
  # vRootsOLD = addVOLROOTS(list(vFRP=       oldCohorts$FRP[Act,]/(root.den/(1-root_por)),
  #                           vCRP=          oldCohorts$CRP[Act,]/(root.den/(1-root_por)),
  #                           # vSRP =         array(rowSums(sapply(1:nSpecies, function(x) (cohorts$SRP[,x,])/(wood.den[x]/(1-root_por[x])))), dim=c(length(elev), ncohorts)),
  #                           vFINEROOTS =   apply(oldCohorts$FINEROOTS[Act,,], FUN=sum, MARGIN=c(1,3))/root.den, #add together the root volumes for each species. They will get moved down together 
  #                           vSMALLROOTS =  apply(oldCohorts$SMALLROOTS[Act,,], FUN=sum, MARGIN=c(1,3))/root.den,
  #                           vCOARSEROOTS = apply(oldCohorts$COARSEROOTS[Act,,], FUN=sum, MARGIN=c(1,3))/root.den),
  #                         #    vSTRUCROOTS =   array(rowSums(sapply(1:nSpecies, function(x) (cohorts$STRUCROOTS[,x,]/wood.den[x])) ), dim=c(length(elev), ncohorts)),
  #                      #    vWOODLITTER=    array(rowSums(sapply(1:nSpecies, function(x) (cohorts$WOODLITTER[,x,]/wood.den[x])) ), dim=c(length(elev), ncohorts)) )
  #                     F)
  
  deltavRoots= vRoots-vRootsOLD                             #change in root volume
  minVoid= vstuff*minPorosity/(1-minPorosity)               #minimum void space to maintain minimum porosity
  # t(array(z.depth, dim=c(ncohorts,length(elev)) ))
  #Vvoid=vstuff*cohorts$POR/(1-cohorts$POR)
  VvoidLessRoots = (vstuff-deltavRoots)*cohorts$POR[Act,]/(1-cohorts$POR[Act,])  #void space of particulates
  cohorts$Vvoid[Act,] = pmax(minVoid, VvoidLessRoots)  #-vRoots            #new void space. roots can reduce pore space up to the minimum void space
  cohorts$VOL[Act,] = round(vstuff+ cohorts$Vvoid[Act,],8)                                #new volume = vol solids + vol voids
  cohorts$POR =  round(pmin(cohorts$POR, cohorts$Vvoid/cohorts$VOL),8)        #new porosity.  can only decrease porosity
  #   cohorts$POR =  round(cohorts$Vvoid/cohorts$VOL,8)  #always update porosity. can increase
  cohorts$POR = ifelse(is.nan(cohorts$POR),maxPorosity, cohorts$POR)
  
  # } else { 
  #   
  #   vstuff=cohorts$VOL-cohorts$Vvoid #volume of solids
  #   cohorts$Vvoid=vstuff*cohorts$POR/(1-cohorts$POR) #update void space based on new porosity
  #   cohorts$VOL=vstuff+cohorts$Vvoid #recalculate volume
  # }
  
  #cohorts$Vvoid/cohorts$VOL - cohorts$POR #check porosity
  
  bdCal =calc.BD(cohorts.vol, 1)
  # #set new BD, based on root growth, mineral deposition, decomposition, etc
  cohorts$BD[Act,1] =  bdCal$BD     
  cohorts$VOL[Act,1] = (cohorts.vol$OM[,1] + cohorts.vol$MIN[,1]) / cohorts$BD[Act,1] + rowSums(t(t(cohorts$WOODLITTER[Act,,1])/wood.den)) + rowSums(t(t(cohorts$STRUCROOTS[Act,,1])/wood.den))
  cohorts$POR[Act,1] = (cohorts$VOL[Act,1]-(cohorts.vol$OM[,1]/OMden) - cohorts.vol$MIN[,1]/(MINden) -rowSums(t(t(cohorts$WOODLITTER[Act,,1])/wood.den)) -rowSums(t(t(cohorts$STRUCROOTS[Act,,1])/wood.den))  )/cohorts$VOL[Act,1]
  cohorts$Vvoid[Act,1] = cohorts$VOL[Act,1]*cohorts$POR[Act,1] 
  
  # cohorts$Vvoid = ifelse(is.nan(cohorts$Vvoid),  t(array(rep(z.depth, nrow(cohorts$Vvoid)), dim=c(length(z.depth),nrow(cohorts$Vvoid)))), cohorts$Vvoid)
  cohorts$Vvoid = ifelse(is.nan(cohorts$Vvoid), 0, cohorts$Vvoid)  
  
  return(cohorts)
}


which_layer<-function(z_in,cz){
  min(which(z_in<=cz))
}

group_layer<-function(stuff,bins,x){
  sum(stuff[bins==x], na.rm=T)
}

order_and_apply <- function(Z_row, bin_to_row) {
  oZ <- order(Z_row)        # Get order of indices for Z_row
  bins <- bin_to_row[oZ]    # Apply this order to bin_to_row
  return(bins)
}

process_row1 <- function(Z_row, oZ_row, STU_row) {
  sorted_Z <- Z_row[oZ_row]
  fr <- diff(c(0, sorted_Z)) * STU_row[oZ_row]
  return(fr)
}

updateMatrial=function(stuff, dz, bins, bin_from, Z, oZ){
  
  STU=cbind(stuff,t(sapply(1:nrow(stuff), function(x) {stuff[x,][bin_from[x,]] } )) )
  
  fr <- t(sapply(1:nrow(Z), function(x) {
    process_row1( Z[x,], oZ[x,], STU[x,])
  }))
  fr=ifelse(is.na(fr),0,fr)
  
  new_stuff= t(sapply(1:nrow(bins), function(x){
    unique_bins=unique(bins[x,])
    sapply(unique_bins, function(i) group_layer(fr[x,], bins[x,], i))/dz[x,]
  } ))
  new_stuff
}

updateDepthPOR = function(cohorts, oldCohorts,wood.den, root.den, root_por,porosityRate,maxPorosity, minPorosity, maxDeltaPor, midPt,DT,Initialize, SPINUP ){
Act = 1:nrow(cohorts$VOL)
  #initialize volume
if(Initialize==T){
  print('initialize cohorts')
}
  cohorts= updateVOLs(cohorts,oldCohorts,porosityRate,maxPorosity, minPorosity, maxDeltaPor, midPt,DT, Initialize, SPINUP)
  
  if(SPINUP){
    print('spinup up cohorts')
    reps=ncohorts
  } else {
    reps=1
  }
  
  for(i in 1:reps){
   # print(i)
    z.depth.l=length(z.depth)
    dz= t(array( rep(z.depth[1:(z.depth.l-1)], nrow(cohorts$VOL)), dim=c(z.depth.l-1,nrow(cohorts$VOL))))  # c(1.2,2,2,2,5,5,5,8,9,10) #vols minus bottom layer z.depth
    ndz= cohorts$VOL#[,1:(ncohorts-1)]  # c(1.2,3.22,9,15,44,20) #vols now - after material is added, roots grow, etc
    #stuff=cbind(ndz*0.2) #fictional material (ROM, LOM)
    
    nz=ncol(dz)
    
    dz=cbind(dz,rowSums(ndz)-rowSums(dz)) #make bottom layer thickness  ***** make catch for negative values
    
    cndz = t(apply(ndz, FUN=cumsum, MARGIN=1))
    cdz = t(apply(dz, FUN=cumsum, MARGIN=1))
    
    bin_from <- t(mapply(function(cdz_row, cndz_row) {
      sapply(cdz_row, function(z_in) which_layer(z_in, cndz_row))
    }, split(cdz, row(cdz)), split(cndz, row(cndz))))
    
    bin_to <- t(mapply(function(cdz_row, cndz_row) {
      sapply(cdz_row, function(z_in) which_layer(z_in, cndz_row))
    }, split(cbind(cndz,cdz), row(cbind(cndz,cdz))), split(cdz, row(cdz))))
    
    Z=cbind(cndz,cdz)
    oZ = t(apply(Z, 1,order))
    
    bins <- t(mapply(function(Z_row, bin_to_row) {
      order_and_apply(Z_row, bin_to_row)
    }, split(Z, row(Z)), split(bin_to, row(bin_to))))
    
    cohorts$LOM = updateMatrial(cohorts$LOM/cohorts$VOL, dz, bins, bin_from, Z, oZ)*dz
    cohorts$ROM = updateMatrial(cohorts$ROM/cohorts$VOL, dz, bins, bin_from, Z, oZ)*dz
    cohorts$MIN = updateMatrial(cohorts$MIN/cohorts$VOL, dz, bins, bin_from, Z, oZ)*dz
    cohorts$SAND = updateMatrial(cohorts$SAND/cohorts$VOL, dz, bins, bin_from, Z, oZ)*dz
    cohorts$Vvoid = updateMatrial(cohorts$Vvoid/cohorts$VOL, dz, bins, bin_from, Z, oZ)*dz
 
   
    for(sp in 1:nSpecies){
      cohorts$WOODLITTER[,sp,]= updateMatrial(cohorts$WOODLITTER[,sp,]/cohorts$VOL, dz, bins, bin_from, Z, oZ)*dz
      cohorts$FRP[,sp,] = updateMatrial(cohorts$FRP[,sp,]/cohorts$VOL, dz, bins, bin_from, Z, oZ)*dz
      cohorts$CRP[,sp,] = updateMatrial(cohorts$CRP[,sp,]/cohorts$VOL, dz, bins, bin_from, Z, oZ)*dz
      cohorts$SRP[,sp,] = updateMatrial(cohorts$SRP[,sp,]/cohorts$VOL, dz, bins, bin_from, Z, oZ)*dz
      cohorts$FINEROOTS[,sp,] =   updateMatrial(cohorts$FINEROOTS[,sp,]/cohorts$VOL, dz, bins, bin_from, Z, oZ)*dz
      cohorts$COARSEROOTS[,sp,] = updateMatrial(cohorts$COARSEROOTS[,sp,]/cohorts$VOL, dz, bins, bin_from, Z, oZ)*dz
      cohorts$STRUCROOTS[,sp,] =  updateMatrial(cohorts$STRUCROOTS[,sp,]/cohorts$VOL, dz, bins, bin_from, Z, oZ)*dz
      cohorts$SMALLROOTS[,sp,] =  updateMatrial(cohorts$SMALLROOTS[,sp,]/cohorts$VOL, dz, bins, bin_from, Z, oZ)*dz
    }
    
    newVols = list(vROM = cohorts$ROM/OMden, vLOM=cohorts$LOM/OMden,  vMIN=cohorts$MIN/MINden, vSAND=cohorts$SAND/MINden,
                    vFRP=          (sapply(1:nSpecies, function(x) (cohorts$FRP[,x,]/(root.den[x]/(1-root_por[x]))))),
                    vCRP=          (sapply(1:nSpecies, function(x) (cohorts$CRP[,x,]/(root.den[x]/(1-root_por[x]))))),
                    vSRP =         (sapply(1:nSpecies, function(x) (cohorts$SRP[,x,])/(wood.den[x]/(1-root_por[x])))),
                    vFINEROOTS =   (sapply(1:nSpecies, function(x) (cohorts$FINEROOTS[,x,]/root.den[x])) ),  #add together the root volumes for each species. They will get moved down together 
                    vSMALLROOTS =  (sapply(1:nSpecies, function(x) (cohorts$SMALLROOTS[,x,]/root.den[x])) ),
                    vCOARSEROOTS = (sapply(1:nSpecies, function(x) (cohorts$COARSEROOTS[,x,]/root.den[x])) ),
                    vSTRUCROOTS =  (sapply(1:nSpecies, function(x) (cohorts$STRUCROOTS[,x,]/wood.den[x])) ),
                    vWOODLITTER=   (sapply(1:nSpecies, function(x) (cohorts$WOODLITTER[,x,]/wood.den[x])) ),
                    Vvoid =         cohorts$Vvoid)
    
    cohorts$POR= newVols$Vvoid/dz
    cohorts$VOL=dz

    vRootsOLD = addVOLROOTS(list(vFRP=          array(rowSums(sapply(1:nSpecies, function(x) (cohorts$FRP[,x,]/(root.den[x]/(1-root_por[x]))))), dim=c(length(elev), ncohorts)),
                              vCRP=          array( rowSums(sapply(1:nSpecies, function(x) (cohorts$CRP[,x,]/(root.den[x]/(1-root_por[x]))))), dim=c(length(elev), ncohorts)),
                              # vSRP =         array(rowSums(sapply(1:nSpecies, function(x) (cohorts$SRP[,x,])/(wood.den[x]/(1-root_por[x])))), dim=c(length(elev), ncohorts)),
                              vFINEROOTS =    array(rowSums(sapply(1:nSpecies, function(x) (cohorts$FINEROOTS[,x,]/root.den[x])) ),  dim=c(length(elev), ncohorts)), #add together the root volumes for each species. They will get moved down together 
                              vSMALLROOTS =   array(rowSums(sapply(1:nSpecies, function(x) (cohorts$SMALLROOTS[,x,]/root.den[x])) ), dim=c(length(elev), ncohorts)),
                              vCOARSEROOTS =  array(rowSums(sapply(1:nSpecies, function(x) (cohorts$COARSEROOTS[,x,]/root.den[x])) ), dim=c(length(elev), ncohorts)) ),
                         #    vSTRUCROOTS =   array(rowSums(sapply(1:nSpecies, function(x) (cohorts$STRUCROOTS[,x,]/wood.den[x])) ), dim=c(length(elev), ncohorts)),
                         #    vWOODLITTER=    array(rowSums(sapply(1:nSpecies, function(x) (cohorts$WOODLITTER[,x,]/wood.den[x])) ), dim=c(length(elev), ncohorts)) )
                         F)
    
    #initialize volume
    if(SPINUP){
      cohorts= updateVOLs(cohorts,oldCohorts,porosityRate,maxPorosity, minPorosity, maxDeltaPor, midPt,DT,F, SPINUP)
    }
    
  }
  
  
  cohorts.vol = list(OM=cohorts$ROM[Act,]+cohorts$LOM[Act,], MIN=cohorts$MIN[Act,] + cohorts$SAND[Act,],
                     liveRoots = cohorts$FINEROOTS[Act,,]+cohorts$SMALLROOTS[Act,,] + cohorts$COARSEROOTS[Act,,],
                     deadRoots = cohorts$FRP[Act,,] + cohorts$CRP[Act,,],
                     woody = cohorts$STRUCROOTS[Act,,] + cohorts$WOODLITTER[Act,,],
                     deadWoody = cohorts$SRP[Act,,])  #for calculate of %OM
  
   if(nSpecies>1){
    liveRootVol = array(rowSums(sapply(1:nSpecies, function(x) (cohorts.vol$liveRoots[,x,]/root.den[x]))), dim=c(length(elev), ncohorts))
    deadRootVol = array(rowSums(sapply(1:nSpecies, function(x) (cohorts.vol$deadRoots[,x,]/(root.den[x]/(1-root_por[x]))))), dim=c(length(elev), ncohorts))
    woodyVol =  array(rowSums(sapply(1:nSpecies, function(x) (cohorts.vol$woody[,x,]/(wood.den[x])))), dim=c(length(elev), ncohorts))
    deadWoodyVol =  array(rowSums(sapply(1:nSpecies, function(x) (cohorts.vol$deadWoody[,x,]/(wood.den[x]/(1-root_por[x]))))), dim=c(length(elev), ncohorts))
    cohorts$BD[Act,]  = (cohorts.vol$OM+cohorts.vol$MIN+ apply(cohorts.vol$liveRoots+cohorts.vol$deadRoots+cohorts.vol$woody+cohorts.vol$deadWoody, FUN=sum, MARGIN=c(1,3) ))/((liveRootVol+deadRootVol+woodyVol+deadWoodyVol+cohorts.vol$OM/OMden+ cohorts.vol$MIN/MINden+ cohorts$POR[Act,]*(cohorts.vol$OM/OMden+ cohorts.vol$MIN/MINden+liveRootVol+deadRootVol+woodyVol+deadWoodyVol)/(1- cohorts$POR[Act,])))
   }else {
     liveRootVol = array(cohorts.vol$liveRoots/root.den, dim=c(length(elev), ncohorts))
     deadRootVol = array(cohorts.vol$deadRoots/(root.den/(1-root_por)), dim=c(length(elev), ncohorts))
     woodyVol =  array(cohorts.vol$woody/wood.den, dim=c(length(elev), ncohorts))
     deadWoodyVol =  array(cohorts.vol$deadWoody/(wood.den/(1-root_por)), dim=c(length(elev), ncohorts)) 
     cohorts$BD[Act,]  = (cohorts.vol$OM+cohorts.vol$MIN+ cohorts.vol$liveRoots+cohorts.vol$deadRoots+cohorts.vol$woody+cohorts.vol$deadWoody)/((liveRootVol+deadRootVol+woodyVol+deadWoodyVol+cohorts.vol$OM/OMden+ cohorts.vol$MIN/MINden+ cohorts$POR[Act,]*(cohorts.vol$OM/OMden+ cohorts.vol$MIN/MINden+liveRootVol+deadRootVol+woodyVol+deadWoodyVol)/(1- cohorts$POR[Act,])))
   } 
  #cohorts$BD[Act,]  = (cohorts.vol$OM+cohorts.vol$MIN+ apply(cohorts.vol$liveRoots+cohorts.vol$deadRoots, FUN=sum, MARGIN=c(1,3) ))/((liveRootVol+deadRootVol+cohorts.vol$OM/OMden+ cohorts.vol$MIN/MINden+ cohorts$POR[Act,]*(cohorts.vol$OM/OMden+ cohorts.vol$MIN/MINden+liveRootVol+deadRootVol)/(1- cohorts$POR[Act,])))
  
  #With structural roots and wood litter
  cohorts$BD[Act,]  = (cohorts.vol$OM+cohorts.vol$MIN+ apply(cohorts.vol$liveRoots+cohorts.vol$deadRoots+cohorts.vol$woody+cohorts.vol$deadWoody, FUN=sum, MARGIN=c(1,3) ))/((liveRootVol+deadRootVol+woodyVol+deadWoodyVol+cohorts.vol$OM/OMden+ cohorts.vol$MIN/MINden+ cohorts$POR[Act,]*(cohorts.vol$OM/OMden+ cohorts.vol$MIN/MINden+liveRootVol+deadRootVol+woodyVol+deadWoodyVol)/(1- cohorts$POR[Act,])))
  
  cohorts
}

