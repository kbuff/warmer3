

  pOM = (mercc$cohorts.state$LOM+mercc$cohorts.state$ROM +apply(mercc$cohorts.state$FINEROOTS+mercc$cohorts.state$FRP+mercc$cohorts.state$COARSEROOTS+mercc$cohorts.state$CRP, sum, MARGIN=1) )/(mercc$cohorts.state$LOM+mercc$cohorts.state$ROM+mercc$cohorts.state$MIN+apply(mercc$cohorts.state$FINEROOTS+mercc$cohorts.state$FRP+mercc$cohorts.state$COARSEROOTS+mercc$cohorts.state$CRP, sum, MARGIN=1) )
  #pOM = (mercc$cohorts.state$LOM+mercc$cohorts.state$ROM )/(mercc$cohorts.state$LOM+mercc$cohorts.state$ROM+mercc$cohorts.state$MIN+ mercc$cohorts$SAND )
  
  pyears= years-1
  BD=OM=depth=z.initial=NULL
  for(i in 1:length(elev)){
    BD = c(BD, mercc$cohorts.state$BD[i,])
    OM = c(OM, pOM[i,]*100)
    depth = c(depth, c(0,cumsum(mercc$cohorts.state$VOL[i,])[1:(ncohorts-1)]))
    z.initial = c(z.initial, rep(round(elev[i],1), length(mercc$cohorts.state$VOL[i,])))
  }
  coreChar = data.frame(BD=BD, OM=OM, depth=depth, Elevation=z.initial)
  coreChar= subset(coreChar, coreChar$depth<=60)
  coreChar$Elevation=as.factor(coreChar$Elevation)
  
  cores$Type = 'Measured'
  
  theLines = data.frame(y=c(rep(-10, 60) ,rep(-20,60)), x=rep(1:60, 2), type = c(rep('dotted', 60), rep('dashed', 60) ) )
  
  plotOM = ggplot() +geom_point(data=cores, aes(x=Depth, y=OM, shape='16')) +
    geom_line(data=coreChar, aes(x=depth, y=OM, col=Elevation), show.legend = T) +
    theme_bw()+ ylab('Percent organic matter') + xlab('Depth (centimeters)')+ xlim(c(0,60)) + ylim(c(0,100))  +
    geom_line(data=theLines, aes(x=x, y=y, linetype=type))+
    geom_point(aes(x=-20, y=-20, shape='16'))+
    scale_shape_manual(name="", values=16, labels='Soil core data',
                       guide=guide_legend(override.aes = list(linetype=0 ), order=2)) +
    scale_linetype_manual(name='', values=c('dotted','dashed'),#labels = c("Published range", 'Published mean'),
                          guide=guide_legend(override.aes = list(linetype=c('dotted','dashed')),  order=1),
                          labels=c("Published range", 'Published mean'))
  
  pBD = ggplot() +geom_point(data=cores, aes(x=Depth, y=BD), show.legend = F)+ geom_line(data=coreChar, aes(x=depth, y=BD, col=Elevation))+# + geom_line(aes(x=depth[1,], y=pOM[1,]*100), col='red') +  geom_line(aes(x=depth[2,], y=pOM[2,]*100), col='darkgreen') +  geom_line(aes(x=depth[3,], y=pOM[3,]*100), col='blue') + xlim(c(0,60)) + ylim(c(0,100)) +
    theme_bw() + ylab('Bulk density (g/cm3)') + xlab('Depth (centimeters)')+ xlim(c(0,60)) + ylim(c(0,1))
  
  
  deadTotal=z.initial=litter=totRoots=liveRoots=yrs=totAGB=NULL
  for(i in 1:length(elev)){
    deadTotal = c(deadTotal,( mercc$deadRoots[i,1:pyears]) /
                    ( mercc$coarseRoots[i,1:pyears]+mercc$fineRoots[i,1:pyears]+mercc$smallRoots[i,1:pyears]+ mercc$deadRoots[i,1:pyears]))
    z.initial = c(z.initial, rep(round(elev[i],1), length(1:pyears)))
    yrs = c(yrs, 1:pyears)
    totRoots = c(totRoots,  mercc$fineRoots[i,1:pyears]+mercc$smallRoots[i,1:pyears]+ mercc$coarseRoots[i,1:pyears] +mercc$deadRoots[i,1:pyears])
    liveRoots = c(liveRoots, mercc$fineRoots[i,1:pyears]+mercc$smallRoots[i,1:pyears]+ mercc$coarseRoots[i,1:pyears])
    totAGB = c(totAGB, colSums(mercc$totABG[i,, 1:pyears]))
    litter =  c(litter, (mercc$leaflitter[i,1:pyears]+mercc$woodlitter[i,1:pyears])*100)
  }
  
  timeseriesOut = data.frame(Year=yrs, Elevation=z.initial, deadTotal=deadTotal,totRoots=totRoots, liveRoots=liveRoots, litter=litter )
  timeseriesOut$Elevation=as.factor(timeseriesOut$Elevation)
  
  pDeadLive=ggplot() +geom_line(data=timeseriesOut, aes(x=Year, y=deadTotal, col=Elevation), show.legend = T) +
    theme_bw()+ ylab('Dead:Total Root Biomass') + xlab('Years')+
    geom_hline(yintercept=calVals$deadTotalRootmin, linetype='dashed')+ geom_hline(yintercept=calVals$deadTotalRootmax, linetype='dashed')
  
  pTotRoots=ggplot() +geom_line(data=timeseriesOut, aes(x=Year, y=totRoots, col=Elevation), show.legend = T) +
    theme_bw()+ ylab('Total Roots (g/cm2)') + xlab('Years')+
    geom_hline(yintercept=calVals$totRootMassmin, linetype='dashed')  + geom_hline(yintercept=calVals$totRootMassmax, linetype='dashed')
  
  pLiveRoots=ggplot() +geom_line(data=timeseriesOut, aes(x=Year, y=liveRoots, col=Elevation), show.legend = T) +
    theme_bw()+ ylab('Live Roots (g/cm2)') + xlab('Years')+
    geom_hline(yintercept=calVals$liveRootsmin, linetype='dashed')+ geom_hline(yintercept=calVals$liveRootsmax, linetype='dashed')
  
  pTotABG=ggplot() +geom_line(data=timeseriesOut, aes(x=Year, y=totAGB*10000, col=Elevation), show.legend = T) +
    theme_bw()+ ylab('Total Aboveground Biomass (g/m2)') + xlab('Years')+
    geom_hline(yintercept=calVals$totAboveBiomassmin, linetype='dashed')  + geom_hline(yintercept=calVals$totAboveBiomassmax, linetype='dashed')
  
  
  pLitter=ggplot() +geom_line(data=timeseriesOut, aes(x=Year, y=litter, col=Elevation), show.legend = T) +
    theme_bw()+ ylab('Litterfall (Mg/ha/yr)') + xlab('Years')+
    geom_hline(yintercept=calVals$littermin, linetype='dashed') + geom_hline(yintercept=calVals$littermax, linetype='dashed')
  
  
  grid_arrange_shared_legend( plotOM, pTotABG, pLiveRoots, pBD,pDeadLive, pLitter ,  ncol=3,nrow=2)  #plotOM, pBD, pLitter
