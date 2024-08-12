#PLOT THINGS

source('src/MULTIPLOT_ggplot.R')

if(sd(SLR)>5){
  slvar='ENSO_Var'
  slrMax=10
  slrMin=-10
} else if(sd(SLR)>1){
  slvar='SL_Var'
  slrMax=4
  slrMin=-4
} else {
  slvar = 'NoSL_Var'
  slrMax=0.75
  slrMin=-.2
}

meanZ = mean(minCal$elevMSL)

getTargDepth = function(acc, yrs, depth) { 
   targDep = acc*yrs
   findCohort(depth, d=targDep)
  }



pyears = years-1
elev.df=data.frame()
yrTrack = NULL
for(i in 1:length(elev)){
 # targDepth=sapply(FUN=getTargDepth, mercc$acc[i,pyears], depth=mercc$cohorts.state$DEPTH[i,], yrs=pyears- spinupYrs )
#  omAccum = sum(mercc$cohorts.state$LOM[i,1:targDepth] + mercc$cohorts.state$ROM[i,1:targDepth]+  mercc$cohorts.state$FINEROOTS[i,,1:targDepth] +  mercc$cohorts.state$MAOM[i,1:targDepth])/pyears- spinupYrs
  
  
  elev.df = rbind(elev.df, cbind(mercc$elev_msl[i,1:pyears], mercc$acc[i,1:pyears], mercc$om[i,1:pyears], mercc$min[i,1:pyears], mercc$carbon[i,1:pyears], mercc$carbon_abv[i,1:pyears],  
                                 rep(1:pyears), rep(round(elev[i],1), pyears), SLR[2:years], c(diff(mercc$allOM[i,])[1:pyears-1],0) , c(diff(mercc$allMIN[i,])[1:pyears-1],0), c(diff(mercc$carbon[i,1:pyears]), 0)  ) )
 #  yrTrack=c(yrTrack, which.min(abs(mercc$elev_msl[i,500:pyears,]-meanZ ))+500)   #500    ### 
   yrTrack=c(yrTrack, which.min(abs(mercc$elev_msl[i,spinupYrs:pyears]-meanZ ))+spinupYrs) 
  }

colnames(elev.df) = c('Elevation_cm_MSL', 'Accretion_cm_yr',"OM_g_60cm", 'MIN_g_60cm', 'Carbon' , 'Carbon_Abv' ,'Years','init.z', 'SLR', 'OM_Accum', 'Min_Accum', 'C_Accum')
#elev.df = data.frame(z= c(mercc$elev_msl[1,,], mercc$elev_msl[2,,], mercc$elev_msl[3,,], mercc$elev_msl[4,,], mercc$elev_msl[5,,], mercc$elev_msl[6,,], mercc$elev_msl[7,,], mercc$elev_msl[8,,]), years=rep(1:years, length(elev)), init=c(rep(elev[1], years),  ) )
p_z = ggplot(elev.df) + geom_line(aes(x=Years, y=Elevation_cm_MSL, col=as.factor(elev.df$init.z), group=as.factor(elev.df$init.z)) ) +theme_classic() + geom_hline(aes(yintercept=mean(minCal$elevMSL), linetype='dotted')) +
  xlim(c(startX, endX)) + theme(legend.position='top') + guides(col=guide_legend(title='Initial elevation (cm MSL)')) + ylab('Elevation (cm, MSL)') +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +  scale_linetype_manual(name='', labels = c("Soil cores average"), values=c("dashed"))
#p_z

p_acc = ggplot(elev.df) + geom_line(aes(y=SLR, x=Years), linetype='dashed', col='darkgrey', alpha=0.75) + geom_hline(aes( yintercept =mean(minCal$acc), linetype='dotted')) + # geom_vline(xintercept = yrTrack, linetype='dashed', alpha=0.25)+
  geom_line(aes(x=Years, y=Accretion_cm_yr, col=as.factor(elev.df$init.z), group=as.factor(elev.df$init.z)) ) +theme_classic() +ylim(c(slrMin, slrMax)) +
  xlim(c(startX, endX)) + theme(legend.position='top') + ylab('Vertical change \n[Acc, SLR] (cm/yr)') +guides(col=guide_legend(title='Initial elevation (cm MSL)')) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +  scale_linetype_manual(name='', labels = c("Soil cores average"), values=c("dashed"))
#p_acc

p_om= ggplot(elev.df) + geom_line(aes(x=Years, y=OM_g_60cm, col=as.factor(elev.df$init.z), group=as.factor(elev.df$init.z)) ) +theme_classic() + 
   xlim(c(startX, endX)) + theme(legend.position='top') + guides(col=guide_legend(title='Initial Elevation (cm MSL)')) + ylab('Organics, top 60 cm (g)') + #geom_hline(aes(yintercept=mean(omarc.input$om_gcm2), linetype='dotted'))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +  scale_linetype_manual(name='', labels = c("Soil cores average"), values=c("dashed"))

p_min= ggplot(elev.df) + geom_line(aes(x=Years, y=MIN_g_60cm, col=as.factor(elev.df$init.z), group=as.factor(elev.df$init.z)) ) +theme_classic() + 
  xlim(c(startX, endX)) + theme(legend.position='top') + guides(col=guide_legend(title='Initial Elevation (cm MSL)')) + ylab('Mineral, top 60 cm (g)') + #geom_hline(aes(yintercept=mean(msarc.input$min_gcm2), linetype='dotted')) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))+  scale_linetype_manual(name='', labels = c("Soil cores average"), values=c("dashed"))

p_omAcc=ggplot(elev.df) + geom_line(aes(x=Years, y=OM_Accum, col=as.factor(elev.df$init.z), group=as.factor(elev.df$init.z)) ) +theme_classic() + 
  ylim(c(-.1, .1))+ xlim(c(startX, endX)) + theme(legend.position='top') + guides(col=guide_legend(title='Initial Elevation (cm MSL)')) + ylab('Organic accumumulation\n(g/cm2/yr)') + geom_hline(aes(yintercept=mean(minCal$omAcc), linetype='dotted')) +# geom_vline(xintercept = yrTrack, linetype='dashed', alpha=0.25) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))+  scale_linetype_manual(name='', labels = c("Soil cores average"), values=c("dashed"))

p_MinAcc=ggplot(elev.df) + geom_line(aes(x=Years, y=Min_Accum, col=as.factor(elev.df$init.z), group=as.factor(elev.df$init.z)) ) +theme_classic() + 
  xlim(c(startX, endX)) + theme(legend.position='top') + guides(col=guide_legend(title='Initial Elevation (cm MSL)')) + ylab('Mineral accumulation\n(g/cm2/yr)') + geom_hline(aes(yintercept=mean(minCal$minAcc), linetype='dotted')) + # geom_vline(xintercept = yrTrack, linetype='dashed', alpha=0.25)+
 theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))+  scale_linetype_manual(name='', labels = c("Soil cores average"), values=c("dashed"))

p_C=ggplot(elev.df) + geom_line(aes(x=Years, y= (elev.df$Carbon_Abv+ elev.df$Carbon)*1e8/1e6, col=as.factor(elev.df$init.z), group=as.factor(elev.df$init.z)) ) +theme_classic() + 
   xlim(c(startX, endX)) + theme(legend.position='top') + guides(col=guide_legend(title='Initial Elevation (cm MSL)')) + ylab('Carbon stock\n(Mg/ha)')  + #geom_hline(aes(yintercept=mean(msarc.input$core_sed), linetype='dotted')) + # geom_vline(xintercept = yrTrack, linetype='dashed', alpha=0.25)+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))+  scale_linetype_manual(name='', labels = c("Soil cores average"), values=c("dashed"))


# p_CAcc=ggplot(elev.df[elev.df$Years==1:(pyears),]) + geom_line(aes(x=Years, y= C_Accum*1e8/1e6, col=as.factor(init.z)) ) +theme_classic() + 
#   xlim(c(startX, endX)) + theme(legend.position='top') + guides(col=guide_legend(title='Initial Elevation (cm MSL)')) + ylab('Carbon sequestration\n(Mg/ha/yr)')  + #geom_hline(aes(yintercept=mean(msarc.input$core_sed), linetype='dotted')) + # geom_vline(xintercept = yrTrack, linetype='dashed', alpha=0.25)+
#   theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) #+  scale_linetype_manual(name='', labels = c("Soil cores average"), values=c("dashed"))

#mean(omarc.input$core_om) - elev.df$OM_Accum[(pyears-300):(pyears)]


p_omAccZ=ggplot(elev.df) + geom_line(aes(x=elev.df$Elevation_cm_MSL, y=OM_Accum, col=as.factor(elev.df$init.z), group=as.factor(elev.df$init.z)) ) +theme_classic() + 
  ylim(c(-.1, .1))+ xlim(c(startX, max(elev.df$Elevation_cm_MSL))) + theme(legend.position='top') + guides(col=guide_legend(title='Initial elevation (cm, MSL)')) + ylab('Organic Accumulation (g/cm2/yr)') + geom_hline(yintercept=mean(minCal$omAcc), linetype='dotted')+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))


cc= sitePars$CC

# coverdf = data.frame(cover = c(mercc$perCover[ncore,1,]/cc, mercc$perCover[ncore,2,]/cc, mercc$perCover[ncore,3,]/cc, mercc$perCover[ncore,4,]/cc, 
#                                mercc$perCover[ncore,5,]/cc, colSums(mercc$perCover[ncore,6:10,])/cc),
#                      Sp = c(rep(speciesList[1],years), rep(speciesList[2], years), rep('SOAL',years), rep('RHAP',years), rep('RHST', years), rep("DEAD",years)),
#                      years= rep(1:years, 6) )

coverdf = data.frame(cover=mercc$perCover[ncore,1,]/cc, Sp = rep(speciesList[1],years), years=1:years)
for(sp in 2:length(speciesList)){
  coverdf = rbind(coverdf, data.frame(cover=mercc$perCover[ncore,sp,]/cc, Sp = rep(speciesList[sp],years), years=1:years))
}
coverdf= rbind(coverdf, data.frame(cover=colSums(mercc$perCover[ncore, (length(speciesList)+1):(length(speciesList)*2),]/cc), Sp = rep('Dead', years), years=1:years))

# coverdf = data.frame(cover = c(mercc$perCover[ncore,1,,1]/cc, mercc$perCover[ncore,2,,1]/cc, colSums(mercc$perCover[ncore,3:4,,1])/cc),
#                      Sp = c(rep('BRGY',years), rep('BRGY2', years), rep("DEAD",years)), 
#                      years= rep(1:years, 3) )




pcov =ggplot(coverdf, aes(x=years, y=cover)) +geom_area(aes( fill=Sp), position='stack') + theme_classic() + xlab('Years')+
  ylab('Proportion')+ xlim(c(startX, endX)) + #ylim(c(0,1))+
  theme(legend.position='top')+ guides(fill=guide_legend(title=paste('Initial elevation: ',elev[ncore], ' (cm MSL)', sep=''), nrow=1))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
#pcov




# relVol = data.frame(lom=mercc$lomVol[ncore,1:pyears]/mercc$topDepth[ncore, 1:pyears],
#                     rom=mercc$romVOL[ncore,1:pyears]/mercc$topDepth[ncore, 1:pyears],
#                     min=mercc$minVol[ncore,1:pyears]/mercc$topDepth[ncore, 1:pyears], 
#                     liveRoots=mercc$liveRootVol[ncore,1:pyears]/mercc$topDepth[ncore, 1:pyears], 
#                     deadRoots=mercc$deadRootVol[ncore,1:pyears]/mercc$topDepth[ncore, 1:pyears])

relVol = data.frame(lom=mercc$lomVol[ncore,1:pyears],#/mercc$topDepth[ncore, 1:pyears],
                    rom=mercc$romVol[ncore,1:pyears],#/mercc$topDepth[ncore, 1:pyears],
                    min=mercc$minVol[ncore,1:pyears],#/mercc$topDepth[ncore, 1:pyears], 
                    liveRoots=mercc$liveRootVol[ncore,1:pyears],#/mercc$topDepth[ncore, 1:pyears], 
                    deadRoots=mercc$deadRootVol[ncore,1:pyears],
                    woodlitter = mercc$woodlitter[ncore,1:pyears])#/mercc$topDepth[ncore, 1:pyears])
               
relVol$total = rowSums(relVol)

relVol_melt = data.frame(percent = c(relVol$lom, relVol$rom, relVol$min, relVol$liveRoots, relVol$deadRoots, relVol$woodlitter)/relVol$total, 
                         material = c(rep('lom',nrow(relVol)), rep('rom',nrow(relVol)), rep('min', nrow(relVol)), rep('liveRoots', nrow(relVol)), rep('deadRoots', nrow(relVol)), rep('woodlitter', nrow(relVol))), year = rep(1:pyears, 6) )

p2 = ggplot(relVol_melt, aes(x=year, y=percent*100)) +geom_area(aes( fill=material), position='stack') + theme_classic() +  xlab('Years')+
  ylab('Percent') + xlim(c(startX, endX))+ theme(legend.position='top')+ guides(fill=guide_legend(title=paste('Elevation: ',elev[ncore], ' core vol. (60 cm)   Material', sep=''), nrow=1))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
#g cm-2 to 60 cm

if(plotCal==T & toFile){
  #scenName = paste('output_Pohnpei\\Figures\\Sensitivity\\', site.name,"_Yrs",startX, "_",endX,"_SLRmm", round(mean(SLR),2)*10,'_',slvar, "_Z", elev[ncore], '_CAL_7_1_20.png', sep='' )
  
  tiff(scenName, width=6.5*2.2, height=7*2, res=300, units='in', compress='lzw')
  multiplot(p_z, pcov, p2, p_acc, p_omAcc,p_MinAcc , cols=2) # p_C
  dev.off()
  
}else if (plotCal & !toFile){
  multiplot(p_z, pcov, p2, p_acc, p_omAcc, p_MinAcc, cols=2) #p_C

} else if(toFile) {
  png(scenName, width=350/2, height=275, res=300, units='mm')
  multiplot(p_z, pcov, p2)
  dev.off()
} else {
  multiplot(p_z, pcov, p2, cols=1)
}