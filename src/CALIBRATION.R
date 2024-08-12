
working.dir = 'C:\\MERCC\\WARMER_v3\\warmer_v3_github_20240812'# 'pathTo_WARMER_v3_folder'
setwd(working.dir)

library(ggplot2)
library(reshape2)
library(raster)
library(foreach)
library(doParallel)
library(snowfall)
library(doSNOW)
library(plyr)
library(tictoc)
require(pracma)
library(abind)

source('src/GeneralFunctions.R')
source('src/MainLoop.R')
source('src/Belowground.R')
source('src/updateDepthPOR.R')
source('src/spCover.R')
source('src/prepSpecies.R')
source('src/SSC_Deposition_Marani2010.R')


out.folder = '\\output'
site.name = 'Restoration'  #'Marsh'

print(site.name)
source(paste('input/', site.name, 'Parameters.R', sep='') )

#if 'user defined, the harmonic constituents must be loaded [HC] in the site parameter file. 
#Otherwise, 1 yr of tidal inundation is simulated based on the provided tidal range and the tide type
tide_type = 'user-defined'# 'mixed semi-diurnal'#  'semi-diurnal' #'diurnal'
  
source('src/GenTidalSignal.R')

## calculate sediment deposition based on tide signal with SSC=1 mg/L -- takes a bit, uses parallel processing
tic()
sedinundFun= genMineralFunction(theLevels)  #deposition fun; input= elevation cm rMSL, output: mineral deposition, assuming SSC=1.
toc()

#input is elevation (cm, MSL) and standing aboveground biomass. Multiply by the SSC [mg/L]
sedinundFun(20, 1000)*30

#source after inundFun and TR are defined 
source('src/spParameters.R')

ztemp = seq(min(theLevels$water)*100, max(theLevels$water)*100 , 1)
plot(ztemp, sapply(ztemp, sedinundFun, 1000)*30, type='l', xlab='Elevation (cm, MSL)', ylab='Sediment Deposition (g/cm2/yr)', main='SSC = 30 mg/L')
plot(ztemp, sapply(ztemp, inundFun)*100, type='l', xlab='Elevation (cm, MSL)', ylab='Flooding Duration (% time)')

##########################################################################################################
#Calibrate SSC against soil core mineral deposition. If no soil core data, set SSC based on best available data

ssc=mean(minCal$minAcc/sapply(minCal$elevMSL, sedinundFun, B1=1000))
plot(ztemp, sapply(ztemp, sedinundFun, 1000)*(ssc), type='l', xlab='Elevation (cm, MSL)', ylab='Sediment deposition g/cm2/yr', main=paste0('SSC= ', round(ssc,1)))
points(minCal$elevMSL, minCal$minAcc)

#Calculate the standard deviation of SSC using rates from multiple soil cores.
if(length(minCal$minAcc)>1){
ssc_temp=NULL
for(i in 1:length(minCal$minAcc)){
  ssc_temp=c(ssc_temp, minCal$minAcc[i]/sapply(minCal$elevMSL[i], sedinundFun, B1=1000))
}
ssc_sd=sd(ssc_temp)
} else {
  ssc_sd=5  
}



##########
#calibration porosity function to soil core BD and %OM
# 
# 
coreIDs = unique(cores$CoreID)
cores2=NULL
for(i in 1:length(coreIDs)){
  coreTemp = cores[cores$CoreID==coreIDs[i],]
  if(!is.na(sum(coreTemp$OM)) & !is.na(sum(coreTemp$BD))){
  int = diff(coreTemp$Depth)
  int=ifelse(int==3, 4, int)
  coreTemp$int = c(int[1], int)
  coreTemp$Mass = coreTemp$BD*coreTemp$int
  coreTemp$MassAbv = c( cumsum(coreTemp$Mass))
  coreTemp$MassAbv = c(0, coreTemp$MassAbv[1:(nrow(coreTemp)-1)])
  cores2=rbind(cores2, coreTemp)
  }
}
cores= cores2

surf=cores[cores$Depth==0,]

POR = 1- (surf$BD*(surf$OM/100/OMden)+(1-surf$OM/100)/MINden)/(surf$OM/100*OMden+(1-surf$OM/100)*MINden)
mean(POR)

calPor=function(pars, cores2){
  cores2=cores2[cores$Depth>3,]
  pOM=cores2$OM/100
  cores2$PORcalc = porosity(cores2$MassAbv, pars[1], pars[2], mean(POR), pOM=pOM, mid=0)
  cores2$BD2 = (cores2$Mass)/((cores2$Mass*cores2$OM/100/OMden+ cores2$Mass*(1-cores2$OM/100)/MINden+cores2$PORcalc*(cores2$Mass*cores2$OM/100/OMden+ cores2$Mass*(1-cores2$OM/100)/MINden)/(1-cores2$PORcalc)))
  cores2 = cores2[complete.cases(cores2$BD2),]
  mean(sqrt((cores2$BD-cores2$BD2)^2))
}

#calibration rate and minimum porosity. Max porosity a constant
out=optim(par=c(0.1, 0.6),fn=calPor, cores2=cores, method='L-BFGS-B',lower=c(0.001, 0.4, 0.4), upper=c(0.3, 0.8, 0.98))

porosityRate=out$par[1]
minPorosity=out$par[2]
maxPorosity=mean(POR)

pOM=cores$OM/100
cores$PORcalc = porosity(cores$MassAbv, porosityRate, minPorosity, maxPorosity, pOM, mid=0)
cores$BD2 = (cores$Mass)/((cores$Mass*cores$OM/100/OMden+ cores$Mass*(1-cores$OM/100)/MINden+cores$PORcalc*(cores$Mass*(cores$OM/100)/OMden+ cores$Mass*(1-cores$OM/100)/MINden)/(1-cores$PORcalc)))

plot(cores$BD, cores$BD2)
abline(a=0, b=1)

plot(cores$Depth, cores$BD, ylim=c(0,1.5))
points(cores$Depth, cores$BD2, col='red')





##############################################################################################################
#RUN CALIBRATION MODELS WITH SPINUP
########################### ###################################################################################

source('src/MULTIPLOT_ggplot.R')

#prepare input parameter list, based on which species are included in speciesList [defined in the site parameter file]
inputPars = prepSpecies(speciesList, years2CC, inundFun)

#Plots to examine species niche ranges, and maximum above & belowground standing biomass. Use to confirm parameters 
#
{
zstar=seq(-1,1.5,0.01)
plot(zstar, sapply(zstar, inputPars$nicheFUN[[1]]), type='l', xlab='Elevation (z*)', ylab='Niche Probability', main='Species niche')
cols = c('black','red','cyan','blue','orange','purple','darkgreen')
if(nSpecies>1){
  for(i in 2:length(speciesList)){
    lines(zstar, sapply(zstar, inputPars$nicheFUN[[i]]), col=cols[i])
  }
}

plot(zstar, sapply(zstar, inputPars$coverAgTransfer[[1]], cc=1)* sapply(zstar, inputPars$nicheFUN[[1]]), type='l', ylab='Biomass (g/cm2)', xlab='z*', main='Maximum Aboveground Biomass')
if(nSpecies>1){
  for(i in 2:length(speciesList)){
    lines(zstar, sapply(zstar, inputPars$coverAgTransfer[[i]], cc=1)*sapply(zstar,inputPars$nicheFUN[[i]] ), col=cols[i])
  }
}

plot(zstar,sapply(zstar, inputPars$coverAgTransfer[[1]], cc=1)* sapply(zstar, inputPars$nicheFUN[[1]])*sapply(inundFun(zstar*(TR/2)*100)*100,inputPars$rsFUN[[1]]), type='l', main='Maximum Belowground Biomass', ylab='Biomass (g/cm2)', xlab='z*' )
if(nSpecies>1){
  for(i in 2:length(speciesList)){
    lines(zstar, sapply(zstar, inputPars$coverAgTransfer[[i]], cc=1)*sapply(zstar,inputPars$nicheFUN[[i]] )*sapply(inundFun(zstar*(TR/2)*100)*100,inputPars$rsFUN[[1]]), col=cols[i])
  }
}
  
}



#initial elevation (cm, relative to MSL). Must be run with at least 2 initial elevations. Using only 1 will fail.
elev=z.init = seq(-TR*100*0.75, TR*100*0.75, 25)  #full range of tidal inundation
#elev=z.init = c(0.5,0.7, 0.8,0.9,1, 1.1, 1.5)*(TR/2)*100 #z from z*

#find initial elevation for a 50 yr calibration run using initial elevation derived from soil core accretion, historic SLR amount, and the NAVD88 elevation of the core 
elev=z.init= c(30, (minCal$elevMSL+MSL) - minCal$acc/10*50 - (MSL-historicSLR*50) )

seaLevel = function(t, SLR_rate)  t*SLR_rate 

 genSLR = function(yrs, SLR_rate){
   sea.Level = seaLevel(1:(yrs+1), SLR_rate)
   slrVar =  diff(sea.Level)
   slrVar
 }
 
 c.density=0.42

########################################
# RE-SOURCE SITE PARAMETERS, RUN MODEL #
########################################

#RUN CALIBRATION

#site parameters that the model is sensitive to. For exploration
# sitePars$kdec=0.0025  #daily decomposition decay rate 
# sitePars$porosityRate=0.03  #should be fast enough so that the last cohort has a porosity = minPorosity.

 
#RUN MODEL
{
  spinupYrs = 100 # Spinup allows soil core & species cover to establish equilibrium while holding elevation constant.  
 
   genSLRArray = function(z.init){
    array(sapply(1:length(z.init), function(x)  c(rep(historicSLR,spinupYrs+51))), dim=c(spinupYrs+51,length(z.init)))
  }
  
  SLR = genSLRArray(z.init)
  years = nrow(SLR)

  step=1
  mercc=NULL
  t.time = Sys.time()
  mercc = MainLoop_wSPINUP(elev, years, SLR, spinupYrs, inputPars,sitePars) 
  print(Sys.time()- t.time)
  
  startX=spinupYrs  #start yr for plotting
  endX= years   #end yr for plotting
  plotCal=T     #include additional plots for assessing calibration
  toFile=F      #output to file? 
  ncore=2       #which initial elevation to plot
  species=speciesList
  scenName = paste('calibration\\', site.name,"_Yrs",startX, "_",endX,"_SLRmm", round(mean(SLR),2)*10,"_Z", elev[ncore],
                   '_CAL_meanAcc_20240716.tif', sep='' )
 source('src/PLOT_THINGS_5sp.R')
}

 
#Plot calibration figures
#only plots the first 3 initial elevations

#cores=NULL
{
pyears= years-1
pOM = (mercc$cohorts.state$LOM+mercc$cohorts.state$ROM +apply(mercc$cohorts.state$FINEROOTS+mercc$cohorts.state$FRP, sum, MARGIN=1) )/(mercc$cohorts.state$LOM+mercc$cohorts.state$ROM+mercc$cohorts.state$MIN+apply(mercc$cohorts.state$FINEROOTS+mercc$cohorts.state$FRP, sum, MARGIN=1) )
depth=NULL
for(i in 1:nrow(mercc$cohorts.state$LOM)){
  depth =rbind(depth, cumsum(mercc$cohorts.state$VOL[i,]) )
}
# # plot(cores$Depth,cores$OM, ylim=c(0,100))
# # lines(pOM[1,]*100~ depth[1,], xlim=c(0,60), col='red')
# # lines(pOM[2,]*100~ depth[2,], xlim=c(0,60), col='green')
# # lines(pOM[3,]*100~ depth[3,],xlim=c(0,60), col='blue')
# # 
# # 
# # plot(cores$Depth, cores$BD, ylab='BD', xlab='Depth')
# # lines(mercc$cohorts.state$BD[1,]~depth[1,], col='red')
# # lines(mercc$cohorts.state$BD[2,]~depth[1,], col='green')
# # lines(mercc$cohorts.state$BD[3,]~depth[1,], col='blue')
 
 

#############
#OM and BD
###############
if(is.null((cores))){
plotOM = ggplot()  + geom_line(aes(x=depth[1,], y=pOM[1,]*100), col='red') +  geom_line(aes(x=depth[2,], y=pOM[2,]*100), col='darkgreen') +  geom_line(aes(x=depth[3,], y=pOM[3,]*100), col='blue') + xlim(c(0,60)) + ylim(c(0,100)) +
  theme_bw() + ylab('Organic Matter (%)') + xlab('Depth (cm)')

pBD= ggplot()  + geom_line(aes(x=depth[1,], y=mercc$cohorts.state$BD[1,]), col='red') +  geom_line(aes(x=depth[2,], y=mercc$cohorts.state$BD[2,]), col='darkgreen') +  geom_line(aes(x=depth[3,], y=mercc$cohorts.state$BD[3,]), col='blue') + xlim(c(0,60)) + #ylim(c(0,1)) +
  theme_bw() + ylab('Bulk Density (g/cm3)') + xlab('Depth (cm)')

} else {
  plotOM = ggplot() +geom_point(aes(x=cores$Depth, y=cores$OM)) + geom_line(aes(x=depth[1,], y=pOM[1,]*100), col='red') +  geom_line(aes(x=depth[2,], y=pOM[2,]*100), col='darkgreen') +  geom_line(aes(x=depth[3,], y=pOM[3,]*100), col='blue') +
    theme_bw() + ylab('Organic Matter (%)') + xlab('Depth (cm)')
  pBD= ggplot() +geom_point(aes(x=cores$Depth, y=cores$BD)) + geom_line(aes(x=depth[1,], y=mercc$cohorts.state$BD[1,]), col='red') +  geom_line(aes(x=depth[2,], y=mercc$cohorts.state$BD[2,]), col='darkgreen') +  geom_line(aes(x=depth[3,], y=mercc$cohorts.state$BD[3,]), col='blue') + xlim(c(0,60)) + #ylim(c(0,1)) +
    theme_bw() + ylab('Bulk Density (g/cm3)') + xlab('Depth (cm)')
  
  #Triangle
  # plotOM1 = ggplot() +geom_point(aes(x=cores$Depth, y=cores$OM)) + geom_line(aes(x=depth[1,], y=pOM[1,]*100), col='red') +   xlim(c(0,100)) + ylim(c(0,50)) +# geom_line(aes(x=depth[2,], y=pOM[2,]*100), col='darkgreen') +  geom_line(aes(x=depth[3,], y=pOM[3,]*100), col='blue')
  #   theme_bw() + ylab('Organic Matter (%)') + xlab('Depth (cm)')
  # pBD1= ggplot() +geom_point(aes(x=cores$Depth, y=cores$BD)) + geom_line(aes(x=depth[1,], y=mercc$cohorts.state$BD[1,]), col='red') +   xlim(c(0,100)) + #ylim(c(0,1)) +#geom_line(aes(x=depth[2,], y=mercc$cohorts.state$BD[2,]), col='darkgreen') +  geom_line(aes(x=depth[3,], y=mercc$cohorts.state$BD[3,]), col='blue') + xlim(c(0,60)) + #ylim(c(0,1)) +
  #   theme_bw() + ylab('Bulk Density (g/cm3)') + xlab('Depth (cm)')
  
  plotOM1 = ggplot() +geom_point(aes(x=cores$Depth, y=cores$OM)) + geom_line(aes(x=depth[1,], y=pOM[1,]*100), col='red') +   xlim(c(0,100)) + ylim(c(0,50)) +# geom_line(aes(x=depth[2,], y=pOM[2,]*100), col='darkgreen') +  geom_line(aes(x=depth[3,], y=pOM[3,]*100), col='blue')
    theme_bw() + ylab('Organic Matter (%)') + xlab('Depth (cm)')
  pBD1= ggplot() +geom_point(aes(x=cores$Depth, y=cores$BD)) + geom_line(aes(x=depth[1,], y=mercc$cohorts.state$BD[1,]), col='red') +   xlim(c(0,100)) + #ylim(c(0,1)) +#geom_line(aes(x=depth[2,], y=mercc$cohorts.state$BD[2,]), col='darkgreen') +  geom_line(aes(x=depth[3,], y=mercc$cohorts.state$BD[3,]), col='blue') + xlim(c(0,60)) + #ylim(c(0,1)) +
    theme_bw() + ylab('Bulk Density (g/cm3)') + xlab('Depth (cm)')
  
 
  }

#png('C:\\MERCC\\MERCC_SFBay\\output\\Triangle_BD-OMcal.png', res=800, units='mm', width=100, height=200, bg='transparent')
# png('C:\\MERCC\\MERCC_SFBay\\output\\Guad_BD-OMcal.png', res=800, units='mm', width=100, height=200, bg='transparent')
# png('C:\\MERCC\\MERCC_SFBay\\output\\Mowry_BD-OMcal.png', res=800, units='mm', width=100, height=200, bg='transparent')
# ggarrange(pBD1,plotOM1)
# dev.off()

#

# ncore=3
# 
# plot(mercc$vbigRoots[ncore,], col='red', type='l')
# lines(mercc$bigRoots[ncore,], col='green')
# lines(mercc$fineRoots[ncore,], col='blue')
# lines(mercc$deadRoots[ncore,], col='brown')

#(mercc$deadRoots[ncore,]/(mercc$bigRoots[ncore,]+mercc$fineRoots[ncore,]+mercc$deadRoots[ncore,]), type='l',
#    ylim=c(0,1), ylab='Dead:Live root biomass')   #0.86 - 0.93
# 
# 
# plot( (mercc$rom[1,]+mercc$lom[1,]+ mercc$deadRoots[1,]) /
#         (mercc$bigRoots[1,]+mercc$fineRoots[1,]+mercc$rom[1,]+mercc$lom[1,]+ mercc$deadRoots[1,]), type='l',
#      ylim=c(0.8,1), ylab='Dead:Live root biomass', col='brown')   #0.86 - 0.93
# lines( (mercc$rom[2,]+mercc$lom[2,]+ mercc$deadRoots[2,]) /
#          (mercc$bigRoots[2,]+mercc$fineRoots[2,]+mercc$rom[2,]+mercc$lom[2,]+ mercc$deadRoots[2,]), col='green')
# lines( (mercc$rom[3,]+mercc$lom[3,]+ mercc$deadRoots[3,]) /
#          (mercc$bigRoots[3,]+mercc$fineRoots[3,]+mercc$rom[3,]+mercc$lom[3,]+ mercc$deadRoots[3,]), col='blue')
# abline(h=0.86, col='red')
# abline(h=0.93, col='red')

##################
#Live:Dead roots 
##################

# core1 = (mercc$rom[1,]+mercc$lom[1,]+ mercc$deadRoots[1,]) /
#   (mercc$bigRoots[1,]+mercc$fineRoots[1,]+mercc$rom[1,]+mercc$lom[1,]+ mercc$deadRoots[1,])
# core2 = (mercc$rom[2,]+mercc$lom[2,]+ mercc$deadRoots[2,]) /
#   (mercc$bigRoots[2,]+mercc$fineRoots[2,]+mercc$rom[2,]+mercc$lom[2,]+ mercc$deadRoots[2,])
# core3 = (mercc$rom[3,]+mercc$lom[3,]+ mercc$deadRoots[3,]) /
#   (mercc$bigRoots[3,]+mercc$fineRoots[3,]+mercc$rom[3,]+mercc$lom[3,]+ mercc$deadRoots[3,])

core1 = ( mercc$deadRoots[1,]) /
  ( mercc$coarseRoots[1,]+mercc$fineRoots[1,]+mercc$smallRoots[1,]+ mercc$deadRoots[1,])
core2 = (mercc$deadRoots[2,]) /
  (mercc$coarseRoots[2,]+mercc$fineRoots[2,]+mercc$smallRoots[2,]+ mercc$deadRoots[2,])
core3 = ( mercc$deadRoots[3,]) /
  (mercc$coarseRoots[3,]+mercc$fineRoots[3,]+mercc$smallRoots[2,]+ mercc$deadRoots[3,])


#core1 = (mercc$deadRoots[1,]) /
#  (mercc$bigRoots[1,]+mercc$fineRoots[1,]+ mercc$deadRoots[1,])

pDeadLive = ggplot() + geom_line(aes(y=core1, x=1:years), col='red')+ geom_line(aes(y=core2, x=1:years), col='darkgreen')+ geom_line(aes(y=core3, x=1:years), col='blue') + theme_bw()+
  ylim(c(0, 1)) + ylab('Dead:Total Root Biomass') + xlab('Year')# + geom_hline(yintercept=calVals$deadTotalRootmin, linetype='dashed')+ geom_hline(yintercept=calVals$deadTotalRootmax, linetype='dashed')


#plot(mercc$roots[ncore,]+mercc$deadRoots[ncore,],type='l' ,ylab='Total Root Biomass g/cm2', main='ALL Roots (fine, big, vbig, dead)')  
# 
# plot(mercc$fineRoots[1,]+mercc$bigRoots[1,]+mercc$deadRoots[ncore,],col='brown', type='l' ,ylab='Total Root Biomass g/cm2', 
#      main='Cored Roots (fine, big, dead)', ylim=c(0,4))  #[1.3-3.7 g/cm2, total root biomass (<=20mm), Corimer et al 2015, Sapwalap]
# lines(mercc$fineRoots[2,]+mercc$bigRoots[2,]+mercc$deadRoots[2,], col='green')
# lines(mercc$fineRoots[3,]+mercc$bigRoots[3,]+mercc$deadRoots[3,], col='blue')
# abline(h=1.3, col='red')
# abline(h=3.7, col='red')

##################
#Total Roots
##################


core11 = mercc$fineRoots[1,]+mercc$smallRoots[1,]+ mercc$coarseRoots[1,] +mercc$deadRoots[1,] 
core22 = mercc$fineRoots[2,]+mercc$smallRoots[2,] + mercc$coarseRoots[2,] +mercc$deadRoots[2,]
core33 = mercc$fineRoots[3,]+mercc$smallRoots[3,]+ mercc$coarseRoots[3,] +mercc$deadRoots[3,] 

pTotRoots = ggplot() +  geom_line(aes(y=core11, x=1:years), col='red')+ geom_line(aes(y=core22, x=1:years), col='darkgreen')+ geom_line(aes(y=core33, x=1:years), col='blue') + theme_bw()+
  ylab('Total Root Biomass (g/cm2)') + xlab('Year') +xlim(c(0,pyears)) #+ geom_hline(yintercept=calVals$totRootMassmin, linetype='dashed')  + geom_hline(yintercept=calVals$totRootMassmax, linetype='dashed') 


# 
# plot(mercc$bigRoots[1,]+mercc$fineRoots[1,], col='brown', type='l', ylab='Living Root Biomass (g/cm2)', ylim=c(0,1))  # g/cm2  [0.04-.26 g/cm2, root (<=20mm) biomass, Corimer et al 2015, Sapwalap]
# lines(mercc$bigRoots[2,]+mercc$fineRoots[2,], col='green')
# lines(mercc$bigRoots[3,]+mercc$fineRoots[3,], col='blue')
# abline(h=0.04, col='red')
# abline(h=0.26, col='red')


##################
#Living Roots, fine and big only [not structual]
##################

core111=mercc$coarseRoots[1,]+mercc$fineRoots[1,]+mercc$smallRoots[1,]+mercc$strucRoots[1,]  #+ mercc$vbigRoots[1,]
core222 =mercc$coarseRoots[2,]+mercc$fineRoots[2,]+mercc$smallRoots[2,]+mercc$strucRoots[2,] #+ mercc$vbigRoots[2,]
core333= mercc$coarseRoots[3,]+mercc$fineRoots[3,]+mercc$smallRoots[3,]+mercc$strucRoots[3,] #+ mercc$vbigRoots[3,]

pLiveRoots = ggplot() +  geom_line(aes(y=core111, x=1:years), col='red')+ geom_line(aes(y=core222, x=1:years), col='darkgreen')+ geom_line(aes(y=core333, x=1:years), col='blue') + theme_bw()+
  ylim(c(0, 2))  +xlim(c(0,pyears))+ ylab('Living Root Biomass (g/cm2)') + xlab('Year')  #+geom_hline(yintercept=calVals$liveRootsmin, linetype='dashed')+ geom_hline(yintercept=calVals$liveRootsmax, linetype='dashed')+ geom_hline(yintercept=0.33, linetype='dotted') 


# 
# 
# plot( (mercc$downedWood_V[ncore,1:pyears]/mercc$woodDen[ncore,1:pyears] +
#          mercc$leaflitter[ncore,1:pyears]+mercc$woodlitter[ncore,1:pyears])*1e8/1000, type='l',
#       ylab='Litter kg ha yr-1')  #Kosrae-- 8800 kg/ha/yr; Gleason and Ewel 2002
# abline(h=8800, col='red')


##################
#Litter fall
##################


pyears = years-1
 #(mercc$downedWood_V[1,1:pyears]/mercc$woodDen[1,1:pyears] +
dcore1 =        (mercc$leaflitter[1,1:pyears]+mercc$woodlitter[1,1:pyears])*100 #+ mercc$downedWood_V[1,1:pyears]/mercc$woodDen[1,1:pyears]*100 # 1e8/1000
 #(mercc$downedWood_V[2,1:pyears]/mercc$woodDen[2,1:pyears] +
dcore2 =            (mercc$leaflitter[2,1:pyears]+mercc$woodlitter[2,1:pyears])*100# + (mercc$downedWood_V[2,1:pyears]/mercc$woodDen[2,1:pyears])*100 #1e8/1000
# (mercc$downedWood_V[3,1:pyears]/mercc$woodDen[3,1:pyears] +
dcore3 =           (mercc$leaflitter[3,1:pyears]+mercc$woodlitter[3,1:pyears])*100# +(mercc$downedWood_V[3,1:pyears]/mercc$woodDen[3,1:pyears])*100#1e8/1000

pLitter = ggplot() +  geom_line(aes(y=dcore1, x=1:pyears), col='red')+ geom_line(aes(y=dcore2, x=1:pyears), col='darkgreen')+ geom_line(aes(y=dcore3, x=1:pyears), col='blue') + theme_bw()+
  ylab('Litter (Mg/ha/yr)') + xlab('Year') +xlim(c(0,pyears)) #+ geom_hline(yintercept=calVals$littermin, linetype='dashed') + geom_hline(yintercept=calVals$littermax, linetype='dashed')    #range from Ding Darling -low:basin, high:fringe. Conrad 2022.

#+ geom_hline(yintercept=7.67, linetype='dashed') + geom_hline(yintercept=10.14, linetype='dashed') 
multiplot( pTotRoots, pLiveRoots,pDeadLive,plotOM, pBD, pLitter ,cols=2)  #plotOM, pBD, pLitter
}


#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################


#Projections

slr = read.csv("input\\SLR_projections_cm_Alameda.csv")

m1 = lm(Int~poly(Year,4), data=slr)
summary(m1)
slr.int=predict(m1, newdata=data.frame(Year=2020:2150))
m2 = lm(IntHigh~poly(Year,4), data=slr)
summary(m2)
slr.intHigh=predict(m2, newdata=data.frame(Year=2020:2150))
m3 = lm(High~poly(Year,4), data=slr)
summary(m3)
slr.High=predict(m3, newdata=data.frame(Year=2020:2150))

slr = data.frame(YEAR=2020:2150, Int=slr.int, IntHigh=slr.intHigh, High=slr.High)


inputPars = prepSpecies(speciesList, years2CC, inundFun)


z.init = seq(0, 200, 10)
#z.init = c(75, 90,100,110, 120, 130, 140,150,160 )

{
  #  source(paste('input/Sites/', site.name, '/', site.name, 'Parameters_bySpp.R', sep='') )
  spinupYrs = 300
  
  genSLRArray = function(z.init, slrProj){
    array(sapply(1:length(z.init), function(x)  c(rep(historicSLR,spinupYrs),slrProj)), dim=c(spinupYrs+131, length(z.init)))
  }
  
  # #fringe accretion rate= 0.11 cm/yr
  # genSLRArray = function(z.init, slr){
  #   array(sapply(1:length(z.init), function(x)  c( rep(0.11, spinupYrs), rep(0.11, 50), genSLR(200, 0.11) )), dim=c(301, length(z.init))) #,  genSLR(100, 0.25), genSLR(100, 0.35), (diff(slr)+0.074) ) ) , dim=c(980, length(z.init)) )
  #   #  c( rep(spinups[x]/1000, 1000), diff(slr)+0.074) ), dim=c(1080, length(z.init)) )
  # }
  
  SLR = cbind(genSLRArray(z.init, diff(slr$Int)), genSLRArray(z.init, diff(slr$IntHigh)), genSLRArray(z.init, diff(slr$High)) )
  elev = rep(z.init,3)
  #SLR = c(rep(.18, 2000) ) #, genSLR(50, 0.375), (diff(slr)+0.074))
  years = nrow(SLR)
  npixels = length(z.init)
  
  km=NULL
  mercc=NULL
  t.time = Sys.time()
  mercc = MainLoop_wSPINUP(elev, years, SLR, spinupYrs, inputPars) 
  print(Sys.time()- t.time)
  
  startX=1 #spinupYrs+1
  endX= years
  plotCal=T  #include additional plots for assessing calibration
  toFile=F
  ncore=3 #which initial elevation to plot
  species=speciesList
#  source('src/PLOT_THINGS_5sp.R')
}


plot(x=2020:2150, y=mercc$elev_msl[4,(spinupYrs+1):(spinupYrs+131)], ylim=c(-200,160),type='l', col='darkblue', ylab='Elevation (cm, MSL)', xlab='Year')
abline(h=0, lty='dotted')
abline(h=70, lty='dashed')
lines(x=2020:2150,mercc$elev_msl[4+9,(spinupYrs+1):(spinupYrs+131)], ylim=c(-200,160), col='blue')
lines(x=2020:2150,mercc$elev_msl[4+18,(spinupYrs+1):(spinupYrs+131)], ylim=c(-200,160), col='lightblue')

z.targ=5

plot(mercc$perCover[z.targ,1,(spinupYrs+1):(spinupYrs+131)], type='l', ylim=c(0,1), lwd=2, xlab='Years', ylab='Cover')
lines(mercc$perCover[z.targ,2,(spinupYrs+1):(spinupYrs+131)], lwd=2, col='red')
lines(mercc$perCover[z.targ,1,(spinupYrs+1):(spinupYrs+131)]+mercc$perCover[z.targ,2,(spinupYrs+1):(spinupYrs+131)]+mercc$perCover[z.targ,3,(spinupYrs+1):(spinupYrs+131)]+mercc$perCover[z.targ,4,(spinupYrs+1):(spinupYrs+131)], lwd=1.5, lty='dotted', col='darkgreen')


saveRDS(mercc, 'output\\Restoration_SLR_Projections.rds')


