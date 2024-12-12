
working.dir = 'C:\\MERCC\\WARMER_v3\\warmer_v3_github_20241109'# 'pathTo_WARMER_v3_folder'
setwd(working.dir)

library(ggplot2)
library(reshape2)
library(foreach)
library(doParallel)
library(snowfall)
library(doSNOW)
library(plyr)
library(tictoc)
require(pracma)
library(abind)
library(egg)
library(grid)
library(openxlsx)
library(gam)
library(abind)
library(sf)
library(testthat)

source('src/GeneralFunctions.R')
source('src/MainLoop.R')
source('src/Belowground.R')
source('src/updateDepthPOR.R')
source('src/spCover.R')
source('src/SSC_Deposition_Marani2010.R')

out.folder = '\\output'
site.name = 'ChinaCamp'# 'PIE'  #'Marsh'

print(site.name)
source(paste('input/',site.name,'/', site.name, '_Parameters.R', sep='') )

#if 'user defined, the harmonic constituents must be loaded [HC] in the site parameter file.
#Otherwise, 1 yr of tidal inundation is simulated based on the provided tidal range and the tide type
tide_type = 'user-defined'# 'mixed semi-diurnal'#  'semi-diurnal' #'diurnal'

source('src/GenTidalSignal.R')

#source after inundFun and TR are defined
source('src/calibrate_sediment.R')
source('src/calibrate_porosity.R')

#update input files with calibrated values
ssc
ssc_sd
minPorosity
maxPorosity
porosityRate
OMden

##############################################################################################################
#RUN CALIBRATION MODELS WITH SPINUP
########################### ###################################################################################

source('src/MULTIPLOT_ggplot.R')

#prepare input parameter list, based on which species are included in speciesList [defined in the site parameter file]
#inputPars = prepSpecies(speciesList, years2CC, inundFun)

#Plots to examine species niche ranges, and maximum above & belowground standing biomass. Use to confirm parameters
#
{
  par(mfrow=c(1,3))
  zstar=seq(-1,1.5,0.01)
  plot(zstar, sapply(zstar, inputPars$nicheFUN[[1]]), type='l', xlab='Elevation (z*)', ylab='Niche Probability', main='Species niche')
  cols = c('black','red','cyan','blue','orange','purple','darkgreen')
  if(nSpecies>1){
    for(i in 2:length(speciesList)){
      lines(zstar, sapply(zstar, inputPars$nicheFUN[[i]]), col=cols[i])
    }
  }
  
  plot(zstar, sapply(sapply(zstar*(TR/2)*100, inundFun), maxB=inputPars$spMaxBiomass[1], FUN= inputPars$biomassFUN[[1]], cc=1) * sapply(zstar, inputPars$nicheFUN[[1]]), type='l', ylab='Biomass (g/cm2)', xlab='z*', main='Maximum Aboveground Biomass')
  if(nSpecies>1){
    for(i in 2:length(speciesList)){
      lines(zstar, sapply(sapply(zstar*(TR/2)*100, inundFun), maxB=inputPars$spMaxBiomass[i], inputPars$biomassFUN[[i]], cc=1)*sapply(zstar,inputPars$nicheFUN[[i]] ), col=cols[i])
    }
  }
  
  
  plot(zstar,sapply(sapply(zstar*(TR/2)*100, inundFun),maxB=inputPars$spMaxBiomass[1], inputPars$biomassFUN[[1]], cc=1)* sapply(zstar, inputPars$nicheFUN[[1]])*sapply(inundFun(zstar*(TR/2)*100)*100, rs=inputPars$root_shoot[1], inputPars$rsFUN[[1]]), type='l', main='Maximum Belowground Biomass', ylab='Biomass (g/cm2)', xlab='z*' )
  if(nSpecies>1){
    for(i in 2:length(speciesList)){
      lines(zstar, sapply(sapply(zstar*(TR/2)*100, inundFun), maxB=inputPars$spMaxBiomass[i], inputPars$biomassFUN[[i]], cc=1)*sapply(zstar,inputPars$nicheFUN[[i]] )*sapply(inundFun(zstar*(TR/2)*100)*100,rs=inputPars$root_shoot[i], inputPars$rsFUN[[i]]), col=cols[i])
    }
  }
}

#initial elevation (cm, relative to MSL). Must be run with at least 2 initial elevations. Using only 1 will fail.
elev=z.init = seq(0, TR*100*0.75, 25)  #full range of tidal inundation
#elev=z.init = c(0.5,0.7, 0.8,0.9,1, 1.1, 1.5)*(TR/2)*100 #z from z*

#find initial elevation for a 50 yr calibration run using initial elevation derived from soil core accretion,
# historic SLR amount, and the NAVD88 elevation of the core
elev=z.init= c((minCal$elevMSL+MSL) - minCal$acc/10*50 - (MSL-historicSLR*50) )

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


#site parameters that the model is sensitive to. For exploration
# sitePars$kdec=0.0025  #daily decomposition decay rate

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
  
  startX=1#spinupYrs  #start yr for plotting
  endX= years   #end yr for plotting
  plotCal=T     #include additional plots for assessing calibration
  toFile=F      #output to file?
  ncore=1       #which initial elevation to plot
  species=speciesList
  scenName = paste('calibration\\', site.name,"_Yrs",startX, "_",endX,"_SLRmm", round(mean(SLR),2)*10,"_Z", elev[ncore],
                   '_CAL_meanAcc_20240716.tif', sep='' )
  source('src/PLOT_THINGS_5sp.R')
}


#Compare against soil core accumulation & accretion rates
#RMSE Accretion rate
accRate=apply(mercc$acc[,(spinupYrs+1):(spinupYrs+50)], FUN=mean, MARGIN=1)
accRateDif=mean(accRate - minCal$acc)

#RMSE Organic Accumulation
omAccum=colMeans(apply(mercc$allOM[,(spinupYrs+1):(spinupYrs+50)], FUN=diff, MARGIN=1))
minAccum=colMeans(apply(mercc$allMIN[,(spinupYrs+1):(spinupYrs+50)], FUN=diff, MARGIN=1))
MinAccDiff=mean(minAccum-minCal$minAcc)*10000
OMAccDiff=mean((omAccum-minCal$omAcc)*10000)

#png(paste0('calibration\\',site.name, '_hindcastskill.png'), res=300, units='mm', height=200, width=150)
par(mfrow=c(3,1))
plot(minCal$minAcc*10000,minAccum*10000, ylim=c(0,2000), xlim=c(0,2000), pch=16 , main=paste0('Min accumulation g/m2/yr ', round(MinAccDiff,2)), ylab='Modeled',xlab='Observed')
abline(a=0,b=1)
plot(minCal$omAcc*10000,omAccum*10000, ylim=c(0,1000), xlim=c(0,1000), pch=16 , main=paste0('OM accumulation g/m2/yr ', round(OMAccDiff,2)), ylab='Modeled',xlab='Observed')
abline(a=0,b=1)
plot(minCal$acc,accRate, ylim=c(0,0.7), xlim=c(0,.7), pch=16 , main=paste0('Accretion cm/yr ', round(accRateDif,3)), ylab='Modeled',xlab='Observed')
abline(a=0,b=1)
# dev.off()

plot(mercc$elev_msl[,spinupYrs], rowSums(mercc$totABG[,,spinupYrs]*10000))
plot(mercc$elev_msl[,spinupYrs+50], rowSums(mercc$totABG[,,spinupYrs+50]*10000))
plot(mercc$elev_msl[,spinupYrs+50], (mercc$totBGB[,spinupYrs+50])*10000, ylim=c(0,4000))


#Plot calibration figures
source('src/PLOT_calibration.R')


#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################


#Generate Projections under SLR. 
source('src/getSLRScenario.R') #IPCC AR6 SLR projections. Uses lat,long of the site to select nearest SLR projection, relative to 2020. 
#Projections available at long-term gauges & interpolated at 1* resolution  

SLR_scenario = c('Int','IntHigh','High')  #vector of the SLR scenarios you want to run. Valid names:  'Low' 'IntLow' 'Int' 'IntHigh' 'High'
SLRYears=2020:2150  #vector of years for SLR projection. Can go up to 2150

slrDF=slrDF[slrDF$Year==SLRYears,]

z.init = seq(-150, 200, 5) # Elevation relative to MSL in cm. For generating projections for interpolating across DEM.
z.init =c(0,0.25, 0.5, 0.75, 1, 1.25)*(TR/2)*100#  Elevation relative to MSL in cm. For testing marsh response to SLR across a range of elevations

#runs all initial elevations and SLR scenarios at the same time.
{
  spinupYrs =100
  slVar=0#(rnorm(nrow(slr), mean=0, sd=0.0297))
  
  genSLRArray = function(z.init, slrProj){
    array(sapply(1:length(z.init), function(x)  c(rep(historicSLR,spinupYrs),slrProj[1],slrProj)), dim=c(spinupYrs+(nrow(slrDF)), length(z.init)))
  }
  
  SLR=NULL
  for(i in 1:length(SLR_scenario)){
    SLR=cbind(SLR, genSLRArray(z.init, diff(slrDF[[SLR_scenario[i]]]+slVar)*100))
  }
  elev = rep(z.init,3)
  years = nrow(SLR)
  mercc=NULL
  t.time = Sys.time()
  mercc = MainLoop_wSPINUP(elev, years, SLR, spinupYrs, inputPars, sitePars)
  print(Sys.time()- t.time)
}

#Collate results for plotting
elevProj=list() 
int=1
for(i in 1:length(SLR_scenario)){
  elevProj[[SLR_scenario[i]]]= cbind(SLRYears, t(mercc$elev_msl[int:(length(z.init)*i),(spinupYrs+1):(spinupYrs+nrow(slrDF))] ))
  colnames( elevProj[[SLR_scenario[i]]])=c('Year',paste0('z',z.init))
  int=length(z.init)*i+1
}

zPlots=list()
for(i in 1:length(SLR_scenario)){
  dat=elevProj[[SLR_scenario[i]]]
  dat2=melt(as.data.frame(dat), id.vars='Year', measured.vars=colnames(dat)[-1],  value.name='z', variable.name = 'Initialz')
  zPlots[[i]]=  ggplot(dat2) + geom_line(aes(x=Year,y=z, col=Initialz)) + theme_classic() + ggtitle(SLR_scenario[i])+ ylab('Elevation cm MSL')
}

#plot elevation results by SLR scenario
ggarrange(plots=zPlots, widths=c(1,1))


z.targ=3+length(z.init)*2

plot(mercc$perCover[z.targ,1,(spinupYrs+1):(spinupYrs+nrow(slrDF))], type='l', ylim=c(0,1), lwd=2, xlab='Years', ylab='Cover')
lines(mercc$perCover[z.targ,2,(spinupYrs+1):(spinupYrs+nrow(slrDF))], lwd=2, col='red')
lines(mercc$perCover[z.targ,3,(spinupYrs+1):(spinupYrs+nrow(slrDF))], lwd=2, col='cyan')
lines(mercc$perCover[z.targ,1,(spinupYrs+1):(spinupYrs+nrow(slrDF))]+mercc$perCover[z.targ,2,(spinupYrs+1):(spinupYrs+nrow(slrDF))]+mercc$perCover[z.targ,3,(spinupYrs+1):(spinupYrs+nrow(slrDF))], lwd=1.5, lty='dotted', col='darkgreen')


theDate=Sys.Date()
saveRDS(mercc, paste0('output/',site.name,'_SLR_Projections_',theDate,'.rds') )


