#!/usr/bin/env Rscript

library(terra)
library(plyr)
library(data.table)


tallgrass=F


site=site.name='yourSiteName'

setwd('pathToWorkingDir')

initialDEM = rast("pathToDEM.tiff")
vals = values(initialDEM)*100
saveRDS(vals, paste0('input/',site.name,'/',site.name,'DEM_values.rds'))

source(paste('input/', site.name, '/', site.name, 'Parameters.R', sep='') )

initial = readRDS( paste0('input/',site.name, '/', site.name, 'DEM_values.rds' ))-MSL


files = list.files(paste0('output/', site.name), full.names=T)


#find number of SLR scenarios, initial elevations, number of species from model output
data = readRDS(files[1])[[1]]
spinupYrs =max(table(data$SLR[,1]))-1
nSLR = length(unique(data$SLR[spinupYrs+2,]))
nZ=length(data$elev_msl[,1])/nSLR
z.init = data$elev_msl[1:nZ,1]
nSpecies =dim(data$perCover)[2]/2
projectionI = spinupYrs:ncol(data$elev_msl)

elevP=carbP=perCovP=vector('list', nSLR)
for (s in 1:nSLR) {
  elevP[[s]] <- list()
  carbP[[s]] <- list()
  perCovP[[s]] <- list()
}
 
k=1
for(i in 1:length(files)){
  print(i)
  data=readRDS(files[i])[[1]]
  if( max(abs(range(data$elev_msl)))<1000 ){
    for(s in 1:nSLR){
      elevP[[s]][[k]]= data$elev_msl[(1:length(z.init))+length(z.init)*(s-1),projectionI]
      carbP[[s]][[k]]= t(apply(data$carbon[(1:length(z.init))+length(z.init)*(s-1),projectionI], MARGIN=1, FUN=diff))
      perCovP[[s]][[k]] = data$perCover[(1:length(z.init))+length(z.init)*(s-1),,projectionI]
    }
    k=k+1
  }
}

calcMean = function(matrix_list, n){
  meanMat = Reduce('+',matrix_list)/n
  sum_sq_diff <- Reduce(`+`, lapply(matrix_list, function(mat) (mat - meanMat)^2))
  sd_matrix <- sqrt(sum_sq_diff / (length(matrix_list)))
  
  return(list(mean=meanMat, sd=sd_matrix))
}

  #Function which takes in an initial elevation, the low and high elevations at which the model was run and the corresponding results for those elevation and 
  #uses a linear interpolation function to determine the new elevation from model results
  
interp.Raster = function(z.init, #initial elevations that the model was run at. cm
                         years,  #vector of years the simulation represents eg, 2020, 2021...
                         speciesList,  #vector of species names
                         demMSL,  # rast file of the DEM, relative to MSL in cm
                         theZData,  #elevation projections, list with mean and SD
                         theCData,  #carbon accumulation projections, list with mean and SD
                         thePCData, #species cover projections, list with mean and SD
                         outputRasterFreq,  #frequnecy of raster output. eg, every 5, 10 20 yrs. if 0, no rasters output
                         outFolderRaster,  #folder path for raster output. NA if not running wanting raster output
                         presenceThres      # proportion cover threshold for species. eg, assume that > 0.15 species are present. Used in aggregated cover results across the DEM.
                         ){  
  
  L.results_Z = as.data.frame(t(theZData$mean))
  L.results_C = as.data.frame(t(theCData$mean))
  ZSD = as.data.frame(t(theZData$sd))
  CSD = as.data.frame(t(theCData$sd))
  
  pcSD = thePCData[[2]]
  thePCData=thePCData[[1]]
  
  nSpecies=length(speciesList)
  
  if(nSpecies>1){
    living = apply(thePCData[,1:nSpecies,]/sitePars$CC, MARGIN=c(1,3), FUN=sum)
    livingSD =  apply(pcSD[,1:nSpecies,]/sitePars$CC, MARGIN=c(1,3), FUN=function(x) sqrt(sum(x^2)))
  
    dead = apply(thePCData[,(nSpecies+1):(nSpecies*2),]/sitePars$CC, MARGIN=c(1,3), FUN=sum)
    deadSD =  apply(pcSD[,(nSpecies+1):(nSpecies*2),]/sitePars$CC, MARGIN=c(1,3), FUN=function(x) sqrt(sum(x^2)))
  } else {
    living = thePCData[,1,]/sitePars$CC
    livingSD =pcSD[,1,]/sitePars$CC
    
    dead = thePCData[,2,]/sitePars$CC
    deadSD =pcSD[,2,]/sitePars$CC
  }
  #######################
  #Code to extrapolate model results to DEM, loops through each SSC and SLR combination
  
  output.table = NULL# matrix(ncol=12+length(speciesList)*2, nrow=length(years))
  
  ################
  outyrs = seq(1, length(years), 1) #201:301   #output only every 10 years
  outyrsRaster=1
  if(outputRasterFreq>0){
    outyrsRaster = seq(1, length(years), outputRasterFreq)
  }
  year = seq(2020, 2019+length(years), 1)

  int_elev =  values(demMSL) #extract values from DEM
  resolution = terra::res(demMSL)[1] #resolution in meters
  totArea_ha=length(na.omit(int_elev))*resolution*resolution/10000  #area in hectares of total domain [not including NAs]
  
  for(k in 1:length(outyrs)){
    print(k)
    
    #relate initial elevation and model projections
    results = approxfun(z.init, L.results_Z[k,], rule=2 )(int_elev)
    resultsSD = approxfun(z.init,ZSD[k,], rule=2 )(int_elev)

    projDead = approxfun(z.init, dead[,k], rule=2)(int_elev)
    projDeadSD = approxfun(z.init, deadSD[,k], rule=2)(int_elev)
    
    projLiving = approxfun(z.init, living[,k], rule=2)(int_elev)
    projLivingSD = approxfun(z.init, livingSD[,k], rule=2)(int_elev)

    spProj=spSDFUN = list()
    for(sp in 1:nSpecies){
      spProj[[sp]] = approxfun(z.init, thePCData[,sp,k]/sitePars$CC, rule=2)(int_elev)
      spSDFUN[[sp]] = approxfun(z.init, thePCData[,sp+nSpecies,k]/sitePars$CC, rule=2)(int_elev)
    }

    #Carbon projections
    if(k<length(years)){  
      resultsC = approxfun(z.init,L.results_C[k,], rule=2)(int_elev)
      resultsCSD = approxfun(z.init, CSD[k,], rule=2)(int_elev)
    } else {
      resultsC=resultsCSD=NA
    }
    
    spCov= spCovSD=NULL
    for(sp in 1:nSpecies){
      spCov= c(spCov,length(which(na.omit(spProj[[sp]])>presenceThres))/length(na.omit(spProj[[sp]])) )
      spCovSD = c(spCovSD,  abs(spCov[[sp]]-length(which(na.omit(spProj[[sp]]-spProjSD[[sp]])>presenceThres ))/length(na.omit(spProj[[sp]]))))
    }
    
    livingPer = sum(ifelse(living>0.05,1,0))/length(living) 
    deadPer = 1-livingPer #sum(ifelse(dead>0.05,1,0))/length(dead) 
    
   # deadPerSDmin = length(which(na.omit(projDead-projDeadSD)>presenceThres ))/length(na.omit(projDead))
    livingPerSDmin = length(which(na.omit(projLiving-projLivingSD)>presenceThres ))/length(na.omit(projLiving))
    
    wetland_ha = sum(ifelse(living>0.05,1,0))/length(living) *totArea_ha
    wetland_ha_sd = ifelse(wetland_ha/totArea_ha-livingPerSDmin>0.05, 1, 0) *totArea_ha
    
    #create and output rasters
    if(outputRasterFreq>0 & k %in% outyrsRaster){
      out_r = setValues(initial, as.vector(unlist(results)))
      out_c = setValues(initial, as.vector(unlist(resultsC))*10000)  #C accumulation in g/m2/yr
      
      spBrick = terra::setValues(x=initial, value=spProj[[1]])
      spSDBrick = setValues(x=initial, value=spProjSD[[1]])
      for(sp in 2:nSpecies){
        spBrick = c(spBrick, terra::setValues(x=initial, value=spProj[[sp]] ))
        spSDBrick = c(spSDBrick, terra::setValues(x=initial, value=spProjSD[[sp]] ))
      }
      
      writeRaster(out_r, paste(outFolderRaster, "ElevationMSL_", year[k], "_SLR_", slr,".tif", sep=""), overwrite=TRUE, format='GTiff')
      writeRaster(out_c, paste(outFolderRaster, "CarbonAccum_",year[k], "_SLR_", slr,".tif", sep=""), overwrite=TRUE, format='GTiff')#, datatype='INT4S')
      writeRaster(spBrick, paste(outFolderRaster, "SpeciesCover_",year[k], "_SLR_", slr,".tif", sep=""), overwrite=TRUE, format='GTiff')
    }
    
     
    output.table = rbind(output.table, c(year[k], mean(unlist(results), na.rm=T),#mean((unlist(results)-(diurnal.TR)), na.rm=T), 
                         sd(unlist(results),na.rm=T), mean(unlist(resultsC), na.rm=T)*100,  # g/cm3 -> g/ha / g->Mg  1e8/1e6,
                         sd(unlist(resultsC)*100, na.rm=T), sum(unlist(resultsC*10000)*res1*res1, na.rm=T)/1e6,
                         spCov, spCovSD, livingPer, deadPer, wetland_ha, wetland_ha_sd))
                           
   

    
  }
    colnames(output.table) = c("Year","MeanElev_MSL","SDElev", 'Carbon_Accum_Mg_ha_yr','Carbon_Accum_Mg_ha_yr_SD','TotalCarbon_Mg_yr',
                              speciesList, paste0(speciesList, '_sd'), 'Dead', 'Living', 'Wetland_ha','Wetland_ha_sd' )#, 'Carbon_Accum_g_m2_yr_CI')
    
   return(output.table)
}

 
  #calculate mean & SD from model results. Also works if no Monte Carlo was done
  elevMean=carbMean=perCovMean=list()
  for(s in 1:nSLR){
    elevMean[[s]]=  calcMean(elevP[[s]],length(elevP[[s]]))
    carbMean[[s]]=  calcMean(carbP[[s]],length(carbP[[s]]))
    perCovMean[[s]]=  calcMean(perCovP[[s]],length(perCovP[[s]]))
  }
  
  #initial_org is the DEM in cm, relative to NAVD88
  
  vlowInterp = interp.Raster(z.init, 2020:2150, speciesList,initial_org-MSL,  elevMean[[1]], carbMean[[1]], perCovMean[[1]], 0,NA, 0.15)
  lowInterp = interp.Raster(z.init, 2020:2150, speciesList,initial_org-MSL,  elevMean[[2]], carbMean[[2]], perCovMean[[2]], 0,NA, 0.15)
  midInterp = interp.Raster(z.init, 2020:2150, speciesList,initial_org-MSL,  elevMean[[3]], carbMean[[3]], perCovMean[[3]], 0,NA, 0.15)
  midHighInterp = interp.Raster(z.init, 2020:2150, speciesList,initial_org-MSL,  elevMean[[4]], carbMean[[4]], perCovMean[[4]],0,NA, 0.15)
  highInterp = interp.Raster(z.init, 2020:2150, speciesList,initial_org-MSL,  elevMean[[5]], carbMean[[5]], perCovMean[[5]], 0,NA, 0.15)
  exInterp = interp.Raster(z.init, 2020:2150, speciesList,initial_org-MSL,  elevMean[[6]], carbMean[[6]], perCovMean[[6]], 0,NA, 0.15)
  
  vlowInterpdf = as.data.frame(vlowInterp)
  lowInterpdf = as.data.frame(lowInterp)
  midInterpdf = as.data.frame(midInterp)
  midHighInterpdf = as.data.frame(midHighInterp)
  highInterpdf = as.data.frame(highInterp)
  exInterpdf = as.data.frame(exInterp)
  
  vlowInterpdf$SLR='SSP 1-1.9'
  lowInterpdf$SLR = 'SSP 1-2.6'
  midInterpdf$SLR = 'SSP 2-4.5'
  midHighInterpdf$SLR  = 'SSP 3-7.0'
  highInterpdf$SLR = 'SSP 5-8.5'
  exInterpdf$SLR = 'SSP 5-8.5 lc'
  
  allResults = rbind(vlowInterpdf, lowInterpdf, midInterpdf, midHighInterpdf, highInterpdf, exInterpdf)
  
  saveRDS(allResults, paste0(site.name, '_interpolatedResults_MEAN_20240814.rds'))
  