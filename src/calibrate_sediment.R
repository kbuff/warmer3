#calibrate sediment


## calculate sediment deposition based on tide signal with SSC=1 mg/L -- takes a bit, uses parallel processing

if(is.na(sitePars$ssc) | !exists("sedinundFun")){
  tic()
  sedinundFun= genMineralFunction(theLevels)  #deposition fun; input= elevation cm rMSL, output: mineral deposition, assuming SSC=1.
  toc()


  #input is elevation (cm, MSL) and standing aboveground biomass. Multiply by the SSC [mg/L]
  ztemp = seq(min(theLevels$water)*100, max(theLevels$water)*100 , 1)
  #plot(ztemp, sapply(ztemp, sedinundFun, 1000)*30, type='l', xlab='Elevation (cm, MSL)', ylab='Sediment Deposition (g/cm2/yr)', main='SSC = 30 mg/L')
  #plot(ztemp, sapply(ztemp, inundFun)*100, type='l', xlab='Elevation (cm, MSL)', ylab='Flooding Duration (% time)')

  ##########################################################################################################
  #Calibrate SSC against soil core mineral deposition. If no soil core data, set SSC based on best available data

  histElev=c((minCal$elevMSL+MSL) - minCal$acc/10*50 - (MSL-historicSLR*50) )

 # ssc=mean(minCal$minAcc/sapply(minCal$elevMSL, sedinundFun, B1=1000))
 # ssc=mean(minCal$minAcc/sapply( (histElev + minCal$elevMSL)/2, sedinundFun, B1=1000))
  ssc=mean(minCal$minAcc/sapply( (histElev), sedinundFun, B1=1000))
  plot(minCal$elevMSL, minCal$minAcc, xlab='Elevation (cm, MSL)', ylab='Sediment deposition g/cm2/yr', main=paste0('SSC= ', round(ssc,1)), xlim=c(min(ztemp), max(ztemp)), ylim=c(0,1))
  lines(ztemp, sapply(ztemp, sedinundFun, 1000)*(ssc))

  #Calculate the standard deviation of SSC using rates from multiple soil cores.
  if(length(minCal$minAcc)>1){
    ssc_temp=NULL
    for(i in 1:length(minCal$minAcc)){
      ssc_temp=c(ssc_temp, minCal$minAcc[i]/sapply(histElev[i], sedinundFun, B1=1000))
    }
    ssc_sd=sd(ssc_temp)
  } else {
    ssc_sd=5
  }

}
