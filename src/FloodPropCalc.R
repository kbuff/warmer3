## calc inundation characteristics on Pohnpei
setwd('C:\\MERCC\\MERCC_FL')

harmonic.constituents = read.csv("input\\Naples_TidalHarmonics_m_MSL.csv" ) #msarc.input$harmonic.constituents

tide <- Vectorize(function(t) with(harmonic.constituents, sum(Amplitude * cos((Speed * pi / 180) * t + Phase))*100 ))     #without MSL -- HC relative to MSL, in cm

getLevels = function(tide){
  k=1
  day.int = 365*24 #round(365/timestep) #365 #interval
  # start = 0  #173 best for week, 50.125 best for day
  end = 1+day.int
  int =  15/(60) ##12 secs  1 minute #0.25 #15 minutes
  interval = seq(1,end,int)  #closest daily interval to full year- 128-128.5, closest monthly interval-- 64-76---
  interval.sec = ((end-1)*60*60)/(length(interval)-1) ## .5 (1 day), 43200 s,
 # sfInit(parallel=T, cpus=7)
#  sfExport('harmonic.constituents')
  level = sapply(interval, tide)
 # sfStop()
  levels=list()
  levels$water=level
  levels$sec=interval.sec
  levels$int = int
  levels
}

theLevels = getLevels(tide)

plot(theLevels$water[1:(4*24*7)], type='l')

zRange <- seq(from=-100, to=200, by=0.05)
inundation.duration<- sapply(zRange, function(z) sum(theLevels$water > z, na.rm=T))
inundDur = inundation.duration/max(inundation.duration)
plot(zRange,inundDur, type='l')
floodingFunDur = approxfun(inundDur~zRange, rule=2)  #flooding duration

saveRDS(floodingFunDur, 'C:\\MERCC\\MERCC_FL\\input\\Naples_FloodingDuration.rds')

# Finds peaks from water level column by finding the max value around a moving window of 90 minutes
peaks = NA #matrix(NA, nrow=nrow(data), ncol=ncol(data))
count = 1
wl_level=2

for(i in 30:(length(theLevels$water)-30)){
  now = as.numeric(theLevels$water[i])
  before = max(as.numeric(theLevels$water[(i-1):(i-30)]))	
  after = max(as.numeric(theLevels$water[(i+1):(i+30)]))
  if(now > before & now > after){
    peaks = rbind(peaks, theLevels$water[i])
    peaks[count,] = theLevels$water[i]
  }
}

zRange <- seq(from=-100, to=200, by=0.05)
inundation.frequencies <- sapply(zRange, function(z) sum(peaks> z, na.rm=T))
inundProp = inundation.frequencies/max(inundation.frequencies)
plot(zRange,inundProp, type='l')
floodingFun = approxfun(inundProp~zRange, rule=2)  #flooding freq

#plot(zRange, I(inundDur/inundProp), type='l')

#max(theLevels$water)

