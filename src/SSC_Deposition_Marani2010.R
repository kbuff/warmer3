
require(foreach)
require(doParallel)

#calculate sediment deposition.  
#deposition sensitive to settling velocity




# Sedimentation function with erosion

genMineralFunction = function(theLevels){
  
  mineralBiomass = function(theLevels){
    timestep = 1 #1 yr
   
    #UNITS: z1= m, SSC.i= g/ml, B=g/m2  theLevels$water=m
    SSC.out <- function(z1, SSC.i, B, theLevels){
      
     SSC.i=SSC.i/1000  #convert mg/l to kg/m3
      ws= 0.2/1000  #m/s    0.2mm/s  
      U = 0.02  #horizontal flow across the marsh, m/s  -- not used
     
      #veg interception-- for marshes. NEED TO ADJUST FOR MANGROVES
      d50=50/1000000 #50 um in m
      beta = 0.382
      alpha = ((U^1.7)*d50^2*(1.02*10^6))
       
      dep.tot = 0
      SSC.inst = SSC.i * .75  #for when the tide function begins with a negative slope
      
      interval.sec=theLevels$sec
      levels_1 = theLevels$water
      SSC.i.temp = SSC.i #  abs(rnorm(length(closed1), mean=SSC.i, sd=0))# 4.8/1000)) # rep(SSC.i, length=length(closed1) )
      
      for(i in 1:(length(levels_1)-1)){  ##.5 = 1 day, 12= 1 month,    0.01 = 180 seconds, 3 min intervals
        #For stations with tidal constituents in cm above MSL
        h=levels_1[i]  #   function is tidal input, water height above MSL, in cm (from tidesandcurrents.noaa.gov)
        h2 = levels_1[i+1]
        depth=h-z1  
        depth1=h2-z1
        depth.avg=(depth1+depth)/2
        
        #underwater with positive slope [flood tide], constant SSC
        if(depth.avg>0.03 & h2>=h ){ 
          SSC.inst = SSC.i.temp                        #remove deposited sediment from water column
          dep= ((-ws*SSC.inst - (alpha*(B^beta)*SSC.inst) )/interval.sec)/depth.avg
          dep.tot = dep.tot+dep*-1*interval.sec  #add deposited sediment to total deposited sediment pool  # 2.65 g/m3 sediment density, convert to g/cm2 with /10000 
          
          
          #negative slope [ebb tide] and a positive inst SSC from the last calculation, SSC decreases  
        }else if(depth.avg>0.03 & h2<h & SSC.inst>0 ){
          dep = ((-ws*SSC.inst - (alpha*(B^beta)*SSC.inst) )/interval.sec)/depth.avg
          ssc.calc= ((-ws*SSC.inst - (alpha*(B^beta)*SSC.inst)   + SSC.inst*(h2-h))/interval.sec)/depth.avg   #removed depth.avg
          SSC.inst = SSC.inst + ssc.calc*interval.sec
          dep.tot = dep.tot+(dep*-1)*interval.sec
        }
      }
      return(dep.tot*0.1) #kg/m2 to  g/cm2
    }
    
    out.Fun = function(inputs){
      SSC.out(inputs$z, inputs$SSC.i,inputs$B, theLevels)#*12 # inputs$B,
    }
    return(out.Fun)
  }  #, z, SSC.i, B
  
 
  #run across a range of elevations and aboveground biomasses
  z = seq(min(theLevels$water)*2, max(theLevels$water)*1.1, TR/25)#/100  #cm, MSL
  B= seq(0,50000, length.out = length(z))
  SSC.i1 = 1  #SSC, mg/l (g/m3)
  
  run_allSSC = mineralBiomass(theLevels)
  inputs = expand.grid(z=z, SSC.i=SSC.i1, B=B) 
  inputs2 = lapply(seq_len(nrow(inputs)), function(i) inputs[i,])
  
  ##
  require(foreach)
  require(doParallel)
  require(pracma)
  
  
  ###############################################
  
  sfInit(parallel=T, cpus=detectCores()-1)
  sfExportAll()
  results = sfLapply(inputs2, fun= run_allSSC)
  sfStop()
    
    output = cbind(inputs,unlist(results))
   
 # plot(z,output[output$B==500,]$`unlist(results)`*10000*20, ylab='Mineral Dep', type='l')
  
 # plot(z,output30[output30$B==500,]$`unlist(results30)`*10000, ylab='Mineral Dep', type='l')
  
  out =NULL
  for(j in 1:length(z)){
   out=cbind(out, output[output$z==z[j],]$`unlist(results)`)
  }
  
  colnames(out)= z
  rownames(out)=B
 # filled.contour(out)
 
  interp1 = function(z1, B1){
    res=pracma::interp2(x=z, y=B, Z=out, xp=z1/100, yp=B1, method='linear')
    if(is.na(res) & z1/100<min(z)){
      res=pracma::interp2(x=z, y=B, Z=out, xp=min(z)/100, yp=B1, method='linear')
    }
    if(is.na(res) & z1/100>max(z)){
      res=0
    }
    if(is.na(res) & B1>max(B)){
      res=interp2(x=z, y=B, Z=out, xp=z1/100, yp=max(B), method='linear')
    }
    return(res)
  }
 # interp1 = interp2(x=z, y=B, Z=out, xp=0.6, yp=1000, method='linear')*5*10000
 #interp1(-50, 500) 

  return(interp1)  #returns a function that provides deposition in g/cm2 provided input elevation (cm) and aboveground biomass (g/m2)
  
}


