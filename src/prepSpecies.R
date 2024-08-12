
prepFlooding = function(floodFun){
  flooding = seq(0,100, 0.1)
  relFlood =pmax(sapply(flooding, floodFun)/max(sapply(flooding, floodFun)),0)
  # relFlood[(which(relFlood==0))[1]:length(relFlood)] = 0
  approxfun(flooding,relFlood, rule=2) 
}


prepSpecies = function(speciesList, years2CC, inundFun){
  
  nicheFUN=list()
  rootDepthFUN=list()
  growth=NA
  cov2AGB=list()
  rootShoot=list()
  pars=list()
  optRSFUN=list()
  
  typeList=NULL
  for(i in 1:length(speciesList)){
    typeList = c(typeList, get(speciesList[i])$Type)
  }
  
  mixedComm= ifelse('Tree' %in% typeList & 'Grass' %in% typeList, 1, 0)
  
  floodTime = sapply(seq(-TR/2,TR/2,0.01)*100, inundFun) #flooding time based on elevation
  
  for(i in 1:length(speciesList)){
    sp= get(speciesList[i])
    # if(sp$Type=='Grass'){
    nicheFUN[[i]] = (sp$nicheFUN)    #prepFlooding
    #if(is.function(sp$root_shoot)){
    rootShoot[[i]] = sp$root_shoot
    #} else {
    
    z1=seq(0,100,0.1)
    optFlood= inundFun(z1[ which.max(sapply(z1, nicheFUN[[i]]))])
    optFlood = ifelse(optFlood<0, 0, optFlood)
    optRSFUN[[i]] = approxfun(floodTime, (1+sp$RSslp*(optFlood-floodTime)), rule=2, ties=mean)
    
    if(mixedComm==1){
      if(sp$Type=='Grass' | sp$Type=='Seagrass'){
        growth[i] = sp$absGrowth#*sitePars$CC#*sp$compA
      } else {
        growth[i] = sp$absGrowth#/sitePars$CC
      }
    } else {
      if(sp$Type=='Seagrass'){
        growth[i] = sp$absGrowth*sitePars$CC#*sp$compA
      } else {
        growth[i] =sp$absGrowth
      }}
    
    cov2AGB[[i]] = sp$biomassFUN
    rootDepthFUN[[i]]=sp$r.depthFUN
    
    if(i==1){
      pars=sp
    } else {
      pars = combineList(pars, sp)
    }
    
  }
  
  nicheFUN[(length(nicheFUN)+1):(length(nicheFUN)*2)] = nicheFUN
  
  #Calculate growth rate based on years2CC and absolute growth rate. 
  #Works best for trees where growth is defined as a change in basal area 
  
  #growth2= 1+(growth-mean(growth))/mean(growth)  #normalized growth rate
  #
  #rate= rep(NA, length(speciesList))
  # for(j in 1:length(speciesList)){
  #   #find optimum rate for Xyears
  #   CC= 1#0.02 # any starting carrying capacity
  #   Cstart=.01*CC
  #   Abs_rate=0.01
  #   Xyears= years2CC[j] #years to 99% CC
  #   x=seq(1,2000)
  #   k=0
  #   test=seq(.02,2,.001)
  #   Y_CC=array(0,length(test))
  #   for( Abs_rate in test){
  #     k=k+1
  #     C=array(Cstart,length(x))
  #     for( i in 2:length(x)){
  #       C[i]=C[i-1]+C[i-1]*Abs_rate*(1-C[i-1]/CC)
  #     }
  #     Y_CC[k]=approx(C/CC,x,.99)$y
  #   }
  #   OUT=rbind(test,Y_CC)
  #   
  #   rate[j]=approx(Y_CC,test,Xyears)$y  #~.1 for 90 years to 99% of CC
  # }
  #
  #pars$r1 = growth2*rate
  #pars$r1 = pars$absGrowth
  
  
  pars$aliveIndex = 1:length(speciesList)
  pars$deadIndex=(length(speciesList)+1):(length(speciesList)*2)
  pars$mixedComm = mixedComm
  pars$grassI= rep(0, length(speciesList)*2)
  pars$grassI[which(c(pars$Type,pars$Type)=='Grass')]=1
  pars$treeI = which(c(pars$Type,pars$Type)=='Tree')
  pars$seagrassI = which(c(pars$Type, pars$Type)=='Seagrass')
  
  # pars$CC=1   #for marsh use CC=1. For forest, use fraction basal area (95th percentile) from inventory plots.
  
  return(list(nicheFUN=nicheFUN, coverAgTransfer=cov2AGB, pars=pars, optRSFUN=optRSFUN, rootDepthFUN=rootDepthFUN, rsFUN=rootShoot))
}
