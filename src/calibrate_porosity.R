#Calibration soil porosity

#calibrate porosity function to soil core BD and %OM


coreIDs = unique(cores$core_id)
if(is.null(coreIDs)){
  print('missing soil core IDs. Check input file')
}

cores2=NULL
for(i in 1:length(coreIDs)){
  coreTemp = cores[cores$core_id==coreIDs[i],]
  if(!is.na(sum(coreTemp$OM)) & !is.na(sum(coreTemp$BD)) & nrow(coreTemp)>1){
    int = diff(coreTemp$Depth)
    int=ifelse(int==3, 4, int)
    coreTemp$int = c(int[1], int)
    coreTemp$Mass = coreTemp$BD*coreTemp$int
    coreTemp$MassAbv = c( cumsum(coreTemp$Mass))
    coreTemp$MassAbv = c(0, coreTemp$MassAbv[1:(nrow(coreTemp)-1)])
    coreTemp$DensAbv = coreTemp$MassAbv/cumsum(coreTemp$int)
    cores2=rbind(cores2, coreTemp)
  }
}
cores= cores2

surf=cores[cores$Depth==0,]
OMden=0.5 #initial organic matter density
POR = 1- (surf$BD*(surf$OM/100/OMden)+(1-surf$OM/100)/MINden)/(surf$OM/100*OMden+(1-surf$OM/100)*MINden)
max(POR)

calPor=function(pars, cores2, maxPor){
  cores2=cores2[cores$Depth>3,]
  pOM=cores2$OM/100
  OMden=pars[3]
  cores2$PORcalc = porosity(cores2$DensAbv, pars[1], pars[2], maxPor, pOM=pOM)
  cores2$BD2 = (cores2$Mass)/((cores2$Mass*cores2$OM/100/OMden+ cores2$Mass*(1-cores2$OM/100)/MINden+cores2$PORcalc*(cores2$Mass*cores2$OM/100/OMden+ cores2$Mass*(1-cores2$OM/100)/MINden)/(1-cores2$PORcalc)))
  cores2 = cores2[complete.cases(cores2$BD2),]
  mean(sqrt((cores2$BD-cores2$BD2)^2))
}

calPorSurf = function(pars, cores2){
  cores2 = cores2[cores$Depth<=2,]
  pOM=cores2$OM/100
  cores2$PORcalc = porosity(cores2$DensAbv,4, 0.6,pars, pOM=pOM)
  cores2$BD2 = (cores2$Mass)/((cores2$Mass*cores2$OM/100/OMden+ cores2$Mass*(1-cores2$OM/100)/MINden+cores2$PORcalc*(cores2$Mass*cores2$OM/100/OMden+ cores2$Mass*(1-cores2$OM/100)/MINden)/(1-cores2$PORcalc)))
  cores2 = cores2[complete.cases(cores2$BD2),]
  mean(sqrt((cores2$BD-cores2$BD2)^2))
}

#calibrate surface porosity
outSurf=optim(par=c(0.9),fn=calPorSurf, cores2=cores, method='L-BFGS-B',lower=c(0.8), upper=c(0.98))
outSurf$par

#calibration rate and minimum porosity. Max porosity a constant
out=optim(par=c(4, 0.6, 0.5),fn=calPor, cores2=cores, maxPor=outSurf$par, method='L-BFGS-B',lower=c(0.05, 0.08, 0.15), upper=c(10, 0.85, 1.3))

porosityRate=out$par[1]
minPorosity=out$par[2]
OMden=out$par[3]
maxPorosity=outSurf$par

pOM=cores$OM/100
cores$PORcalc = porosity(cores$DensAbv, porosityRate, minPorosity, maxPorosity, pOM)
cores$BD2 = (cores$Mass)/((cores$Mass*cores$OM/100/OMden+ cores$Mass*(1-cores$OM/100)/MINden+cores$PORcalc*(cores$Mass*(cores$OM/100)/OMden+ cores$Mass*(1-cores$OM/100)/MINden)/(1-cores$PORcalc)))

bd.rmse=rmserr(cores$BD, cores$BD2)$rmse #local calibration
rmserr(cores$BD, sapply(pOM, bulk.density))$rmse #Ideal mixing model

par(mfrow=c(1,2))
plot(cores$BD, cores$BD2, ylim=c(0,1), xlim=c(0,1), xlab='Observed bulk density', ylab='Modeled bulk density', main=paste0('RMSE= ',round(bd.rmse,3)))
abline(a=0, b=1)

#points(cores$BD, sapply(pOM, bulk.density), col='blue')

plot(cores$Depth, cores$BD, ylim=c(0,1.5), pch=16, ylab='Bulk density g/cm3', xlab='Depth (cm)')
points(cores$Depth, cores$BD2, col='red', pch=16)
#points(cores$Depth, sapply(pOM, bulk.density), col='blue')
