
#Functions from marsh organ experiments in San Francisco Bay. janousek et al 2016. doi: 10.3354/meps11683

#Bolbocenous maritimas
boma.B = function(inundation.time){
  #inundation.time= inundation.time*100
  if(inundation.time<30){ #5th percentile inundation from veg surveys at Petaluma
    return( max(0, -0.0022	*inundation.time^2+	0.1210*inundation.time +	4.13))
  } else {
    return(0)
  }
  
}

#Salicornia pacifica
sapa.B = function(inundation.time){
  if(inundation.time<=50){  # 22 min inundation of SpaFol at Petaluma-- need distributions to connect
    0.0036 * inundation.time^2 -	0.4629*inundation.time	+	14.49
  } else {
    return(0)
  }
}

#Juncus balticus
juba.B = function(inundation.time){
  max(0, -0.02807	*inundation.time^2+	1.10969*inundation.time	+	18.5437 )
}

#Scheonplectus americanus
scam.B = function(inundation.time){
  max(0, -0.0802*inundation.time^2+	3.5626*inundation.time	+	40.203)
}

#Spartina foliosa
spfo.B = function(inundation.time){
  if(inundation.time>1){  #10
    return( max(0, -0.00798*inundation.time^2+	0.60243*inundation.time	+	1.08638))
  }  else {
    return(0)
  }
}
#Scheonoplectus actus
scac.B = function(inundation.time){
  #  -0.0149*inundation.time + 1.4842  #SCAC only, Lisa Schile
  max(0, -0.028*inundation.time^2 + 1.8289*inundation.time + 28.289 ) #SCAC from Lisa Schile, and highest elevation SCAM from CJN
}

scacFun= approxfun(0:100/100, sapply(0:100, scac.B)/max(sapply(0:100, scac.B)), rule=2)
spfoFun= approxfun(0:100/100, sapply(0:100, spfo.B)/max(sapply(0:100, spfo.B)), rule=2)
scamFun= approxfun(0:100/100, sapply(0:100, scam.B)/max(sapply(0:100, scam.B)), rule=2)
jubaFun= approxfun(0:100/100, sapply(0:100, juba.B)/max(sapply(0:100, juba.B)), rule=2)
sapaFun= approxfun(0:100/100, sapply(0:100, sapa.B)/max(sapply(0:100, sapa.B)), rule=2)
bomaFun= approxfun(0:100/100, sapply(0:100, boma.B)/max(sapply(0:100, boma.B)), rule=2)

saveRDS(scacFun, 'input/speciesFuns/Biomass/ScheonoplectusActus_relProd_inundTime.rds')
saveRDS(spfoFun, 'input/speciesFuns/Biomass/SpartinaFoliosa_relProd_inundTime.rds')
saveRDS(scamFun, 'input/speciesFuns/Biomass/ScheonoplectusAmericanus_relProd_inundTime.rds')
saveRDS(jubaFun, 'input/speciesFuns/Biomass/JuncusBalticus_relProd_inundTime.rds')
saveRDS(sapaFun, 'input/speciesFuns/Biomass/SalicorniaPacifica_relProd_inundTime.rds')
saveRDS(bomaFun, 'input/speciesFuns/Biomass/BolboschoenusMaritimus_relProd_inundTime.rds')





