

#Morris et al 2013 SpaAlt
#2008

spaAlt_Prod = function(x) {
  out= (-5.3*x+1856)/1856  #x is cm NAVD88. Estimated y-intercept from figure
  ifelse(out<0,0,ifelse(out>1,1,out))

}

elevMSL= seq(-50, 200, 1)
inund=sapply(elevMSL, inundFun)
plot(elevMSL, inund)

MSL= -13

plot(inund, sapply(elevMSL+MSL, spaAlt_Prod))

spaAlt_Prob_inund = approxfun(inund, sapply(elevMSL+MSL, spaAlt_Prod), rule=2)

spaAlt_Prob_inund(0.01)

saveRDS(spaAlt_Prob_inund,'input/speciesFuns/Biomass/SpartinaAlterniflora_relProd_inundTime.rds' )

plot(elevMSL, sapply(elevMSL, spaAlt_Prod))
