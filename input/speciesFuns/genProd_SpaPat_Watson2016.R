

#Watson et al 2016 SpaPatens
#TR site

spaPat_Prod = function(x) {
  out= ( -0.35*x + 26.18)/26.18  #x is inundation time
  ifelse(out<0,0,ifelse(out>1,1,out))

}

plot(seq(0,100,1), sapply(seq(0,100,1),spaPat_Prod))

spaPat_Prob_inund = approxfun(seq(0,1,0.01),  sapply(seq(0,100,1),spaPat_Prod), rule=2)

saveRDS(spaPat_Prob_inund,'input/SpeciesFuns/Biomass/SpartinaPatens_relProd_inundTime.rds' )


