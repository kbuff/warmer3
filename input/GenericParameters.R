# WARMER input
speciesList =  c( 'mangrove', 'marsh','upland') 
nSpecies = length(speciesList)
years2CC = c(25,10,60)

MSL = 80
mangroveAgeThres = 15 #age when young mangroves mature to old mangroves

#values from root analysis
calVals=list(
  deadTotalRootmin=0.54,  #ratio of dead:total root mass -- not structural
  deadTotalRootmax=0.6,   #ratio of dead:total root mass -- not structural
  totRootMassmin=0.238,  #g/cm2 -- not structural
  totRootMassmax=0.679,  #g/cm2 -- not structural  
  liveRootsmin=0.15,     #g/cm2 -- not structural
  liveRootsmax=0.35,     #g/cm2 -- not structural
  littermin = 2.82,      #Mg/ha/yr
  littermax=3.26         #Mg/ha/yr  
)


ssc=10  #calibration, for inundation function mg/L
#SSC= 5# 21.039
#SSC_sd= 15.53

t_slr = list( ###???
  target_slr=c(53, 115, 183), #2100
  yr30 = c(15, 22, 29),
  yr50 = c(27, 42,58) )


#initial_org = raster("input/Sites/Ding/LEAN_DEM_DingDarling_11_8_21.tif")*100

# SLRC.input
slrc.input <- list(
  slr.target = 142,
  historical.slr =0.193   #Lovelock et al. 2015 -- 1984 to 2010 https://doi.org/10.1007/s10021-011-9443-9
)

# OMARC input #not used as input, only for the calibration (dotted line on graphs)
omarc.input <- list( ## SALTMARSH -- Santini et al. 2019 Table S1. https://doi.org/10.1007/s10021-019-00373-x
 core_om  = c(273.9925, 340.9086, 149.0098)/10000, ## g/cm2 #J Drexler cores #DAB, DBB, DCB 
 elev_om     =c(96.4, 114.4, 104.5)-MSL, #from Geosciences 5m DEM, Nadia/Cath may have exact data but unlikely. #June 2021 elev surveys    #c(23.3, 22.19, 20.9) - MSL, # cm (LEAN DEM)
	acc = c(0.34, 0.22, 0.18),   #accretion rate cm/yr (Pb, 50 yrs) 
# acc = c(0.46, 0.44, 0.36), #accretion rate cm/yr (Pb, 100 yrs)  
  zTolerance  = 1,
  granularity = 8#,
 # file.name   = 'input/Alohkapw/omarc'
)

c.density=0.46 #average of mang/mel/cas is 0.48 (Robyn, Adame, Kelleway), sm 0.38-0.42
kdec= 0.002
AGE=F  #track age of trees?  no needed
dynRS = F #dynamic root:shoot ratio, based on optimum elevation

# MSARC input
msarc.input <- list(   ##use subtropical Aus dated cores 
#	harmonic.constituents =  read.csv(paste('input\\Naples_TidalHarmonics_m_MSL.csv',sep=''), header=T),
 core_sed = c(200.176, 177.38, 99.084)/10000,  #J Drexler cores #DAB, DBB, DCB 
  elev_sed  = c(96.4, 114.4, 104.5)-MSL #June 2021 elev surveys    #c(23.3, 22.19, 20.9) - MSL # cm (LEAN DEM)
#  file.name   = 'input_Pohnpei/Sites/Alohkapw/msarc',
#  minDep.file.name = 'input/Mineral_Deposition_5_SSC_ws02.csv'
)


#Overall parameters for the study site
sitePars = list(
  CC =0.0036,# 0.0036, #from ding darling
  #ALPH = 0.01, 
 # k.labile =2,   #   controls rate of loss with inundation freq.  HIGHER = slower decomp & more LOM.  Largely responsible for attractive domain. Lower vals=more separation between elevations
 # k.labile.mx = 0.5, #0.8  0.9,  #controls maximum rate of loss per yr.
#k.labile = 15, #controls rate of loss with root density. higher values mean more decomposition
# k.labile.mx = 0.05, #minimum rate of decomposition
#   k.refrac.mx.coef = 100, #125
  initCover = 1e-5,
  spILim= 0.6  #inundation time limit for all species
)



mangrove = list(
  Type='Tree',
  floodFun =  species1_nicheFun,  #niche curve based on %time inundated
  #maxBiomass= function(cov,t) cov*(101/0.42)*exp(-29.6/t)/100,   ##blueCAM model-- includes age
#  maxBiomass= function(cov) ifelse(cov==0, 0, exp(5.919098 + log(cov+1e-10)*0.969934) ) , #cover only, no age      #   #for all 3 species combined (4/6/2023)   #agb - cover function
 
  maxBiomass = function(cov)     ifelse(cov<=0, 0,exp(6.48842 + log(cov+1e-10)*1.14705)  ), #g/cm2 #Brugeria
   #absGrowth = function(t) 0.0001*exp(-0.109*t), ## Xiong et al. global review, plantations, up to 25 years old stands
  absGrowth = 0.00018, ## non-plantations Moreton Bay avergae #0.00035,  #m2/yr
  strucRoots= 0.8,    #very large proportion of root biomass -- control over observerable live root biomass 
  coarseRoots= 0.71,  # 5-20 mm s --of remainng (1-strucRoots)
  smallRoots = 0.16, #2-5 mm  --of remainng (1-strucRoots)
  fineRoots = 0.13,   #<2 mm.   --of remainng (1-strucRoots)
  deadFallRate = 0.2,# 0.2,   # rate that standing dead fall   
  r.roots = 85, #function(cov) 7.56*cov+0.5, #defines max rooting depth, cm  ##Ohira et al. 2013 https://link.springer.com/article/10.1007/s00468-012-0782-8
  root_shoot = 0.58, #  calibrated to living root biomass range- Conrad 2022.    
  root_por = 0.27,   #Cheng et al 2012
  r.depthFUN = function(mxdepth, t) 10+mxdepth*exp(-10/t) , #rooting depth increases with age.
  kma = 0.1, # MAOM production rate-function of living root mass [not structural] -- calibration
  kr = 0.2, #0.43, #0.1,    #fine root turnover -
  km = 0.1, #0.08,      #coarse root turnover   - 
  kt = 0.05, #0.015,     #very large root turnover (structural, cable)  -- Castaneda-Moya et al 2011, same at  both sites
  ks = 0.16, # small root turnover
  krr = 0.2, #fine root rate to particulate OM  
  kmm = 0.05,  #coarse root rate to particulate OM
  ktt= 0.02,  #very large root rate to particulate OM
  fc1 = 0.8,  #0.75          #labile fraction of fine +small root
  fc2= 0.86,    #0.85         #labile fraction of coarse + structural root
  leafProp = 0.89,
  leafRefrac = 0.15,     #leaf refractory fraction 
  litterDep = 0.005,   # 0.015, #0.15,      #litter deposition rate
  litterExport = 0.25,#.8,      #maximum proportion of litter that remains on-site ### Mean litter export during tidal inundation is approximately 50% of the mean C litterfall observed in mangroves (Saenger & Snedaker 1993)
 # r1 =0.136,# 0.205,            #growth rate (basal area/yr)
  seeding = 0.1,           #initial seeding rate
  compA = 2,             #scalar for space competition 
  wood.den=0.725,         # g/cm3 A. marina; Mackey 1993 https://www-publish-csiro-au.ezproxy.library.uq.edu.au/mf/pdf/MF9930721 
  root.den = 0.24,        #root density, g/cm3. Derived from Fig 3 Ola et al. 2018 -- 0.04/0.17 = 0.235
  initCover = 0.02,  #% CC
  SRslp=-2, #4
  mann= 0.03 # made up
  
)



#uncomment to include young mangroves

# if('mangrove' %in% speciesList){
#   mangroveY = mangrove
#   speciesList = c('mangroveY', speciesList)
#   mangroveY$absGrowth = mangroveY$absGrowth*1
#   mangrove$initCover=0                #old mangroves cant initialize
#   mangroveY$deadFallRate=0.9  #young mangrove fall fast
#   mangroveY$strucRoots=0.2
#   nSpecies = length(speciesList)
#   years2CC = c(years2CC[which('mangrove' %in% speciesList)], years2CC)
#   }



# r.depthFUN = function(mxdepth, t) -10+(mxdepth+10)*exp(-1.5/t)
# plot(1:100, sapply(1:100,r.depthFUN, mxdepth=300), xlab='Age', ylab='Rooting Depth', type='l')
# 
# 
# r.depthFUN1 = function(cov) 7.56*cov+0.5
# plot(seq(0,10,0.1), sapply(seq(0,10,0.1),r.depthFUN1), xlab='DBH', ylab='Rooting Depth')


# maxBiomass= function(cov) cov[,1]*101*exp(-29.6/cov[,2])/100
# 
# cov=cbind(c(0.1, 0.25, .5, 1), c(2, 5, 10, 50))
# 
# maxBiomass(cov)

#plot(seq(0,1,0.01),sapply(seq(0,1,0.01), mangrove$maxBiomass), xlab='Cover', ylab='Aboveground Biomass g/cm2', type='l')

# maxBiomass = function(cov, t)   cov*(0.75)
# 
# mangAgeB = function(cov, t) 101*exp(-29.6/t) #*(350/10000)
# mangAgeBCov=  function(cov, t) cov*101*exp(-29.6/t) #*(350/10000)


# plot(seq(0,1,0.01),sapply(seq(0,1,0.01), maxBiomass), xlab='Cover', ylab='Aboveground Biomass', type='l')
# 
# plot(seq(0,100,1),sapply(seq(0,100,1), mangAgeBCov, cov=1)/100, xlab='Age', ylab='Aboveground Biomass', type='l', col='red')
# lines(seq(0,100,1),sapply(seq(0,100,1), mangAgeBCov, cov=0.5)/100, xlab='Age', ylab='Aboveground Biomass', type='l', col='blue')
# lines(seq(0,100,1),sapply(seq(0,100,1), mangAgeBCov, cov=0.75)/100, xlab='Age', ylab='Aboveground Biomass', type='l', col='cyan')
# 
# plot(seq(0,1,0.01),sapply(seq(0,1,0.01), mangAgeBCov, t=100)/100, xlab='Cover', ylab='Aboveground Biomass', type='l', col='black')
# lines(seq(0,1,0.01),sapply(seq(0,1,0.01), mangAgeBCov, t=25)/100, xlab='Age', ylab='Aboveground Biomass', type='l', col='blue')
# lines(seq(0,1,0.01),sapply(seq(0,1,0.01), mangAgeBCov, t=50)/100, xlab='Age', ylab='Aboveground Biomass', type='l', col='red')
# lines(seq(0,1,0.01),sapply(seq(0,1,0.01), mangAgeBCov, t=75)/100, xlab='Age', ylab='Aboveground Biomass', type='l', col='cyan')
# lines(seq(0,1,0.01),sapply(seq(0,1,0.01), mangAgeBCov, t=5)/100, xlab='Age', ylab='Aboveground Biomass', type='l', col='gold')



marsh = list(
  Type='Grass',
  floodFun =  species2_nicheFun,  #based on basin WL
 # maxBiomass = function(cov,t) cov*(1.36/0.42)*exp(-1/t)/100 ,  #blueCAM a, but k not stated, copy mangroves   #maximum aboveground biomass
   maxBiomass= function(cov) ifelse(cov==0, 0, exp(3 + log(cov+1e-10)*0.9) ) , #cover only, no age      #   #for all 3 species combined (4/6/2023)   #agb - cover function
  
  absGrowth = 0.0002,  #m2/yr, ## 6.23ha =62300 m2 grew in 3 yrs (BlueCAM salt pond case study but report pulled??), assume stem density
  strucRoots= 0,    #very large proportion of root biomass
  coarseRoots= 0,  # 5-20 mm s
  smallRoots = 0.2, #2-5 mm
  fineRoots = 0.8,   #<2 mm. 
  deadFallRate = 0.8,# 0.2,   #* rate that standing dead fall   
  r.roots = 30,         #defines max rooting depth, cm
  root_shoot = 0.179 ,# for pickleweed; Janousek 2016   
  root_por = 0.15, 
  r.depthFUN = function(mxdepth, t) mxdepth ,
  kma = 0.05, # MAOM production rate-function of living root mass [not structural]
  kr = 0.5, #0.43, #0.1,    #fine root turnover - Conrad 2022 Table 2.2
  #kWHATALPHABET = 0.2,  #small root turnover
  km = 0.2, #0.08,      #coarse root turnover   - Conrad 2022 Table 2.2
  kt = 0, #0.015,     #very large root turnover (structural, cable)  -- Castenya-Moya et al 2011, same at  both sites
  ks = 0.3, #small root turnover
  krr = 0.3, #fine root rate to particulate OM  
  kmm = 0.15,  #coarse root rate to particulate OM
  ktt= 0.1,  #very large root rate to particulate OM
   fc1 = 0.8,  #0.75          #labile fraction of fine root
  fc2= 0.86,    #0.85         #labile fraction of coarse root
  leafProp = 1,
  leafRefrac = 0.15,     #leaf refractory fraction 
  litterDep =0.02,# 0.015, #0.15,    ## kg/m2  #litter deposition rate ##Kelleway et. al 2017 https://bg.copernicus.org/articles/14/3763/2017/bg-14-3763-2017.pdf 
  litterExport = 0.6,#.8,      #maximum proportion of litter that remains on-site
  #r1 =0.136,# 0.205,            #growth rate (basal area/yr)
  seeding = 0.1,           #initial seeding rate
  compA = 1,             #scalor for space competition 
  wood.den=0.2,         #
  root.den = 0.2,        #root density, g/cm3
  initCover = 0.02,  #% CC
  SRslp=0, #4
  mann= 0.03 # made up
  
)

# plot(seq(0,1,0.01), sapply(seq(0,1,0.01), marsh$maxBiomass, t=1)*10000, col="black", type='l', ylim=
#        c(0,500))
# lines(seq(0,1,0.01), sapply(seq(0,1,0.01), marsh$maxBiomass, t=5)*10000, col="blue")
# lines(seq(0,1,0.01), sapply(seq(0,1,0.01), marsh$maxBiomass, t=10)*10000, col="red")
# lines(seq(0,1,0.01), sapply(seq(0,1,0.01), marsh$maxBiomass, t=20)*10000, col="green")

upland = list(Type='Tree',
              floodFun=species3_nicheFun,
             # maxBiomass=  function(cov,t) cov*(100/0.42)*exp(-29.6/t)/100, #g/cm2    #~10000 g/m2 at full CC
              maxBiomass= function(cov) ifelse(cov==0, 0, exp(5 + log(cov+1e-10)*0.7) ) ,
              absGrowth = 0.0027,  #m2/yr
              strucRoots= 0.6,    #very large proportion of root biomass
              coarseRoots= 0.56,  # 5-20 mm s
              smallRoots = 0.24, #2-5 mm
              fineRoots = 0.2,   #<2 mm. fine-coarse root proportion of remaining root bioass (1-vBigRoots)
              deadFallRate = 0.15,  #rate that standing dead fall
              r.roots = 300, # cm 0.1,         #defines rooting depth 
              root_shoot = 0.378,  #Shoot:Root ratio; Qu et al 2019 https://doi.org/10.1016/j.gecco.2019.e00606
              root_por = 0.1,
              r.depthFUN= function(mxdepth, t) -10+(mxdepth+10)*exp(-1.5/t),
              kma = 0.05, # MAOM production rate-function of living root mass [not structural]
              kr = 0.2,             #fine root turnover
              km = 0.1,             #coarse root turnover
              kt = 0.02,             #very large root turnover (structural, cable)
              ks = 0.15,   #small root turnover
              krr = 0.25, #fine root rate to particulate OM  
              kmm = 0.15,  #coarse root rate to particulate OM
              ktt= 0.1,  #very large root rate to particulate OM
              fc1 =0.8,             #labile fraction of fine root
              fc2= 0.86,             #labile fraction of coarse root
              leafProp = 0.89,
              leafRefrac = 0.2,     #leaf refractory fraction
              litterDep = 0.2 ,      #litter deposition rate
              litterExport = 0.6,      #maximum proportion of litter that remains on-site
              # r1 = 0.08536,            #growth rate (basal area/yr)
              seeding = 0.1,           #initial seeding rate
              compA = 1,             #scalor for space competition-- lower=more competitive
              wood.den= 0.75,         #wood density, g/cm3; https://www.woodsolutions.com.au/wood-species/hardwood/tea-tree-broad-leaved 
              root.den = 0.2,         #root density, g/cm3
              initCover=0.02,
              SRslp= 0,
              mann = 0.03   #made up
)







