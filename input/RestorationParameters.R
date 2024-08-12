
bulk.densityMass = function(perOM, Mass){  #morris et al 2016
  k1=0.085
  k2=1.99
  min(1/(perOM/k1 + (1-perOM)/k2),1.9)
}


#organic and mineral densities
OMden=1.14 #g/cm3
MINden=2.65 #g/cm3
c.density=0.42  #proportion of biomass that is carbon

HC = read.csv('input/RichmondBr_TidalConst.csv')
MSL = 98.2 #cm  redwood city
msl =0 # MSL/100 #m
TR = 2 #m  diurnal tidal range
historicSLR = 0.235 #cm/yr, eg, NOAA tides and currents

cores = rbind( cbind(read.csv('input/CC_Callaway_AHigh.csv'), CoreID='AHigh'),
               cbind(read.csv('input/CC_Callaway_BHigh.csv'),CoreID='BHigh'),
               cbind(read.csv('input/CC_Callaway_AMid.csv'),CoreID='AMid'),
               cbind(read.csv('input/CC_Callaway_BMid.csv'),CoreID='BMid'),
               cbind(read.csv('input/CC_Callaway_ALow.csv'),CoreID='ALow'),
               cbind(read.csv('input/CC_Callaway_BLow.csv'), CoreID='BLow'))

sitePars = list(
  CC =1,                #maximum carrying capacity. Should be 1 for marsh species. 
  spILim= 0.6,          #flooding time [0-1] limit for species
  perOMCohorts = .15,   #the initial proportion of organic material in the soil cohorts
  surfOMDep = 0.1,     #proportion of mineral deposition that is organic material. Gets added to refractory organic matter pool.  
  wood2soil = 0.01,     #rate that wood litter is integrated into soil. should be slow [~0.01]
  kdec=0.001,          #daily decomposition decay constant. Strong control over organic accumulation rate
  decompSens=1.001,     #controls the sensitivity of decomposition to inundation - values >1 allow decomposition to occur when fully flooded [1-1.5]. 
  maxPorosity = 0.84,   #porosity at the surface [from soil cores]
  minPorosity=0.73,      #porosity at depth [from soil cores]
  porosityRate = 2,  #rate that porosity changes as a function of mass of above material. Controls accretion rate
  midPtDepth=0,          #g, controls the shape of the porosity function. higher value will slow the rate of consolidation
  maxDeltaPor = 0.005   #maximum change in cohort porosity per year. Controls soil consolidation rate 
)

speciesList =  c( 'SalPac', 'SpaFol') 
nSpecies = length(speciesList)
years2CC = rep(c(15,5), nSpecies) # not used in current iteration of species growth



#From soil cores -- China  Camp. Callaway et al 2012
minCal = list(minAcc =  c(327.00, 527.7824,  1165.045, 183.4084, 1810.126)/10000, #g/cm2
              elevMSL = ( c(186.8, 181.5,  180, 179, 195.5) - MSL),  # NAVD88 elevation cm, made relative to MSL
              omAcc =c(115.1619, 121.6546,  163.69, 50.89028, 145.1573)/10000, #g/cm2
              acc = c(.13, .22, .29, .1, .24)*10)  #mm/yr

