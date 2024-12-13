#input

#read in parameter file
genPars = openxlsx::read.xlsx( paste0('input/', site.name, '/', site.name, '_Parameters.xlsx'), 'Site')
spPars =  openxlsx::read.xlsx( paste0('input/', site.name, '/', site.name, '_Parameters.xlsx'), 'Species')

#parse the parameter file, assigning to appropriate list
genParsGlobal = genPars[genPars$List == 'global', ]
genParsSite = genPars[genPars$List == 'sitePars', ]
genParsPlot = genPars[genPars$List == 'plotting', ]
for (i in 1:nrow(genParsGlobal)) {
  assign(genParsGlobal$Parameters[i], genParsGlobal$Value[i])
}

sitePars = list()
for (i in 1:nrow(genParsSite)) {
  sitePars[[genParsSite$Parameters[i]]] = genParsSite$Value[i]
}

calVals = list()
for (i in 1:nrow(genParsPlot)) {
  calVals[[genParsPlot$Parameters[i]]] = genParsPlot$Value[i]
}

bulk.density = function(perOM) {
  #morris et al 2016
  min(1 / (perOM / k1 + (1 - perOM) / k2), 1.9)
}

#calibration soil cores
cores = rbind( cbind(read.csv('input/ChinaCamp/CC_Callaway_AHigh.csv'), core_id='AHigh'),
               cbind(read.csv('input/ChinaCamp/CC_Callaway_BHigh.csv'),core_id='BHigh'),
               cbind(read.csv('input/ChinaCamp/CC_Callaway_AMid.csv'),core_id='AMid'),
               cbind(read.csv('input/ChinaCamp/CC_Callaway_BMid.csv'),core_id='BMid'),
               cbind(read.csv('input/ChinaCamp/CC_Callaway_ALow.csv'),core_id='ALow'),
               cbind(read.csv('input/ChinaCamp/CC_Callaway_BLow.csv'), core_id='BLow'))

coreRates = read.csv(paste0('input/', site.name, '/', site.name, '_soil_core_summary.csv'))
minCal = list()
minCal$minAcc = coreRates$MAccum_gm2yr / 10000
minCal$elevMSL = coreRates$elevation * 100 - MSL
minCal$omAcc = coreRates$OAccum_gm2yr / 10000
minCal$acc = coreRates$Accretion_cmyr

#tidal harmonics
HC = read.csv(paste0("input/", site.name, '/', site.name, '_tidalconstituents.csv'),
              header = T)

#combine species parameters to 'inputPars' list
spIndex = which(grepl('sp.', colnames(spPars)))
speciesList = sapply(
  colnames(spPars)[spIndex],
  FUN = function(x)
    strsplit(x, 'sp.')[[1]][2]
)
nSpecies = length(speciesList)

inputPars = list()
for (i in 1:nrow(spPars)) {
  #assign(genPars$Parameters[i], genPars$Value[i])
  if (spPars$Type[i] == 'numeric') {
    inputPars[[spPars$Parameters[i]]] = as.numeric(spPars[i, spIndex])
  } else {
    inputPars[[spPars$Parameters[i]]] = (spPars[i, spIndex])
  }
}

theRow=which(spPars$Parameters=='spBiomassFUNpath')
for (i in 1:nSpecies) {
  if (!is.na(spPars[theRow,spIndex[i]])) {
    biomass_fun = readRDS(as.character(spPars[theRow,spIndex[i]]))
    inputPars$biomassFUN[[i]] <- local({
      species_fun <- biomass_fun  # set the biomass function for this species
      function(inundation, cc, maxB) {
        cc * (maxB / 10000) * species_fun(inundation)
      }
    })
  } else { #niche function is also used as biomass function
    inputPars$biomassFUN[[i]] = function(inundation, cc, maxB)
      cc * (maxB / 10000)
  }
}


#Generate root:shoot ratio function relative to inundation time [0-1]
#SalPac RS
SAPA_root_shoot100 = function(i100) exp(min(-0.22, -1.5649- i100*0.0005929 + 0.0001965*i100^2))# Janousek et al 2016, function takes in inundation duration % [0-100%]
sapaRS = approxfun( (0:100)/100, sapply(0:100, SAPA_root_shoot100)/max( sapply(0:100, SAPA_root_shoot100)), rule=2 )

#SpaFol RS
SPFO_root_shoot100=function(i100) exp(min(1,(0.5149- i100*0.03892 + 0.0005896*i100^2)))  # Janousek et al 2016, function takes in inundation duration % [0-100%]
spfoRS =  approxfun( (0:100)/100, sapply(0:100, SPFO_root_shoot100)/max(sapply(0:100, SPFO_root_shoot100)), rule=2 )

for (i in 1:nSpecies) {
  inputPars$nicheFUN[[i]] = readRDS(as.character(inputPars$nicheFunPath[i]))
  if (inputPars$rsFUN_Inund[[i]] == 'TRUE') {
    if (speciesList[i] == 'SpaFol') {
      inputPars$rsFUN[[i]] = function(h, rs)
        spfoRS(h)*rs
    } else if (speciesList[i] == 'SalPac'){
      inputPars$rsFUN[[i]] = function(h, rs)
        sapaRS(h)*rs
    } else {
      inputPars$rsFUN[[i]] = function(h, rs)
        h * rs
    }
  } else {
    inputPars$rsFUN[[i]] = function(h, rs)
      rs
  }
}

inputPars$aliveIndex = 1:nSpecies
inputPars$deadIndex = (nSpecies + 1):(nSpecies * 2)

