

#Species flooding functions
{
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
    #inundation.time= inundation.time*100
    if(inundation.time<=50){  # 22 min inundation of SpaFol at Petaluma-- need distributions to connect
      0.0036 * inundation.time^2 -	0.4629*inundation.time	+	14.49
    } else {
      return(0)
    }
  }

  
  #Juncus balticus
  juba.B = function(inundation.time){
    #inundation.time= inundation.time*100
    max(0, -0.02807	*inundation.time^2+	1.10969*inundation.time	+	18.5437 )
  }
  
  #Scheonplectus americanus
  scam.B = function(inundation.time){
    #inundation.time= inundation.time*100
    max(0, -0.0802*inundation.time^2+	3.5626*inundation.time	+	40.203)
  }
  
  #Spartina foliosa
  spfo.B = function(inundation.time){
    #inundation.time= inundation.time*100
    if(inundation.time>1){  #10
      return( max(0, -0.00798*inundation.time^2+	0.60243*inundation.time	+	1.08638))
    }  else {
      return(0)
    }
  }
  
  #Scheonoplectus actus
  scac.B = function(inundation.time){
    #inundation.time= inundation.time*100
    #  -0.0149*inundation.time + 1.4842  #SCAC only, LS
    max(0, -0.028*inundation.time^2 + 1.8289*inundation.time + 28.289 ) #SCAC from LS, and highest elevation SCAM from CJN
  }
  
  #Spartina patens; Kirwan and Guntenspergen 2015, for TR site
  sppa.B = function(inundation.time){
    # inundation.time= inundation.time*100
    if(inundation.time<0.5){  #0.5%
     return(0.02)
    } else {
      return(max(0,(-0.35*inundation.time+26.18)/26.18))
    }
  }
}

#convert to z*, based on site-specific TR and inundation function

inund2zStar = function(inundFun, spInundFun, TR){
  z = seq(-TR,TR, 0.01)*100
  inund = sapply(z, inundFun)
  zStar = z/(TR*100/2)
  spFUN=sapply(inund*100, spInundFun)
  spFUN=spFUN/max(spFUN)
  approxfun(x=zStar, y=spFUN, rule=2)
}


inund2zStar_SAPA = function(inundFun, spInundFun, TR){
  z = seq(-TR,TR/2, 0.01)*100
  inund = sapply(z, inundFun)
  zStar = z/(TR*100/2)
  spFUN=sapply(inund*100, spInundFun)
  spFUN=spFUN/max(spFUN)
  approxfun(x=zStar, y=spFUN, rule=2)
}

spfo.zStar = inund2zStar(inundFun, spfo.B, TR)
#plot(seq(-2,1,0.01), sapply(seq(-2,1,0.01), spfo.zStar), xlab='Elevation (zStar)', ylab='Niche', main='SpaFol', type='l')

sapa.zStar = inund2zStar_SAPA(inundFun, sapa.B, TR)
# plot(seq(-1,2,0.01), sapply(seq(-1,2,0.01), sapa.zStar), type='l', xlab='z*', ylab='Aboveground Biomass')
# lines(seq(-1,2,0.01), sapply(seq(-1,2,0.01), spfo.zStar)*0.7, col='red')
# 
# plot(seq(-1,2,0.01), sapply(seq(-1,2,0.01), sapa.zStar)*1.5/0.37, type='l', xlab='z*', ylab='Belowground Biomass')
# lines(seq(-1,2,0.01), sapply(seq(-1,2,0.01), spfo.zStar)*0.7/1.33, col='red')



#nicheFUN - the elevation range where the species can live. Controls % cover. should range from 0-1 and be relative to zStar [ zStar = (z-MSL)/(MHHW-MSL)  or (z-MSL)/(tideRange/2)]
#biomassFUN - the relationship between maximum biomass and elevation (zStar), while accounting for the % of carrying capacity. Returns standing stock of aboveground biomass in g/cm2


#Species input parameters
{
  
  #Salicornia pacifica
  SalPac = list(Type='Grass',
                nicheFUN=readRDS('input\\sapaNiche_zStar.rds'),  #functional relationship between z* and species presence [0-1] 
                salT=40,
                biomassFUN  =function(zStar, cc) {sapa.zStar(zStar)*2500/10000*cc },  #aboveground biomass. functional scalor should be in g/cm2    #3500  #Mowry=1700  #Guad=4000, tri=4500
                absGrowth = 0.15,  #maximum growth rate increase
                strucRoots= 0.05, #0     #very large proportion of root biomass  
                        #remaining root class proportions must add to 1
                coarseRoots= 0.3,  
                smallRoots= 0,
                fineRoots = 0.7,
                deadFallRate = 0.4,  #rate that standing dead fall   
                r.roots = 40,         #defines rooting depth (cm)
                rooting_shape = 'lin',  #'exp' or 'lin'  for exponential or linear. linear best for herbaceous veg, exponential best for mangroves
                rooting_shape_coef = NA, #decay coefficient for distributing roots. For 'exp' shape only. Try -3, -5, -10
                #  root_shoot = 0.3, #function(h) 0.37, #0.174,    #Shoot:Root ratio  #Janousek et al 2016   0.5149- i100*0.03892 + 0.0005896*i100^2
                root_shoot = function(i100) exp(min(-0.22, -1.5649- i100*0.0005929 + 0.0001965*i100^2)), # Janousek et al 2016, function takes in inundation duration % [0-100%]
                root_por = 0.1,       #root porosity -- defines amount of root collapse upon plant death. Higher value=more collapse [0-1] GET BETTER ESTIMATE!
                kf = 0.95,     #0.9    #fine root turnover
                km = 0,    #0.4       #small root turnover
                kc = 0.5,    #0.2     #coarse root turnover  
                ks = 0.05,    #0.02   #very large root turnover (structural, cable)
                kfp = 0.8,   #0.5     #dead fine root to particulate 
                fc1 = 0.4,  #0.4,     #labile fraction of fine and small roots
                fc2 = 0.6,  #0.5,     #labile fraction of coarse & structural roots  
                leafRefrac = 0.15,    #leaf refractory fraction 
                leafProp=0.75,         #proportion of litter that is leaves [vs branches]
                litterDep =0.1,       #litter deposition rate
                litterRetain = 1,     #maximum proportion of litter that remains on-site 
                seeding = 0.01,      #initial seeding rate
                inundLimReseed=0.5,  #seeding inundation limit. proportion inundation time [0-1]
                compA = 1,           #scalor for space competition. Should average to 1 across species considered 
                wood.den=0.4,        #wood density, g/cm3. for structural roots/rhizomes
                root.den = 0.2,      #root density, g/cm, for fine, small and coarse roots
                initCover=0.01,      #proportion of carrying capacity
                RSslp=1,             #root:shoot slope - how sensitive is the root:shoot to distance from optimum conditions
                mann= 0.06 # not used
  )
  
  #particulate rates get slower with increased root size. Small roots get put into the Fine root particulate pool
  SalPac$kcp=SalPac$kfp*0.66
  SalPac$ksp = SalPac$kfp*.33
  
  
  #Spartina foliosa
  SpaFol = list(Type='Grass',
                nicheFUN=readRDS('input\\spfoNiche_zStar.rds'), #where can the species grow;  
                salT = 40,
                biomassFUN  =function(zStar, cc) {spfo.zStar(zStar)*1800/10000*cc}, #based on inundation tolerances, marsh organ information
                absGrowth = 0.2,
                strucRoots= 0.05, # #very large proportion of root biomass    
                coarseRoots= 0.3,
                smallRoots= 0,
                fineRoots = 0.7,
                deadFallRate = 0.95,  #rate that standing dead falls   
                r.roots = 30,         #defines rooting depth (cm)
                rooting_shape = 'lin',
                rooting_shape_coef = NA,
                #  root_shoot = 1.33,#function(h) 1.33,    #Shoot:Root ratio  #Janousek et al 2016
                root_shoot =  function(i100) exp(min(1,(0.5149- i100*0.03892 + 0.0005896*i100^2))),  # Janousek et al 2016, function takes in inundation duration % [0-100%]
                root_por = 0.1,    #root porosity -- defines amount of root collapse upon tree death. GET BETTER
                kf = 0.8,     #0.9        #fine root turnover
                km = 0,    #0.4         #small root turnover
                kc = 0.5,    #0.2         #coarse root turnover  
                ks = 0.05,    #0.02        #very large root turnover (structural, cable, rhizome)
                kfp= 0.8,                  #conversion rate of dead fine roots to particulate organic matter 
                fc1 =0.4,  #0.4           #labile fraction of fine root
                fc2= 0.6, #0.5            #labile fraction of coarse root  
                leafRefrac = 0.15,    #leaf refractory fraction 
                leafProp=1,
                litterDep =0.1,     #litter deposition rate
                litterRetain = 1,     #maximum proportion of litter that remains on-site. Scalor for considering herbivory 
                seeding = 0.01,          #initial seeding rate
                inundLimReseed = 0.6,  #seeding inundation limit. proportion inundation time
                compA = 1,          #scalor for space competition 
                wood.den= 0.4,       #wood/rhizome density, g/cm3
                root.den = 0.2,       #root density, g/cm3
                initCover=0.01,       #proportion of carrying capacity
                RSslp=1,
                mann= 0.04 #MADE UP, FIND BETTER VALUE
  )
  
  SpaFol$kcp=SpaFol$kfp*0.66
  SpaFol$ksp=SpaFol$kfp*.33
  
  
  #Bolbaceneous maritima
  BolMar = list(Type='Marsh',
                nicheFUN=boma.B, 
                salT=25,
                biomassFUN = function(zStar, cc) {cc*(1976/10000)}, #aboveground production g/cm2 -from Buffington clip plots at Petaluma. 75th percentile from clip plots. multiply by 0.5 to account for perennial nature
                absGrowth = 0.3/2, #maximum change in cover per year
                strucRoots= 0.1,  #0     #very large proportion of root biomass    
                coarseRoots= 0.4,
                smallRoots= 0,
                fineRoots = 0.6,   #fine-coarse root proportion of remaining root bioass (1-vBigRoots)
                deadFallRate = 0.95,  #rate that standing dead fall   
                r.roots = 40,         #defines rooting depth (cm)
                rooting_shape = 'lin',
                rooting_shape_coef = NA,
                root_shoot =function(h) 4.54,    #root:shoot ratio  #Janousek et al 2016
                root.por = 0.01,    #root porosity -- defines amount of root collapse upon tree death
                kf = 0.5,             #fine root turnover
                km=0,                 #small root tunnover
                kc = 0.2,             #coarse root turnover  
                ks = 0,               #very large root turnover (structural, cable)
                kfp=0.8,              #dead fine root rate to particulate OM
                kfp= 0.4,             #dead small root rate to particulate OM
                ktt=0.0,              #dead structural root rate to particulate OM
                fc1 =0.75,            #labile fraction of fine root
                fc2= 0.85,            #labile fraction of coarse root  
                leafRefrac = 0.15,    #leaf refractory fraction 
                leafProp=1,
                litterDep =0.01 ,     #litter deposition rate
                litterRetain = 0.75,     #maximum proportion of litter that remains on-site 
                seeding = 0,          #initial seeding rate
                compA = 0.2,          #scalor for space competition 
                wood.den=0.2,         #wood density, g/cm3
                root.den = 0.2,       #root density, g/cm3
                initCover=0.01,
                RSslp=1,
                mann = 0.03 #not used
  )
  
  BolMar$kcp = BolMar$kfp*0.66
  BolMar$ksp = BolMar$kfp*0.33
  

  #Juncas balticus
  JunBal = list(Type='Marsh',
                nicheFUN=juba.B, 
                salT=35,
                biomassFUN = function(zStar, cc) {cc*(200/10000)}, #g/cm2 - MADE UP
                absGrowth = 0.2,
                strucRoots= 0.0,     
                coarseRoots= 0.4,
                smallRoots= 0,
                fineRoots = 0.6, 
                deadFallRate = 0.95,  #rate that standing dead fall   
                r.roots = 30,         #defines rooting depth (cm)
                rooting_shape = 'lin',
                rooting_shape_coef = NA,
                root_shoot = function(h) 2.13,    #Shoot:Root ratio  #Janousek et al 2016
                root.por = 0.01,    #root porosity -- defines amount of root collapse upon tree death. GET BETTER
                kf = 0.5,             #fine root turnover
                km = 0,               #small root turnover
                kc = 0.2,             #coarse root turnover  
                ks = 0,               #very large root turnover (structural, cable)
                kfp=0.8,
                fc1 =0.75,            #labile fraction of fine root
                fc2= 0.85,            #labile fraction of coarse root  
                leafRefrac = 0.15,    #leaf refractory fraction
                leafProp=1,
                litterDep =0.01,     #litter deposition rate
                litterRetain = 1,     #maximum proportion of litter that remains on-site 
                seeding = 0,          #initial seeding rate
                inundLimReseed = 0.5, #proportion inundation time
                compA = 0.2,          #scalor for space competition 
                wood.den=0.2,         #wood density, g/cm3
                root.den = 0.2,       #root density, g/cm3
                initCover=0.01,
                RSslp=1,
                mann=0.03
  )
  
  JunBal$kcp = JunBal$kfp*0.66
  JunBal$ksp = JunBal$kfp*0.33
  
  
  SchAme =  list(Type='Marsh',
                 nicheFUN=scam.B, 
                 biomassFUN = function(zStar, cc) {cc*(500/10000)}, #g/cm2 - MADE UP
                 absGrowth = 0.5,
                 strucRoots= 0,  #0     #very large proportion of root biomass    
                 coarseRoots= 0.4,
                 smallRoots= 0,
                 fineRoots = 0.6, 
                 deadFallRate = 0.95,  #rate that standing dead fall -- graminiods falls fast   
                 r.roots = 40,         #defines rooting depth (cm)
                 rooting_shape = 'lin',
                 rooting_shape_coef = NA,
                 root_shoot = function(h) 1.67,    #Shoot:Root ratio  #Janousek et al 2020
                 root.por = 0.01,    #root porosity -- defines amount of root collapse upon tree death. GET BETTER
                 kf = 0.5,             #fine root turnover
                 km = 0,
                 kc = 0.2,             #coarse root turnover  
                 ks = 0,               #very large root turnover (structural, cable)
                 kfp=0.8,
                 fc1 =0.75,            #labile fraction of fine root
                 fc2= 0.85,            #labile fraction of coarse root  
                 leafRefrac = 0.15,    #leaf refractory fraction 
                 leafProp=1,
                 litterDep =0.01 ,     #litter deposition rate
                 litterRetain = 1,     #maximum proportion of litter that remains on-site 
                 seeding = 0,          #initial seeding rate
                 inundLimReseed = 0.75, #proportion inundation time
                 compA = 0.2,          #scalor for space competition 
                 wood.den=0.2,         #wood density, g/cm3
                 root.den = 0.2,       #root density, g/cm3
                 initCover=0.01,
                 RSslp=1,
                 mann=0.03
  )
  
SchAme$kcp = SchAme$kfp*0.66
SchAme$ksp = SchAme$kfp*0.33
  
  SchAcu =  list(Type='Marsh',
                 nicheFUN=scac.B, 
                 biomassFUN = function(zStar,cc) {cc*(500/10000)}, #g/cm2 - MADE UP
                 absGrowth = 0.4,
                 
                 vBigRoots= 0,       #very large proportion of root biomass     
                 fineBigRoots = 0.6,   #fine-coarse root proportion of remaining root bioass (1-vBigRoots)
                 deadFallRate = 0.95,  #rate that standing dead fall   
                 r.roots = 40,         #defines rooting depth (cm)
                 rooting_shape = 'lin',
                 rooting_shape_coef = NA,
                 root_shoot =function(h) 1.67,    #Shoot:Root ratio  #Janousek et al 2020
                 root.por = 0.01,    #root porosity -- defines amount of root collapse upon tree death. GET BETTER
                 kma=0.05,
                 kr = 0.5,             #fine root turnover
                 km = 0.2,             #coarse root turnover  
                 kt = 0,               #very large root turnover (structural, cable)
                 krr=0.8,
                 kfp= 0.4,
                 ktt=0.0,
                 fc1 =0.75,            #labile fraction of fine root
                 fc2= 0.85,            #labile fraction of coarse root  
                 leafRefrac = 0.15,    #leaf refractory fraction 
                 leafProp=1,
                 litterDep =0.01 ,     #litter deposition rate
                 litterRetain = 1,     #maximum proportion of litter that remains on-site 
                 seeding = 0,          #initial seeding rate
                 compA = 0.2,          #scalor for space competition 
                 wood.den=0.2,         #wood density, g/cm3
                 root.den = 0.2,       #root density, g/cm3
                 initCover=0.01,
                 RSslp=1,
                 mann=0.03
  )
  
  SpaPat =  list(Type='Marsh',
                 nicheFUN=sppa.B,  #CHANGE
                 salT = 40,
                 biomassFUN = function(zStar,cc) {cc/0.0048*(0.025)}, #g/cm2 -  #Kirwin and Guntenspergen 2015-- marsh organ results: 0.162 g/cm2; natural marsh max: ~750 g/m2 (fig 8)
                 absGrowth = 0.8,
                 vBigRoots= 0,       #very large proportion of root biomass     
                 fineBigRoots = 0.6,   #fine-coarse root proportion of remaining root bioass (1-vBigRoots)
                 deadFallRate = 0.95,  #rate that standing dead fall   
                 r.roots = 30,         #defines rooting depth (cm)
                 root_shoot = function(h) 3,    #Shoot:Root ratio  #Kirwin and Guntenspergen 2015, est from fig 7
                 root.por = 0.01,    #root porosity -- defines amount of root collapse upon tree death. GET BETTER
                 kma=0.05,
                 kr = 0.8,             #fine root turnover
                 km = 0.5,             #coarse root turnover  
                 kt = 0.001,               #very large root turnover (structural, cable)
                 krr=0.8,
                 kfp= 0.4,
                 ktt=0.0,
                 fc1 =0.75,            #labile fraction of fine root
                 fc2= 0.85,            #labile fraction of coarse root  
                 leafRefrac = 0.15,    #leaf refractory fraction 
                 leafProp=1,
                 litterDep =0.01 ,     #litter deposition rate
                 litterRetain = 1,     #maximum proportion of litter that remains on-site 
                 seeding = 0,          #initial seeding rate
                 compA = 1,          #scalor for space competition 
                 wood.den=0.2,         #wood density, g/cm3
                 root.den = 0.2,       #root density, g/cm3
                 initCover=0.01,
                 RSslp=1,
                 mann = 0.04 #made up; get better number
  )
  
  SpaAlt =  list(Type='Marsh',
                 nicheFUN=sppa.B, 
                 salT = 40,
                 biomassFUN = function(zStar,cc) {cc*(0.2)}, #g/cm2, MADE UP, get reference
                 absGrowth = 0.02,
                 vBigRoots= 0,       #very large proportion of root biomass     
                 fineBigRoots = 0.6,   #fine-coarse root proportion of remaining root bioass (1-vBigRoots)
                 deadFallRate = 0.95,  #rate that standing dead fall   
                 r.roots = 40,         #defines rooting depth (cm)
                 root_shoot = function(h) 4,    #Shoot:Root ratio  #Kirwin and Guntenspergen 2015
                 root.por = 0.01,    #root porosity -- defines amount of root collapse upon tree death. GET BETTER
                 kr = 0.5,             #fine root turnover
                 km = 0.2,             #coarse root turnover  
                 kt = 0,               #very large root turnover (structural, cable)
                 krr=0.8,
                 kfp= 0.4,
                 ktt=0.0,
                 fc1 =0.75,            #labile fraction of fine root
                 fc2= 0.85,            #labile fraction of coarse root  
                 leafRefrac = 0.15,    #leaf refractory fraction 
                 leafProp=1,
                 litterDep =0.01 ,     #litter deposition rate
                 litterRetain = 1,     #maximum proportion of litter that remains on-site 
                 seeding = 0,          #initial seeding rate
                 compA = 0.2,          #scalor for space competition 
                 wood.den=0.2,         #wood density, g/cm3
                 root.den = 0.2,       #root density, g/cm3
                 initCover=0.01,
                 RSslp=1,
                 mann = 0.04  #MADE UP, FIND BETTER VALUE
  )

}