
spfo_sigmoid <- function(x) {
  if (x <= 0.57) {
    return(1 / (1 + exp(-(30 * (x - 0.35)))) )
  } else {
    return(1 / (1 + exp(15 * (x - 0.85))))
  }
}

# Generate x values -- zStar
x <- seq(-1, 2, length.out = 200)

# Calculate y values using the combined sigmoid function
ySPFO <- sapply(x, spfo_sigmoid)

# Plot the combined sigmoid curve
plot(x, ySPFO, type = "l", col = "blue", xlim = c(-1, 1), ylim = c(0, 1), xlab = "x", ylab = "y", main = "Combined Sigmoid Functions")


sapa_sigmoid <- function(x) {
  if (x <= 0.95) {
    return(1 / (1 + exp(-(30 * (x - 0.6)))) )
  } else {
    return(1 / (1 + exp(25 * (x - 1.35))))
  }
}

# Generate x values -- zStar
x <- seq(-1, 2, length.out = 200)
#x=hSPFO$mids

# Calculate y values using the combined sigmoid function
ySAPA <- sapply(x, sapa_sigmoid)

# Plot the combined sigmoid curve
plot(x, ySAPA, type = "l", col = "blue", xlim = c(-1, 2), ylim = c(0, 1), xlab = "zStar", ylab = "Niche", main = " Cover Functions")
lines(x, ySPFO, type = "l", col = "red", xlim = c(-1, 2), ylim = c(0, 1), xlab = "x", ylab = "y")

ySAPA= ifelse(ySAPA<1e-3,0,ySAPA)
ySPFO= ifelse(ySPFO<1e-3,0,ySPFO)

sapaFUN = approxfun(x=x, y=ySAPA, rule=2)
spfoFUN = approxfun(x=x, y=ySPFO, rule=2)

saveRDS(sapaFUN, 'input/speciesFuns/Niche/sapaNiche_zStar.rds')
saveRDS(spfoFUN, 'input/speciesFuns/Niche/spfoNiche_zStar.rds')

#https://www.sciencebase.gov/catalog/item/660b32efd34e4df16bd58cd4
plantSurvey= read.csv('input/speciesFuns/baydelta_plants_cover_height.csv')
plantSurveySub=subset(plantSurvey, Region=='San Pablo Bay'  | Region=='Central SF Bay')

spfoALL = subset(plantSurveySub, plantSurveySub$Spartina.foliosa.Cover>=50)
sapaALL  = subset(plantSurveySub, plantSurveySub$Salicornia.pacifica.Cover>=50)
par(bg='white')

hist(spfoALL$Zstar, breaks=seq(-.6,2.3,0.05),freq=T, main='Spartina foliosa Presence', xlab='z*', col='grey50')
hist(sapaALL$Zstar, breaks=seq(-.6,2.3,0.05), freq=T,main='Sarcocornia pacifica Presence', xlab='z*', col='grey50')

hSPFO = hist(spfoALL$Zstar, breaks=seq(-0.6,3.5,0.05) )
hSAPA = hist(sapaALL$Zstar, breaks=seq(-0.6,3.5,0.05))
hALL =  hist(plantSurveySub$Zstar, breaks=seq(-0.6,3.5,0.05))

relSPFO = hSPFO$counts/max(hSPFO$counts)
relSAPA = hSAPA$counts/max(hSAPA$counts)
relSAPA = ifelse(is.nan(relSAPA), 0, relSAPA)
relSPFO = ifelse(is.nan(relSPFO), 0, relSPFO)

dat = data.frame(x=hSPFO$mids, SPFO = relSPFO, SAPA=relSAPA)
datLines = data.frame(x=x, ySAPA=ySAPA, ySPFO=ySPFO)

pCOV=ggplot(data=dat, aes(x=x)) + geom_histogram( aes(y=SPFO), stat='identity' , alpha=0.2, fill='darkblue') + 
  theme_bw()+ ylab('Relative Cover')+ xlab('Elevation (z*)')+ xlim(c(-.5,2))+
  geom_histogram( aes(y=SAPA), stat='identity' , alpha=0.2, fill='darkgreen') +
  geom_line(data=datLines, aes(x=x, y=ySAPA), col='darkgreen', size=1)+ geom_line(data=datLines, aes(x=x, y=ySPFO), size=1,col='darkblue')+
  theme(rect = element_rect(fill = "transparent",   colour = NA_character_), text=element_text(size=20)  ) 
pCOV

png('input/speciesFuns/NicheFun_Hist_SFBay_SPFO-SAPA.png', units='mm', bg='transparent', res=800, height=100, width=140)
pCOV
dev.off()
