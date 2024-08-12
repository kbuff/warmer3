# Define the sigmoid function
sigmoid <- function(x, a, b, c, d) {
  y <- a + (b - a) / (1 + exp(-c * (x - d)))
  return(y)
}

# Define the parameters for the first sigmoid function (0-1)
a1 <- 0   # Lower asymptote
b1 <- 1   # Upper asymptote
c1 <- 10   # Steepness
d1 <- 0.5   # Midpoint

# Define the parameters for the second sigmoid function (1-0)
a2 <- 1   # Lower asymptote
b2 <- 0   # Upper asymptote
c2 <- 20   # Steepness
d2 <- 0.5   # Midpoint

# Define the range of y values
x <- seq(0, 1, length.out = 100)

# Calculate the y values for each sigmoid function
y1 <- sigmoid(x, a1, b1, c1, d1)
y2 <- sigmoid(x, a2, b2, c2, d2)

# Plot the combined curve
plot(x, y1, type = "l", col = "blue", xlim = c(0, 1), ylim = c(0, 1), xlab = "x", ylab = "y", main = "Combined Sigmoid Functions")
lines(x, y2, col = "red")
legend("topright", legend = c("0-1", "1-0"), col = c("blue", "red"), lty = 1)



##
# spfo_sigmoid <- function(x) {
#   if (x <= 0.65) {
#     return(1 / (1 + exp(-(20 * (x - 0.4)))) )
#   } else {
#     return(1 / (1 + exp(15 * (x - 0.95))))
#   }
# }

spfo_sigmoid <- function(x) {
  if (x <= 0.57) {
    return(1 / (1 + exp(-(30 * (x - 0.35)))) )
  } else {
    return(1 / (1 + exp(15 * (x - 0.85))))
  }
}

# Generate x values -- zStar
x <- seq(-1, 1.5, length.out = 200)

# Calculate y values using the combined sigmoid function
ySPFO <- sapply(x, spfo_sigmoid)

# Plot the combined sigmoid curve
plot(x, ySPFO, type = "l", col = "blue", xlim = c(-1, 1), ylim = c(0, 1), xlab = "x", ylab = "y", main = "Combined Sigmoid Functions")

  
##
# sapa_sigmoid <- function(x) {
#   if (x <= 0.95) {
#     return(1 / (1 + exp(-(25 * (x - 0.65)))) )
#   } else {
#     return(1 / (1 + exp(20 * (x - 1.25))))
#   }
# }

sapa_sigmoid <- function(x) {
  if (x <= 0.95) {
    return(1 / (1 + exp(-(30 * (x - 0.6)))) )
  } else {
    return(1 / (1 + exp(25 * (x - 1.35))))
  }
}

# Generate x values -- zStar
x <- seq(-1, 1.5, length.out = 200)
x=hSPFO$mids

# Calculate y values using the combined sigmoid function
ySAPA <- sapply(x, sapa_sigmoid)

# Plot the combined sigmoid curve
plot(x*(TR/2), ySAPA, type = "l", col = "blue", xlim = c(-1, 3.5), ylim = c(0, 1), xlab = "zStar", ylab = "Niche", main = " Cover Functions")
lines(x*(TR/2), ySPFO, type = "l", col = "red", xlim = c(-1, 3.5), ylim = c(0, 1), xlab = "x", ylab = "y")

ySAPA= ifelse(ySAPA<1e-3,0,ySAPA)
ySPFO= ifelse(ySPFO<1e-3,0,ySPFO)

sapaFUN = approxfun(x=x, y=ySAPA, rule=2)
spfoFUN = approxfun(x=x, y=ySPFO, rule=2)

saveRDS(sapaFUN, 'C:\\MERCC\\MERCC_SFBay\\input\\sapaNiche_zStar.rds')
saveRDS(spfoFUN, 'C:\\MERCC\\MERCC_SFBay\\input\\spfoNiche_zStar.rds')


DonEdAll <- read.csv("C:\\Kevin\\1_Projects\\1_SFBay_Region\\2_DonEdwards\\Data\\Survey_Data_Final.csv")

DESub <- subset(DonEdAll, Marsh.Site == 'Mowry Marsh Slough' |  Marsh.Site== 'Mowry Marsh South' | Marsh.Site=='Mowry Marsh North' |
                        Marsh.Site=='Guadalupe Mouth' |  Marsh.Site=='Guadalupe to Stevens'  | Marsh.Site=='Laumeister' |  Marsh.Site=='Faber' | Marsh.Site=='Calaveras')

spfoALL = subset(DonEdAll, DonEdAll$SPFO..>30)
sapaALL  = subset(DonEdAll, DonEdAll$SAPA..>30)
par(bg='white')

hist((spfoALL$Z_MHW_m+1.03)/(2.5/2), breaks=seq(0,1.5,0.05),freq=T, main='Spartina foliosa Presence', xlab='z*', col='grey50')
hist((sapaALL$Z_MHW_m+1.03)/(2.5/2), breaks=seq(0,1.5,0.05), freq=T,main='Sarcocornia pacifica Presence', xlab='z*', col='grey50')


hSPFO = hist((spfoALL$Z_MHW_m+1.03)/(2.5/2), breaks=seq(-0.05,2.5,0.05) )
hSAPA = hist((sapaALL$Z_MHW_m+1.03)/(2.5/2), breaks=seq(-0.05,2.5,0.05))
hALL =  hist((DonEdAll$Z_MHW_m+1.03)/(2.5/2), breaks=seq(-0.05,2.5,0.05))
hALL$counts/max(hALL$counts)

#hSPFO$breaks
relSPFO=hSPFO$counts/hALL$counts
relSAPA = hSAPA$counts/hALL$counts
relSAPA = ifelse(is.nan(relSAPA), 0, relSAPA)
relSPFO=ifelse(is.nan(relSPFO), 0, relSPFO)
relSAPA[9]=0
relSAPA = relSAPA/max(relSAPA)
relSPFO= relSPFO/max(relSPFO)



dat = data.frame(x=hSPFO$mids, SPFO = relSPFO, SAPA=relSAPA)

pCOV=ggplot(data=dat, aes(x=x)) + geom_histogram( aes(y=SPFO), stat='identity' , alpha=0.2, fill='darkblue') +  theme_bw()+ ylab('Relative Cover')+ xlab('Elevation (z*)')+ xlim(c(0,1.6))+
  geom_histogram( aes(y=SAPA), stat='identity' , alpha=0.2, fill='darkgreen') + geom_line( aes(x=x, y=ySAPA), col='darkgreen', size=1)+ geom_line( aes(x=x, y=ySPFO), size=1,col='darkblue')+
   theme(rect = element_rect(fill = "transparent",   colour = NA_character_), text=element_text(size=20)  ) 
pCOV

png('C:\\MERCC\\MERCC_SFBay\\output\\RelativeCover_Hist_lgerText.png', units='mm', bg='transparent', res=800, height=100, width=140)
pCOV
dev.off()
