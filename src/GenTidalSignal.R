
#generates inundation duration as a function of elevation based on tidal range and type of tide [diurnal, semi-diurnal, mixed semi-diurnal]

library(ggplot2)

generate_tidal_signal <- function(tidal_range, time_range, time_interval, constituents) {
  # Calculate the number of data points
  num_points <- ceiling((time_range[2] - time_range[1]) / time_interval) + 1
  
  # Generate time points
  time <- seq(time_range[1], time_range[2], by = time_interval)
  
  # Generate tidal signal based on harmonic constituents
  water <- rep(0, num_points)
  for (i in 1:nrow(constituents)) {
    constituent <- constituents[i, ]
    amplitude <- tidal_range/2 * constituents$amplitude[i] 
    phase <- constituents$phase[i]
    speed <- constituents$speed[i]
   # tidal_signal <- tidal_signal + amplitude * sin(2 * pi * (time - phase) * speed)
    water = water + (amplitude * cos((speed * pi / 180) * time + phase))
  }
  
  # Create a data frame with time and tidal signal
  data <- list(time = time, water = water)
  
  return(data)
}

# Define parameters
#tidal_range <- 4 # Input tidal range (in meters)
time_range <- c(0, 24*365) # Time range (in hours)
time_interval <-0.1#  0.25 # Time interval (in hours)
#tide_type <- "mixed semi-diurnal" # Choose tide type: "diurnal", "semi-diurnal", or "mixed semi-diurnal"
##tide_type <- "diurnal" # Choose tide type: "diurnal", "semi-diurnal", or "mixed semi-diurnal"
#tide_type <- "semi-diurnal" # Choose tide type: "diurnal", "semi-diurnal", or "mixed semi-diurnal"

# Define constituent properties
constituent_properties <- data.frame(constituent = c("M2", "S2", "N2", "K2"),
                                     amplitude = c(0.7, 0.3, 0.2, 0.1),
                                     phase = c(110, 121, 102, 118),
                                     speed = c(28.984106, 30, 28.43973, 30.082138))

# Subset constituents based on tide type
constituents <- switch(tide_type,
                       "diurnal" = constituent_properties[1, ],
                       "semi-diurnal" = constituent_properties[1:2, ],
                       "mixed semi-diurnal" = constituent_properties[1:4, ],
                       #'user-defined' = data.frame(constituent=HC$Constituent, amplitude=HC$Amplitude, phase=HC$Phase, speed=HC$Speed),
                       'user-defined' = data.frame(amplitude=HC$Amplitude, phase=HC$Phase, speed=HC$Speed),
                       stop("Invalid tide type. Please choose 'diurnal', 'semi-diurnal', or 'mixed semi-diurnal'."))

# Generate tidal signal
if(tide_type=='user-defined'){
theLevels <- generate_tidal_signal(2, time_range, time_interval, constituents)
} else {
  theLevels <- generate_tidal_signal(TR, time_range, time_interval, constituents)
}

theLevels$sec=time_interval*60*60

# Plot the tidal signal
ggplot()  +
  geom_line( aes(x = theLevels$time, y = theLevels$water)) +
  labs(x = "Time", y = "Tidal Signal") +
  ggtitle(paste("Tidal Signal with Tidal Range =", TR, "meters and Tide Type =", tide_type))


zRange=seq(min(theLevels$water)-0.1, max(theLevels$water)+0.1, 0.01)
inundation.duration<- sapply(zRange, function(z) sum(theLevels$water > z, na.rm=T))
inundDur = inundation.duration/max(inundation.duration)
inundFun = approxfun((inundDur)~(zRange*100), rule=2)  #flooding freq

plot(zRange,inundDur, type='l', ylab='Inundation Time',xlab='Elevation rMSL m', main=tide_type)




#IF =  sapply(zRange, function(z) sum(rleid( ifelse(diff(theLevels$water)>0 & theLevels$water>z,1,0))  ))


