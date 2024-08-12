library(ggplot2)

generate_niche <- function(lower_bound, upper_bound) {
  elevation <- seq(lower_bound, upper_bound, length.out = 100)
  niche <- -0.001 * (elevation - (lower_bound + upper_bound) / 2)^2 + 1
  niche <- (niche - min(niche)) / (max(niche) - min(niche))
  
  data <- data.frame(elevation = elevation, niche = niche)
  
  return(data)
}

# Define species niche boundaries  units: z* [MSL/TR]; multiply by tide range to convert to elevation relative to MSL.
species1_lower_bound <- 0
species1_upper_bound <- 1.03

species2_lower_bound <- 0.73
species2_upper_bound <- 1.05

species3_lower_bound <- 0.85
species3_upper_bound <- 2

# Generate species niches
species1_niche <- generate_niche(species1_lower_bound, species1_upper_bound)
species2_niche <- generate_niche(species2_lower_bound, species2_upper_bound)
species3_niche <- generate_niche(species3_lower_bound, species3_upper_bound)


species1_nicheFun = approxfun(x=species1_niche$elevation*(TR/2)*100, y=species1_niche$niche, rule=2)
species2_nicheFun = approxfun(x=species2_niche$elevation*(TR/2)*100, y=species2_niche$niche, rule=2)
species3_nicheFun = approxfun(x=species3_niche$elevation*(TR/2)*100, y=species3_niche$niche, rule=2)


# Plot the species niches
p1=ggplot() +
  geom_line(data = species1_niche, aes(x = elevation, y = niche, color = "Mangrove")) +
  geom_line(data = species2_niche, aes(x = elevation, y = niche, color = "Marsh")) +
  geom_line(data = species3_niche, aes(x = elevation, y = niche, color = "Upland")) +
  labs(x = "Elevation (z*)", y = "Niche Value") +
  ggtitle("Species Niches Based on Elevation") + theme_classic()+
  scale_color_manual(values = c("Mangrove" = "red", "Marsh" = "green", "Upland" = "blue"))

print(p1)