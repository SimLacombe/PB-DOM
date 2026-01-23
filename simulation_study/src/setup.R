library(tidyverse)

rm(list = ls())

r <- terra::rast(xmin = 0, xmax = 300, ymin = 0, ymax = 300, resolution = c(10, 10))
d <- as.matrix(dist(terra::crds(r)))

walk(1:10, r = r, d = d, sigma = 10, function(.x, r, d, sigma){
  tmp <- r
  terra::values(tmp) <- c(scale(MASS::mvrnorm(n = 1, mu = rep(0,length(terra::values(r))),
                                              Sigma =  1 / sigma**2 * exp(-d/(2 * sigma**2)))))
  terra::writeRaster(tmp, paste0("simulation_study/data/X/X", .x,".tiff"), overwrite = TRUE)
})

save(d, file = "simulation_study/data/distanceMatrix.Rdata")