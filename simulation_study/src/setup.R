library(tidyverse)

rm(list = ls())

grid <- sf::st_make_grid(x = sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 1, ymax = 1)),
                         cellsize = 10, n = c(30,30))

r <- terra::rast(xmin = 0, xmax = 300, ymin = 0, ymax = 300, resolution = c(10, 10))
d <- as.matrix(dist(terra::crds(r)))

walk(1:10, r = r, d = d, sigma = 15, function(.x, r, d, sigma){
  tmp <- r
  terra::values(tmp) <- MASS::mvrnorm(n = 1, mu = rep(0,length(terra::values(r))),
                                      Sigma =  1 / sigma**2 * exp(-d/(2 * sigma**2)))
  terra::writeRaster(tmp, paste0("simulation_study/dat/X/X", .x,".tiff"), overwrite = TRUE)
})

rm(r)

save.image("simulation_study/dat/grid.Rdata")