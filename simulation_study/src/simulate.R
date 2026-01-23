library(purrr)

rm(list = ls())

source("functions/utilities.R")
source("functions/simulate_PBDOM.R")
load("simulation_study/data/distanceMatrix.Rdata")

N <- dim(d)[1]
T <- 6
K <- 4

params <- data.frame(omega = 0.33, rho = 0.66, alpha = 7.5,
                     sigma = c(12.5,25,12.5,25),
                     alpha_gam = logit(0.5),
                     beta_gam = c(0,0,5,5))

X <- terra::rast(list.files("simulation_study/data/X", full.names = TRUE)[runif(1,1,10)])

simul <-  pmap(params, simulate_PBDOM, 
               N = N, T = T, d = d, K = K, X = c(scale(terra::values(X))), A = 100)

names(simul) <- paste0("Scenario", 1:nrow(params))

### SAVE -----------------------------------------------------------------------

id_simul <- floor(runif(1,0,100000))
dir.create(paste0("simulation_study/out/sim_", id_simul))
walk(names(simul), ~ dir.create(paste0("simulation_study/out/sim_", id_simul, "/", .x)))
iwalk(simul, ~ saveRDS(.x, file=paste0("simulation_study/out/sim_", id_simul,"/", .y, "/data.rds")))

### PLOTS ----------------------------------------------------------------------

png(paste0("simulation_study/out/sim_", id_simul, "/X.png"),
    width=21, height=21, units = "cm", res = 300)
terra::plot(X)
dev.off()

png(paste0("simulation_study/out/sim_", id_simul, "/Z.png"),
    width=29.7, height=21, units = "cm", res = 300)
par(mfrow = (c(nrow(params),T)))
pwalk(expand.grid(t  = 1:T, s = 1:nrow(params)), function(t, s){
  r <- terra::rast(xmin = 0, xmax = 300, ymin = 0, ymax = 300, resolution = c(10, 10))
  terra::values(r) <- simul[[s]]$z[t, ]
  terra::plot(r, legend = FALSE)})
dev.off()

png(paste0("simulation_study/out/sim_", id_simul, "/gamma.png"),
    width=29.7, height=21, units = "cm", res = 300)
par(mfrow = (c(nrow(params),T-1)))
pwalk(expand.grid(t  = 1:(T-1), s = 1:nrow(params)), function(t, s){
  r <- terra::rast(xmin = 0, xmax = 300, ymin = 0, ymax = 300, resolution = c(10, 10))
  terra::values(r) <- simul[[s]]$gamma[t, ]
  terra::plot(r, legend = FALSE)})
dev.off()