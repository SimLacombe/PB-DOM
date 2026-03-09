library(purrr)
library(nimble)

rm(list = ls())

source("functions/utilities.R")
source("functions/simulate_PBDOM.R")
source("functions/nimbleModel_simulations_indepTrials.R")

load("simulation_study/data/distanceMatrix.Rdata")
dsparse <- getSparse(d, thr = 60)

Nsim <- 5

N <- dim(d)[1]
T <- 6
K <- 4

params <- data.frame(omega = 0.33, 
                     rho = 0.66, 
                     alpha = 4,
                     sigma = 20,
                     alpha_gam = logit(0.5),
                     beta_gam = 0)

params <- map_dfr(1:Nsim, ~ params)

X <- terra::rast(list.files("simulation_study/data/X", full.names = TRUE)[4])

simul <-  pmap(params, simulate_PBDOM, 
               N = N, T = T, d = d, K = K, X = c(scale(terra::values(X))), A = 100)

par(mfrow = (c(nrow(params),T)))
pwalk(expand.grid(t  = 1:T, s = 1:nrow(params)), function(t, s){
  r <- terra::rast(xmin = 0, xmax = 300, ymin = 0, ymax = 300, resolution = c(10, 10))
  terra::values(r) <- simul[[s]]$z[t, ]
  terra::plot(r, legend = FALSE)})

for(id in 1:Nsim){
  
  mod <- nimbleModel(
    code = myModel,
    data = list(y = simul[[id]]$y),
    constants = list(nSites = N, 
                     nSeasons = T,
                     nSurveys = 4, 
                     dmatP = dsparse$p,
                     dmatI = unname(dsparse$i),
                     d2 = dsparse$d**2,
                     X = X),
    inits = initial.values(zst = array(1, dim = c(myConstants$nSeasons, myConstants$nSites))),
    calculate = FALSE
  )
  
  Cmod <- compileNimble(mod)
  mod.Conf <- configureMCMC(mod, enableWAIC = FALSE)
  mod.MCMC <- buildMCMC(mod.Conf, useConjugacy = FALSE)
  Cmod.MCMC <- compileNimble(mod.MCMC, project = mod)
  
  mod.out <- runMCMC(
    mcmc = Cmod.MCMC,
    niter = 2500,
    nburnin = 1000,
    nchains = 2
  )
  
  id_simul <- floor(runif(1, 1,100000))
  dir.create(paste0("simulation_study/out/_sim_", id_simul))
  saveRDS(mod.out, file = paste0("simulation_study/out/sim_",id_simul, "/MCMC_samples.rds"))
  
}
