library(purrr)
library(nimble)

rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
id_simul <-as.integer(args[1])
id_scenar <-as.integer(args[2])

source("functions/utilities.R")
source("functions/nimbleModel_simulations.R")
load("simulation_study/data/distanceMatrix.Rdata")

dat <- readRDS(paste0("simulation_study/out/sim_",id_simul, "/Scenario", id_scenar, "/data.rds"))

for(t in 2:4){
  cat("\n", "RUNNING simulation ", id_simul, " - Scenario ", id_scenar, " - thr ", t, " sigma", "\n")
  dsparse <- getSparse(d, thr = t * dat$sigma)
  
  myConstants <- list(nSites = ncol(dat$y), 
                      nSeasons = nrow(dat$y),
                      nSurveys = 4, 
                      dmatP = dsparse$p,
                      dmatI = unname(dsparse$i),
                      d2 = dsparse$d**2,
                      X = dat$X,
                      pi = pi,
                      A = 100)
  
  mod <- nimbleModel(
    code = myModel,
    data = list(y = dat$y),
    constants = myConstants,
    inits = initial.values(zst = array(1, dim = c(myConstants$nSeasons, myConstants$nSites))),
    calculate = FALSE
  )
  
  Cmod <- compileNimble(mod)
  mod.Conf <- configureMCMC(mod, enableWAIC = FALSE)
  mod.MCMC <- buildMCMC(mod.Conf, useConjugacy = FALSE)
  Cmod.MCMC <- compileNimble(mod.MCMC, project = mod)
  
  mod.out <- runMCMC(
    mcmc = Cmod.MCMC,
    niter = 7500,
    nburnin = 5000,
    nchains = 2
  )
  
  saveRDS(mod.out,
          file = paste0("simulation_study/out/sim_",id_simul, "/Scenario",
          id_scenar, "/MCMC_samples_thr-", t, "-sigma.rds"))
  
}
