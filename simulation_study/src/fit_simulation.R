library(purrr)
library(nimble)

rm(list = ls())

source("functions/utilities.R")
source("functions/nimbleModel_simulations.R")
load("simulation_study/data/distanceMatrix.Rdata")

path <- list.dirs(list.dirs("simulation_study/out", recursive = FALSE), recursive = FALSE)
path <- path[!file.exists(paste0(path, "/done.txt"))][1]
cat(path)
cat("done", file = paste0(path,"/done.txt"))

dat <- readRDS(paste0(path,"/data.rds"))

for(t in 2:4){
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
  
  saveRDS(mod.out, file = paste0(path, "/MCMC_samples_thr", t, "-sigma.rds"))
          
  # png(paste0(path, "/MCMC_traceplot_", t,"sigma.png"), 
  #     width=21, height=21, units = "cm", res = 300)
  # MCMCvis::MCMCtrace(
  #   Rhat = TRUE,
  #   n.eff = TRUE,
  #   object = mod.out,
  #   pdf = FALSE,
  #   ind = TRUE)
  # dev.off()
  
}
