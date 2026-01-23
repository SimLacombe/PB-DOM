library(purrr)
library(sf)
library(tidyverse)
library(nimble)
library(terra)

rm(list = ls())

source("functions/utilities.R")
source("functions/nimbleModel_wolfeCaseStudy.R")

filenames <- list(y = "wolfe_case_study/data/wolf_louvrier.rds", 
                  effort = "wolfe_case_study/data/effort_louvrier.rds",
                  envCovs = "wolfe_case_study/data/envcov_louvrier.rds",
                  coords = "wolfe_case_study/data/coordcells_wolf.rds")


### 1. Get data and covariates -------------------------------------------------

dat <- readRDS(filenames$y)

Nsites <- nrow(dat)
Nyears <- 23
Nseasons <- 4

dat <- map(1:Nsites, function(i){matrix(dat[i,], nrow = Nseasons, ncol = Nyears)}) %>%
  simplify2array() %>%
  aperm(c(2,3,1))

y <- apply(dat, c(1,2), function(x){sum(sign(x), na.rm = T)})

effort <- readRDS(filenames$effort)

envCovs <- readRDS(filenames$envCovs)

X_det  <- list(t(effort),
               matrix(envCovs[, "p_road"], 23, 3450, byrow = TRUE)) %>%
  simplify2array()
  
### 2. Get distance matrix -----------------------------------------------------

coords <- readRDS(filenames$coords) %>% 
  st_as_sf(coords = c('X','Y'), crs = st_crs(grid)) %>%
  select()

dmat <- st_distance(coords)
dmat <- matrix(as.numeric(dmat)/ 1000, Nsites, Nsites)

sparsemat <- getSparse(dmat, thr = 50)

### 3. Set up model ------------------------------------------------------------

myConstants <- list(nSites = Nsites, 
                    nSeasons = Nyears,
                    nSurveys = Nseasons, 
                    ncovs_col = 4,
                    ncovs_det = 2,
                    X_col = envCovs[, c("p_forest", "p_agri", "p_alti", "p_halt")],
                    X_det = X_det,
                    dmatP = sparsemat$p,
                    dmatI = sparsemat$i,
                    d2 = sparsemat$d**2,
                    A = 100,
                    pi = pi)

mod <- nimbleModel(
  code = myModel,
  data = list(y = y),
  constants = myConstants,
  inits = initial.values(zst = array(1, dim = c(myConstants$nSeasons, myConstants$nSites)),
                         ncovs_col = 4, ncovs_det = 2),
  calculate = FALSE
)

Cmod <- compileNimble(mod)
mod.Conf <- configureMCMC(mod, enableWAIC = FALSE)
mod.Conf$addMonitors("z")
mod.MCMC <- buildMCMC(mod.Conf, useConjugacy = FALSE)
Cmod.MCMC <- compileNimble(mod.MCMC, project = mod)

### Burn in --------------------------------------------------------------------

chain <- floor(runif(1, 100000, 999999)) #get Unique identifier for the MCMC chain
Cmod.MCMC$run(niter = 5000, nburnin = 5000, reset = TRUE)

### Sample ---------------------------------------------------------------------
nChunks <- 10
nIter <- 500
for(c in 1:nChunks){
  cat("Running chunk", c, "of", nChunks, "\n")
  Cmod.MCMC$run(niter = nIter, nburnin = 0, reset = FALSE, resetMV = TRUE)
  samples <- as.matrix(Cmod.MCMC$mvSamples)
  saveRDS(samples,
          file = paste0("wolfe_case_study/out/wolfModelFitted_chain", chain, "_chunk", c, ".rds"))
  rm(samples)
  gc()
}
