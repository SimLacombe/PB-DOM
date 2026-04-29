library(purrr)
library(sf)
library(tidyverse)
library(nimble)
library(terra)

rm(list = ls())

source("functions/utilities.R")
source("functions/nimbleModel_wolfCaseStudy.R")

filenames <- list(y = "wolf_case_study/data/wolf_louvrier.rds", 
                  effort = "wolf_case_study/data/effort_louvrier.rds",
                  envCovs = "wolf_case_study/data/envcov_louvrier.rds",
                  coords = "wolf_case_study/data/coordcells_wolf.rds",
                  france = "wolf_case_study/data/countries_shp/")

RJMCMC <- TRUE

col_covs <- c("p_forest", "p_agri", "p_alti", "p_dbarr", "p_halt")
  
### 1. Get data and covariates -------------------------------------------------

coords <- readRDS(filenames$coords) %>% 
  st_as_sf(coords = c('X','Y'), crs = 2154) %>%
  select()

sites_kept <- !sapply(st_intersects(coords, filter(read_sf(filenames$france), NAME == "France")),
                      is_empty)

coords <- coords [sites_kept,]
dat <- readRDS(filenames$y)[sites_kept,]

Nsites <- nrow(dat)
Nyears <- 23
Nseasons <- 4

dat <- map(1:Nsites, function(i){matrix(dat[i,], nrow = Nseasons, ncol = Nyears)}) %>%
  simplify2array() %>%
  aperm(c(2,3,1))

y <- apply(dat, c(1,2), function(x){sum(sign(x), na.rm = T)})

effort <- readRDS(filenames$effort)[sites_kept,]

envCovs <- readRDS(filenames$envCovs)[sites_kept,]

X_det  <- list(t(effort),
               matrix(envCovs[, "p_road"], 23, Nsites, byrow = TRUE)) %>%
  simplify2array()
  
### 2. Get distance matrix -----------------------------------------------------

dmat <- st_distance(coords)
dmat <- matrix(as.numeric(dmat)/ 1000, Nsites, Nsites)

sparsemat <- getSparse(dmat, thr = 50)

### 3. Set up model ------------------------------------------------------------

myConstants <- list(nSites = Nsites, 
                    nSeasons = Nyears,
                    nSurveys = Nseasons, 
                    ncovs_col = length(col_covs),
                    ncovs_det = 2,
                    X_col = envCovs[, col_covs],
                    X_det = X_det,
                    dmatP = sparsemat$p,
                    dmatI = sparsemat$i,
                    d2 = sparsemat$d**2,
                    A = 100,
                    pi = pi,
                    RJMCMC = RJMCMC)

mod <- nimbleModel(
  code = myModel,
  data = list(y = y),
  constants = myConstants,
  inits = initial.values(zst = array(1, dim = c(myConstants$nSeasons, myConstants$nSites)),
                         ncovs_col = myConstants$ncovs_col, ncovs_det = myConstants$ncovs_det),
  calculate = FALSE
)

Cmod <- compileNimble(mod)
mod.Conf <- configureMCMC(mod, enableWAIC = FALSE)
mod.Conf$addMonitors("z")

if(RJMCMC){
  mod.Conf$addMonitors("rj")
  
  configureRJ(conf = mod.Conf,
              targetNodes = "beta_ksi",
              indicatorNodes = "rj",
              control = list(mean = 0, scale = 2))
}

mod.MCMC <- buildMCMC(mod.Conf, useConjugacy = FALSE)
Cmod.MCMC <- compileNimble(mod.MCMC, project = mod)

### Burn in --------------------------------------------------------------------

chain <- floor(runif(1, 100000, 999999)) #get Unique identifier for the MCMC chain
cat(chain, "\n")
dir.create(paste0("wolf_case_study/out/CHAIN_", chain))

Cmod.MCMC$run(niter = 10000, nburnin = 10000, reset = TRUE)

### Sample ---------------------------------------------------------------------

nChunks <- 10
nIter <- 500

for(c in 1:nChunks){
  cat("Running chunk", c, "of", nChunks, "\n")
  Cmod.MCMC$run(niter = nIter, nburnin = 0, reset = FALSE, resetMV = TRUE)
  samples <- as.matrix(Cmod.MCMC$mvSamples)
  saveRDS(samples,
          file = paste0("wolf_case_study/out/CHAIN_", chain, "/MCMC_samples_wolf_model_chunk", c, ".rds"))
  rm(samples)
  gc()
}
