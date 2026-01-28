library(tidyr)
library(dplyr)
library(purrr)
library(sf)
library(nimble)

rm(list = ls())

source("functions/utilities.R")

filenames <- list(loutres = "otter_case_study/data/v_ref_presence_loutre_par_an_2000_25_utf8.csv",
                  bg = "otter_case_study/data/background_grid.gpkg",
                  eff = "otter_case_study/data/effort.rds")

annees <- 2000:2023
thr <- 50

gamma_covs <- c("riverLength","ripForest", "CYP", "ANG", "ECR", "SAL")

### 1. DONNEES  ----------------------------------------------------------------

donnees_loutres <- read.csv(filenames$loutres, sep = ";") %>%
  filter(annee %in% annees) %>%
  rename(n = nombre.de.données.dans.la.maille)
donnees_loutres <- st_as_sf(donnees_loutres, wkt = "wkt_geom_4326", crs = 4326) %>%
  st_transform(crs = 2154)

bg <- read_sf(filenames$bg) %>%
  filter(cellArea >= 3.33e07)

effort <- readRDS(filenames$eff)

### 2. FORMATAGE DES DONNEES ---------------------------------------------------

K <- 4
N <- nrow(bg)
T <- length(annees) / K

y_df <- split(st_centroid(donnees_loutres), donnees_loutres$annee) %>%
  map2(annees, ~ mutate(st_join(bg, .x), annee = .y)) %>%
  reduce(rbind) %>%
  st_drop_geometry %>%
  mutate(t = (annee - annees[1]) %/% K + 1,
         k = (annee - annees[1]) %% K + 1,
         y = sign(replace_na(n, 0))) %>%
  select(id, t, k, y) %>%
  arrange(id)

y <- array(y_df$y, dim = c(K, T, N)) %>%
  aperm(c(2, 3, 1))

d <- st_distance(st_centroid(bg))
d <- matrix(as.numeric(d)/1000, N, N)
d_sparse <- getSparse(d, thr = thr)
rm(d)

bg <- bg %>% 
  mutate(across(all_of(c("ANG", "CYP", "ECR", "SAL")),
                ~ log(.x + 1))) %>%
  mutate(across(all_of(gamma_covs), ~ c(scale(.x)))) %>%
  arrange(id)

yearMat <- matrix(1:(T*K), T, K, byrow = TRUE)
  
### 3. Fit NIMBLE --------------------------------------------------------------

source("functions/nimbleModel_otterCaseStudy.R")

chain <- floor(runif(1, 100000, 999999)) #get Unique identifier for the MCMC chain
cat(chain, "\n")
dir.create(paste0("otter_case_study/out/CHAIN_", chain))

my.constants <- list(nSites = N, 
                     nSeasons = T,
                     nSurveys = K, 
                     ncovs_gamma = length(gamma_covs),
                     X = as.matrix(st_drop_geometry(bg[,gamma_covs])),
                     year = yearMat,
                     effort = effort,
                     dmatP = d_sparse$p,
                     dmatI = d_sparse$i ,
                     d2 = d_sparse$d**2,
                     pi = pi)

mod <- nimbleModel(
  code = myModel,
  data = list(y = y),
  constants = my.constants,
  inits = initial.values(zst = array(1, dim = c(T, N)),
                         ncovs_gamma = my.constants$ncovs_gamma),
  calculate = FALSE
)

Cmod <- compileNimble(mod)
mod.Conf <- configureMCMC(mod, enableWAIC = FALSE)
mod.Conf$addMonitors("z")

mod.MCMC <- buildMCMC(mod.Conf, useConjugacy = FALSE)
Cmod.MCMC <- compileNimble(mod.MCMC, project = mod)

Cmod.MCMC$run(niter = 5000, nburnin = 5000, reset = TRUE)

nChunks <- 5
nIter <- 500
for(c in 1:nChunks){
  cat("Running chunk", c, "of", nChunks, "\n")
  Cmod.MCMC$run(niter = nIter, nburnin = 0, reset = FALSE, resetMV = TRUE)
  samples <- as.matrix(Cmod.MCMC$mvSamples)
  saveRDS(samples,
          file = paste0("otter_case_study/out/CHAIN_", chain, "/MCMC_samples_otter_model_chunk", c, ".rds"))
  rm(samples)
  gc()
}