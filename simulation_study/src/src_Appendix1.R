library(purrr)
library(nimble)
library(terra)

rm(list = ls())

nX = 50
N = nX **2
T = 4
K = 4

source("functions/utilities.R")
source("functions/simulate_PBDOM.R")
source("functions/load_chains.R")

source("functions/nimbleModel_simulations_indepTrials.R")
nimbleModIndepTrials <- myModel
initIndepTrials <- initial.values

source("functions/nimbleModel_simulations.R")
nimbleModDisp <- myModel
initDisp <- initial.values

res <- 10
sigma <- 20

r <- rast(xmin = 0, xmax = nX*res, ymin = 0, ymax = nX*res, resolution = c(res, res))
d <- as.matrix(dist(terra::crds(r)))

X <- c(scale(MASS::mvrnorm(n = 1, mu = rep(0,length(terra::values(r))),
                         Sigma =  1 / sigma**2 * exp(-d/(2 * sigma**2)))))

params <- data.frame(omega = 0.33, 
                     rho = 0.66, 
                     lambda = 5,
                     sigma = 20,
                     alpha_ksi = logit(0.5),
                     beta_ksi = 0,
                     alpha_psi = -5,
                     beta_psi = 3)

simul_indepTrials <-  pmap(params, simulate_PBDOM, N = N, T = T, d = d, K = K,
                           X = X, A = 100, model = "indepTrials", initZ = "3")
simul_disp <-  pmap(params, simulate_PBDOM, N = N, T = T, d = d, K = K,
                    X = X, A = 100, initZ = "3")

id_simul <- floor(runif(1, 1,100000))
dir.create(paste0("simulation_study/out/Appendix1_sim_", id_simul))

dsparse <- getSparse(d, thr = 60)

### fit1 -----------------------------------------------------------------------

mod <- nimbleModel(
  code = nimbleModIndepTrials,
  data = list(y = simul_indepTrials[[1]]$y),
  constants = list(nSites = N,
                   nSeasons = T,
                   nSurveys = 4,
                   dmatP = dsparse$p,
                   dmatI = unname(dsparse$i),
                   d2 = dsparse$d**2,
                   X = X),
  inits = initIndepTrials(zst = array(1, dim = c(T, N))),
  calculate = FALSE)

Cmod <- compileNimble(mod)
mod.Conf <- configureMCMC(mod, enableWAIC = FALSE)
mod.MCMC <- buildMCMC(mod.Conf, useConjugacy = FALSE)
Cmod.MCMC <- compileNimble(mod.MCMC, project = mod)

out <- runMCMC(
  mcmc = Cmod.MCMC,
  niter = 5000,
  nburnin = 2500,
  nchains = 2)

saveRDS(out, file = paste0("simulation_study/out/Appendix1_sim_", id_simul,
                           "/MCMC_samples_simIndepTrials_fitIndepTrials.rds"))

rm(out)
gc()

### fit2 -----------------------------------------------------------------------

mod <- nimbleModel(
  code = nimbleModDisp,
  data = list(y = simul_indepTrials[[1]]$y),
  constants = list(nSites = N,
                   nSeasons = T,
                   nSurveys = 4,
                   dmatP = dsparse$p,
                   dmatI = unname(dsparse$i),
                   d2 = dsparse$d**2,
                   X = X,
                   pi = pi,
                   A = 100),
  inits = initDisp(zst = array(1, dim = c(T, N))),
  calculate = FALSE)

Cmod <- compileNimble(mod)
mod.Conf <- configureMCMC(mod, enableWAIC = FALSE)
mod.MCMC <- buildMCMC(mod.Conf, useConjugacy = FALSE)
Cmod.MCMC <- compileNimble(mod.MCMC, project = mod)

out <- runMCMC(
  mcmc = Cmod.MCMC,
  niter = 7500,
  nburnin = 5000,
  nchains = 2)

saveRDS(out, file = paste0("simulation_study/out/Appendix1_sim_", id_simul,
                           "/MCMC_samples_simIndepTrials_fitDisp.rds"))

rm(out)
gc()

### fit3 -----------------------------------------------------------------------

mod <- nimbleModel(
  code = nimbleModIndepTrials,
  data = list(y = simul_disp[[1]]$y),
  constants = list(nSites = N, 
                   nSeasons = T,
                   nSurveys = 4, 
                   dmatP = dsparse$p,
                   dmatI = unname(dsparse$i),
                   d2 = dsparse$d**2,
                   X = X),
  inits = initIndepTrials(zst = array(1, dim = c(T, N))),
  calculate = FALSE)

Cmod <- compileNimble(mod)
mod.Conf <- configureMCMC(mod, enableWAIC = FALSE)
mod.MCMC <- buildMCMC(mod.Conf, useConjugacy = FALSE)
Cmod.MCMC <- compileNimble(mod.MCMC, project = mod)

out <- runMCMC(
  mcmc = Cmod.MCMC,
  niter = 5000,
  nburnin = 2500,
  nchains = 2)

saveRDS(out, file = paste0("simulation_study/out/Appendix1_sim_", id_simul,
                           "/MCMC_samples_simDisp_fitIndepTrials.rds"))

rm(out)
gc()
### fit3 -----------------------------------------------------------------------

mod <- nimbleModel(
  code = nimbleModDisp,
  data = list(y = simul_disp[[1]]$y),
  constants = list(nSites = N, 
                   nSeasons = T,
                   nSurveys = 4, 
                   dmatP = dsparse$p,
                   dmatI = unname(dsparse$i),
                   d2 = dsparse$d**2,
                   X = X,
                   pi = pi,
                   A = 100),
  inits = initDisp(zst = array(1, dim = c(T, N))),
  calculate = FALSE)

Cmod <- compileNimble(mod)
mod.Conf <- configureMCMC(mod, enableWAIC = FALSE)
mod.MCMC <- buildMCMC(mod.Conf, useConjugacy = FALSE)
Cmod.MCMC <- compileNimble(mod.MCMC, project = mod)

out <- runMCMC(
  mcmc = Cmod.MCMC,
  niter = 7500,
  nburnin = 5000,
  nchains = 2)

saveRDS(out, file = paste0("simulation_study/out/Appendix1_sim_", id_simul,
                           "/MCMC_samples_simDisp_fitDisp.rds"))
rm(out)
gc()

##### plot estimates -----------------------------------------------------------

# true.vals <- data.frame(par = c("lambda", "alpha_ksi", "alpha_omega",
#                                 "alpha_psi", "alpha_rho", "beta_ksi", "sigma"),
#                         true.val = c(5, logit(0.5), logit(.33), NA , logit(.66), 0, 20))
# 
# sim_paths <- list.files("simulation_study/out", full.names = TRUE, pattern = "Appendix1") %>%
#   list.files(full.names = TRUE)
# 
# mod_summ <- map(sim_paths, readRDS) %>%
#   map(reduce, rbind) %>%
#   map(apply, 2, my_quantile) %>%
#   map(t) %>%
#   map2(sim_paths, ~ mutate(as.data.frame(.x), path = .y)) %>%
#   map(rownames_to_column, var = "par") %>%
#   reduce(rbind)
# 
# mod_summ$simID <- as.numeric(as.factor(reduce(map(str_split(mod_summ$path, "/"), 3), c)))
# mod_summ$sim <- ifelse(grepl("simDisp", mod_summ$path), "disp", "indepTrials") 
# mod_summ$fit <- ifelse(grepl("fitDisp", mod_summ$path), "disp", "indepTrials") 
# mod_summ$match <- mod_summ$fit == mod_summ$sim
# 
# mod_summ <- select(mod_summ, -path)
# mod_summ <- left_join(mod_summ, true.vals)
# 
# ggplot(mod_summ) +
#   geom_pointrange(aes(x=med, xmin = inf, xmax = sup, y = par,
#                       group = simID, pch = paste0("sim:", sim, " - fit:", fit)),
#                   position = position_dodge2(width = 1)) +
#   geom_vline(aes(xintercept = true.val), linetype = "dashed", color = "red") +
#   scale_shape_manual(values = c(16, 2, 1, 17), name = "model")+
#   facet_wrap(~ par, scales = "free") +
#   xlab("") + ylab("") +
#   theme_bw() +
#   theme(text = element_text(face = "bold"),
#         axis.text.y = element_blank(),
#         strip.placement = "outside",
#         strip.background = element_blank())
# 
##### plot kernels -----------------------------------------------------------
#
# A <- 100
# sigma <- 20
# x <- seq(1,60,0.05)
# 
# delta <- function(sigma){
#   function(x){
#     exp(-x**2/(2*sigma**2))
#   }
# }
# 
# delta_<- function(sigma, A, alpha){
#   function(x){
#     1 - exp(-A*alpha/(2*pi*sigma**2)*exp(-x**2/(2*sigma**2)))
#   }
# }
# 
# plot(x, delta(sigma)(x), type = "l", xlab = "d(i,j)", ylab = "delta(i,j)")
# lines(x, delta_(sigma, A, alpha = 5)(x), col = "green")
# lines(x, delta_(sigma, A, alpha = 10)(x), col = "blue")
# lines(x, delta_(sigma, A, alpha = 25)(x), col = "orange")
# lines(x, delta_(sigma, A, alpha = 50)(x), col = "red")
# legend("topright",
#        legend = c("indep-trial", "lambda = 50", "lambda = 25", "lambda = 10",
#                   "lambda = 5"),
#        col = c("black", "red", "orange", "blue", "green"),
#        lwd =1,
#        cex = 0.8,
#        text.col = "black",
#        horiz = F ,
#        inset = c(0.05, 0.05))
