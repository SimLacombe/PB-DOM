library(purrr)

rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)
id_simul <-ifelse(is_empty(args),  floor(runif(1,0,100000)), as.integer(args[1]))

source("functions/utilities.R")
source("functions/simulate_PBDOM.R")

# set.seed(9553)

nX <- 40
N <- nX * nX
T <- 5
K <- 4

r <- terra::rast(xmin = 0, xmax = nX * 10, ymin = 0, ymax = nX * 10, resolution = c(10, 10))
d <- as.matrix(dist(terra::crds(r)))

X <- r
terra::values(X) <- c(scale(MASS::mvrnorm(n = 1, mu = rep(0,length(terra::values(r))),
                                          Sigma =  1 / 7.5**2 * exp(-d/(2 * 7.5**2)))))

params <- data.frame(omega = 0.2, rho = 0.66,
                     lambda = 5,
                     sigma = c(12.5,25,12.5,25),
                     alpha_ksi = logit(c(0.50,0.50, 0.5,0.5)),
                     beta_ksi = c(0,0,5,5))

simul <-  pmap(params, simulate_PBDOM, 
               N = N, T = T, d = d, K = K, X = c(scale(terra::values(X))), A = 100)

names(simul) <- paste0("Scenario", 1:nrow(params))

### SAVE -----------------------------------------------------------------------

dir.create(paste0("simulation_study/out/sim_", id_simul))
walk(names(simul), ~ dir.create(paste0("simulation_study/out/sim_", id_simul, "/", .x)))
iwalk(simul, ~ saveRDS(.x, file=paste0("simulation_study/out/sim_", id_simul,"/", .y, "/data.rds")))
 
### ROUTINE PLOTS --------------------------------------------------------------

png(paste0("simulation_study/out/sim_", id_simul, "/X.png"),
    width = 1200, height = 960, units = "px")
terra::plot(X)
dev.off()

png(paste0("simulation_study/out/sim_", id_simul, "/Z.png"),
    width = 1200, height = 960, units = "px")
par(mfrow = (c(nrow(params),T)))
pwalk(expand.grid(t  = 1:T, s = 1:nrow(params)), function(t, s){
  r <- terra::rast(xmin = 0, xmax = 300, ymin = 0, ymax = 300, resolution = c(10, 10))
  terra::values(r) <- simul[[s]]$z[t, ]
  terra::plot(r, legend = FALSE, mar = c(.1,.1,.1,.1), axes = NA, col = c("lightgrey", "chartreuse4"))})
dev.off()

png(paste0("simulation_study/out/sim_", id_simul, "/gamma.png"),
    width = 1200, height = 960, units = "px")
par(mfrow = (c(nrow(params),T-1)))
pwalk(expand.grid(t  = 1:(T-1), s = 1:nrow(params)), function(t, s){
  r <- terra::rast(xmin = 0, xmax = 300, ymin = 0, ymax = 300, resolution = c(10, 10))
  terra::values(r) <- simul[[s]]$gamma[t, ]
  terra::plot(r, legend = FALSE)})
dev.off()

### FIG 1 ----------------------------------------------------------------------

# library(tidyverse)
# 
# X.df <- as.data.frame(X, xy = TRUE)
# X.df <- map(simul, ~ cbind(X.df, t(.x$z))) %>%
#   imap(~ mutate(.x, scenario = .y)) %>%
#   reduce(rbind)
# 
# colnames(X.df) <- c("x", "y", "X", paste0("t=", 1:5), "scenario")
# X.df <- pivot_longer(X.df, cols = all_of(paste0("t=", 1:5)), names_to = "t", values_to = "z")
# 
# ggplot(X.df) + 
#   geom_tile(aes(x=x,y=y, fill = X), 
#             col = NA) +
#   geom_tile(data = X.df %>% filter(z == 1), 
#             aes(x=x,y=y), col = "black",fill = "darkblue", alpha = .5) +
#   paletteer::scale_fill_paletteer_c("grDevices::RdYlGn")+
#   facet_grid(scenario~t) +
#   theme_void() +
#   theme(text = element_blank()) +
#   coord_fixed()
# 
# ggsave("figures/simuls1.png", height = 1920, width = 2560, units = "px")
# 
# x <- seq(-2,2,0.1)
# y <- invlogit(5*x)
# y2 <- 0.80
# d <- seq(0,75,1)
# Lambda1 <- 100 * 5 / (2 * pi * 12.5**2) * exp(-d**2/(2*12.5**2))
# Lambda2 <- 100 * 5 / (2 * pi * 25**2) * exp(-d**2/(2*25**2))
# 
# plot1 <- rbind(
#   data.frame(d = d, y = Lambda1, scenario = "Short dispersal distance", panel = "Dispersal distance"),
#   data.frame(d = d, y = Lambda2, scenario = "Long dispersal distance", panel = "Dispersal distance")) %>%
#   ggplot() +
#   geom_line(aes(x=d, y=y, col = scenario), linewidth = .75) +
#   xlab("Dispersal distance (km)") + ylab("Dispersal pressure (\u039b)")+
#   theme_minimal() +
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         text = element_text(face = "bold"))
# plot2 <- rbind(
#   data.frame(x = x, y = y, scenario = "Habitat specialist", panel = "Specialist VS generalist"),
#   data.frame(x = x, y = y2, scenario = "Habitat generalist", panel = "Specialist VS generalist")) %>%
#   ggplot() +
#     geom_line(aes(x=x, y=y, col = scenario), linewidth = .75) +
#     xlab("Scaled covariate value") + ylab("Installation probability (\u03be)")+
#     theme_minimal() +
#     theme(legend.position = "bottom",
#           legend.title = element_blank(),
#           text = element_text(face = "bold"))
# 
# 
# ggpubr::ggarrange(plot1, plot2, nrow = 2)
# ggsave("figures/simuls2.png", width = 1280, height = 1920, units = "px")
