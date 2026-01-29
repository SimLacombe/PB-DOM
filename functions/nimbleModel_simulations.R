myModel <- nimbleCode({
  
  ### PRIORS ###
  
  alpha_psi   ~ dnorm(0, 1/2.5^2)   # initial occupancy intercept
  alpha_gam   ~ dnorm(0, 1/2.5^2)   # colonization intercept
  beta_gam    ~ dnorm(0, 1/2.5^2)   # colonization intercept
  alpha_omega ~ dnorm(0, 1/2.5^2)   # extinction intercept
  alpha_rho   ~ dnorm(0, 1/2.5^2)   # detection intercept
  
  alpha  ~ dlnorm(log(1), 1)
  sigma ~ dlnorm(log(30), 1)
  
  ### ECOLOGICAL MODEL ###
  
  # loop over sites
  for (i in 1:nSites) {
    
    # first year
    z[1, i] ~ dbern(psi)
    
    # following years
    for (t in 2:nSeasons) {
      # loop over the neighbors of i
      for (j in (dmatP[i] + 1):dmatP[i + 1]) {
        kappa[t - 1, j] <- A * alpha / (2 * pi * sigma **2) * exp(-d2[j] / (2 * sigma ** 2)) * z[t - 1, dmatI[j]]
        }
      gamma[t - 1, i] <- gamma0[i] * (1 - exp(- sum(kappa[t - 1, (dmatP[i] + 1):dmatP[i + 1]])))
      z[t, i] ~ dbern(z[t - 1, i] * (1 - omega) +
                        (1 - z[t - 1, i]) * gamma[t - 1, i])
    }
  }
  
  ### OBSERVATION MODEL ###
  
  for (i in 1:nSites) {
    for (t in 1:nSeasons) {
      y[t, i] ~ dbinom(size = nSurveys, prob = z[t, i] * rho)
      
    }
  }
  
  ### LOGIT LINEAR PREDICTORS ###
  
  omega <- ilogit(alpha_omega)
  rho <- ilogit(alpha_rho)
  psi <- ilogit(alpha_psi)
  for(i in 1:nSites){
    gamma0[i] <- ilogit(alpha_gam + beta_gam * X[i])
  }
})

initial.values <- function(zst) {
  list(
    z = zst,
    alpha_psi = rnorm(1, 0, 1/2.5**2),
    alpha_gam = rnorm(1, 0, 1/2.5**2),
    beta_gam = rnorm(1, 0, 1/2.5**2),
    alpha_omega = rnorm(1, 0, 1/2.5**2),
    alpha_rho = rnorm(1, 0, 1/2.5**2),
    sigma = rlnorm(1, log(30), 1),
    alpha = rlnorm(1, log(1), 1)
  )
}