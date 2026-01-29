require(nimble)

myModel <- nimbleCode({
  
  ### PRIORS ###
  
  # Intercepts
  alpha_gam   ~ dnorm(0, 1/2.5^2)         # colonization intercept
  alpha_psi   ~ dnorm(0, 1/2.5^2)         # initial occu intercept
  alpha_omega ~ dnorm(0, 1/2.5^2)         # extinction intercept
  alpha_rho[1] ~ dnorm(0, 1/2.5^2)       # detection intercept for opport visits
  alpha_rho[2] ~ dnorm(0, 1/2.5^2)       # detection intercept for stdzed data
  
  # Slopes 
  for(c in 1:ncovs_gamma){
    beta[c] ~ dnorm(0, 1/2.5^2)
  }
  beta_rho_t[1] ~ dnorm(0, 1/2.5^2)  # detection time trend for opport visits
  beta_rho_t[2] ~ dnorm(0, 1/2.5^2)  # detection time trend for stdzed visits
  
  # Scale and intensity parameters of the colonisation kernel 
  sigma ~ dlnorm(log(30), 1)       
  alpha ~ dlnorm(log(0.5), 1)
  
  ### ECOLOGICAL MODEL ###
  
  # loop over sites
  for (i in 1:nSites) {
    # first year
    z[1, i] ~ dbern(psi)
    
    # following years
    for (t in 2:nSeasons) {
        # loop over the neighbors of i
        for (j in (dmatP[i] + 1):dmatP[i + 1]) {
          delta[t - 1, j] <- exp(-d2[j] / (2 * sigma ** 2)) * z[t - 1, dmatI[j]]
        }
        gamma_ti[t - 1, i] <- gamma[i] * (1 - exp(- 100 * alpha/ (2 * pi * sigma ** 2) * sum(delta[t - 1, (dmatP[i] + 1):dmatP[i + 1]])))
        z[t, i] ~ dbern(z[t - 1, i] * (1 - omega) + (1 - z[t - 1, i]) * gamma_ti[t - 1, i])
      }
    }
  
  ### OBSERVATION MODEL ###
  
  for (i in 1:nSites) {
    for (t in 1:nSeasons) {
      for (k in 1:nSurveys){
        y[t, i, k] ~ dbern(prob = z[t, i] * rho[t, i, k])
      }
    }
  }
  
  ### LOGIT LINEAR PREDICTORS ###
  
  omega <- ilogit(alpha_omega)
  psi <- ilogit(alpha_psi)

  for (i in 1:nSites) {
    gamma[i] <- ilogit(alpha_gam + inprod(beta[1:ncovs_gamma], X[i, 1:ncovs_gamma]))
    for (t in 1:nSeasons) {
      for (k in 1:nSurveys){
        rho[t, i, k] <- ilogit(alpha_rho[effort[year[t,k], i] + 1] +
                                 beta_rho_t[effort[year[t,k], i] + 1] * (year[t,k] - 12))
      }
    }
  }
})

initial.values <- function(zst, ncovs_gamma) {
  list(
    z = zst,
    alpha_psi = rnorm(1, 0, 1/2.5**2),
    alpha_gam = rnorm(1, 0, 1/2.5**2),
    alpha_omega = rnorm(1, 0, 1/2.5**2),
    alpha_rho = rnorm(2, 0, 1/2.5**2), 
    beta = rnorm(ncovs_gamma, 0, 1/2.5**2),
    beta_rho_t = rnorm(2, 0, 1/2.5**2),
    sigma = rlnorm(1, log(30), 1),
    alpha = rlnorm(1, log(0.2), 1)
  )
}
