require(nimble)

myModel <- nimbleCode({
  
  ### PRIORS ###
  
  # Intercepts
  
  alpha_psi   ~ dnorm(0, 1/2.5^2)         # initial occu intercept
  alpha_ksi   ~ dnorm(0, 1/2.5^2)         # colonization intercept
  alpha_omega ~ dnorm(0, 1/2.5^2)         # extinction intercept
  alpha_rho[1] ~ dnorm(0, 1/2.5^2)        # detection intercept for opport visits
  alpha_rho[2] ~ dnorm(0, 1/2.5^2)        # detection intercept for stdzed data
  alpha_lam_min ~ dnorm(logit(0.05), 1) # long distance colonization probability
  
  # Slopes 
  p_accept ~ dbeta(1, 1)
  for(c in 1:ncovs_col){
    beta_ksi[c] ~ dnorm(0, 1/2.5^2)
    if(RJMCMC) rj[c] ~ dbern(p_accept)
    if(!RJMCMC) rj[c] <- 1
    beta_ksi_eff[c] <- beta_ksi[c] * rj[c] 
  }
  
  beta_rho_t[1] ~ dnorm(0, 1/2.5^2)  # detection time trend for opport visits
  beta_rho_t[2] ~ dnorm(0, 1/2.5^2)  # detection time trend for stdzed visits

  # Scale and intensity parameters of the colonisation kernel 
  sigma ~ dlnorm(log(30), 1)       
  lambda ~ dlnorm(log(0.5), 1)
  
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
        gamma[t - 1, i] <- ksi[i] * (1 - exp(- (lam_min + 100 * lambda/ (2 * pi * sigma ** 2) * sum(delta[t - 1, (dmatP[i] + 1):dmatP[i + 1]]))))
        z[t, i] ~ dbern(z[t - 1, i] * (1 - omega) + (1 - z[t - 1, i]) * gamma[t - 1, i])
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
  psi   <- ilogit(alpha_psi)
  lam_min <- ilogit(alpha_lam_min)
  for (i in 1:nSites) {
    ksi[i] <- ilogit(alpha_ksi + inprod(beta_ksi_eff[1:ncovs_col], X[i, 1:ncovs_col]))
    for (t in 1:nSeasons) {
      for (k in 1:nSurveys){
        rho[t, i, k] <- ilogit(alpha_rho[effort[year[t,k], i] + 1] + beta_rho_t[effort[year[t,k], i] + 1] * (year[t,k] - 12))
      }
    }
  }
})

initial.values <- function(zst, ncovs_col) {
  list(
    z = zst,
    alpha_psi = rnorm(1, 0, 1/2.5**2),
    alpha_ksi = rnorm(1, 0, 1/2.5**2),
    alpha_omega = rnorm(1, 0, 1/2.5**2),
    alpha_rho = rnorm(2, 0, 1/2.5**2), 
    beta_ksi = rnorm(ncovs_col, 0, 1/2.5**2),
    beta_rho_t = rnorm(2, 0, 1/2.5**2),
    alpha_lam_min = rnorm(1, logit(0.05), 1), 
    sigma = rlnorm(1, log(30), 1),
    alpha = rlnorm(1, log(0.2), 1),
    rj = rep(1, ncovs_col),
    p_accept = 1
  )
}
