require(nimble)

myModel <- nimbleCode({
  
  ### PRIORS ###
  
  alpha_psi  ~ dnorm(0, 1/2.5^2)   # initial occupancy intercept
  alpha_gam  ~ dnorm(0, 1/2.5^2)   # colonization intercept
  alpha_omega~ dnorm(0, 1/2.5^2)   # extinction intercept
  alpha_rho  ~ dnorm(0, 1/2.5^2)   # detection intercept
  alpha_ginf ~ dnorm(logit(0.05), 1) # long distance colonization probability
  
  alpha  ~ dlnorm(log(1), 1)
  sigma ~ dlnorm(log(30), 1)
  
  for (c in 1:ncovs_col){
    beta_gam[c]  ~ dnorm(0, 1/2.5^2)   
  }   
  for (c in 1:ncovs_det){
    beta_rho[c]  ~ dnorm(0, 1/2.5^2)   
  }   
  
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
      gamma[t - 1, i] <- gamma0[i] * (1 - exp(- (sum(kappa[t - 1, (dmatP[i] + 1):dmatP[i + 1]]) + gam_inf)))
      
      z[t, i] ~ dbern(z[t - 1, i] * (1 - omega) +
                        (1 - z[t - 1, i]) * gamma[t - 1, i])
    }
  }
  
  ### OBSERVATION MODEL ###
  
  for (i in 1:nSites) {
    for (t in 1:nSeasons) {
      y[t, i] ~ dbinom(size = nSurveys, prob = z[t, i] * rho[t,i])
      
    }
  }
  
  ### LOGIT LINEAR PREDICTORS ###
  
  omega <- ilogit(alpha_omega)
  psi <- ilogit(alpha_psi)
  gam_inf <- ilogit(alpha_ginf)
  
  for(i in 1:nSites){
    gamma0[i] <- ilogit(alpha_gam + inprod(beta_gam[1:ncovs_col], X_col[i, 1:ncovs_col]))
    for(t in 1:nSeasons){
      rho[t, i] <- ilogit(alpha_rho + inprod(beta_rho[1:ncovs_det], X_det[t, i, 1:ncovs_det]))
    }
  }
  
  
})

initial.values <- function(zst, ncovs_col, ncovs_det) {
  list(
    z = zst,
    alpha_psi = rnorm(1, 0, 1/2.5**2),
    alpha_gam = rnorm(1, 0, 1/2.5**2),
    alpha_omega = rnorm(1, 0, 1/2.5**2),
    alpha_rho = rnorm(1, 0, 1/2.5**2),
    beta_gam = rnorm(ncovs_col, 0, 1/2.5**2),
    beta_rho = rnorm(ncovs_det, 0, 1/2.5**2),
    alpha_ginf = rnorm(1, logit(0.05), 1), 
    sigma = rlnorm(1, log(30), 1)
  )
}
