simulate_PBDOM <- function(N, T, K, d, X, A, rho, omega, sigma, alpha, alpha_gam,
                           beta_gam, alpha_psi = NULL, beta_psi = NULL, model = "propagule"){
  
 kappa <- exp(-d**2/(2*sigma**2))
 if(model == "propagule") kappa <- A * alpha / (2 * pi * sigma **2) * kappa
 gamma0 <- invlogit(alpha_gam + beta_gam * X)
 
 z <- matrix(0, T, N)
 y <- matrix(0, T, N)
 gamma <- matrix(0, T - 1 , N)
 
if(is.null(alpha_psi)) {
  z[1, ] <- sign(as.numeric(d[sample(1:nrow(d),1),] < 50) + rbinom(N, 1, 0.005))
} else {
  z[1, ] <- rbinom(N, 1, invlogit(alpha_psi + beta_psi * X))
}
   
 y[1,] <- rbinom(N, K, rho * z[1,])
 
 for(t in 2:T){
   if(model == "propagule") gamma[t-1, ] <-  gamma0 * apply(kappa, 1, function(.x){(1 - exp (- sum(.x*z[t-1,])))})
   if(model == "indepTrials") gamma[t-1, ] <-  gamma0 * apply(kappa, 1, function(.x){(1 - prod(1 - .x*z[t-1,]))})
   z[t, ] <- rbinom(N, 1, z[t-1, ] * (1 - omega) + (1 - z[t-1, ]) * gamma[t-1, ] )
   y[t,] <- rbinom(N, K, rho * z[t,])
 }
 
 return(list(sigma = sigma, 
             X = X,
             y = y,
             z = z,
             gamma = gamma))
}
