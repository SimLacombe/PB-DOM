simulate_PBDOM <- function(N, T, K, d, X, A, rho, omega, sigma, alpha, alpha_gam, beta_gam){
  
 kappa <- A * alpha / (2 * pi * sigma **2) * exp(-d**2/(2*sigma**2))
 gamma0 <- invlogit(alpha_gam + beta_gam * X)
 
 z <- matrix(0, T, N)
 y <- matrix(0, T, N)
 gamma <- matrix(0, T - 1 , N)
 
 # z[1, ] <- as.numeric(d[which(X == sort(X)[floor(0.50 * length(X))]),] < 50)
 z[1, ] <- as.numeric(d[sample(1:nrow(d),1),] < 50)
 y[1,] <- rbinom(N, K, rho * z[1,])
 
 for(t in 2:T){
   gamma[t-1, ] <-  gamma0 * apply(kappa, 1, function(.x){(1 - exp (- sum(.x*z[t-1,])))})
   z[t, ] <- rbinom(N, 1, z[t-1, ] * (1 - omega) + (1 - z[t-1, ]) * gamma[t-1, ] )
   y[t,] <- rbinom(N, K, rho * z[t,])
 }
 
 return(list(sigma = sigma, 
             X = X,
             y = y,
             z = z,
             gamma = gamma))
}
