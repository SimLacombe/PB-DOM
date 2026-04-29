simulate_PBDOM <- function(N, T, K, d, X, A, rho, omega, sigma, lambda, alpha_ksi,
                           beta_ksi, alpha_psi = NULL, beta_psi = NULL, initZ = "1", model = "propagule"){
  
 delta <- exp(-d**2/(2*sigma**2))
 if(model == "propagule") delta <- A * lambda / (2 * pi * sigma **2) * delta
 
 ksi <- invlogit(alpha_ksi + beta_ksi * X)
 
 z <- matrix(0, T, N)
 y <- matrix(0, T, N)
 gamma <- matrix(0, T - 1 , N)
 
z[1, ] <-as.numeric(d[N/2+sqrt(N)/2,] < 50)
if(initZ == "2"){
  z[1, ] <- sign(as.numeric(d[1,] < 100) + rbinom(N, 1, 0.005))
}
 if(initZ == "3"){
   z[1, ] <- rbinom(N, 1, invlogit(alpha_psi + beta_psi * X))
 }
   
 y[1,] <- rbinom(N, K, rho * z[1,])
 
 for(t in 2:T){
   if(model == "propagule") gamma[t-1, ] <-  ksi * apply(delta, 1, function(.x){(1 - exp (- sum(.x*z[t-1,])))})
   if(model == "indepTrials") gamma[t-1, ] <-  ksi * apply(delta, 1, function(.x){(1 - prod(1 - .x*z[t-1,]))})
   z[t, ] <- rbinom(N, 1, z[t-1, ] * (1 - omega) + (1 - z[t-1, ]) * gamma[t-1, ] )
   y[t,] <- rbinom(N, K, rho * z[t,])
 }
 
 return(list(sigma = sigma, 
             d = d,
             X = X,
             y = y,
             z = z,
             gamma = gamma))
}
