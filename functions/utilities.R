logit <- function(p) log(p/(1-p))
invlogit <- function(x) 1/(1 + exp(-x))
delta <- function(x, lambda, sigma){return(lambda * 100 / (2*pi*sigma**2) * exp(-x**2/(2*sigma**2)))}

getSparse <- function(x, thr){
  # NB: I can replace 0 in the following lines by Inf and (> by !=). It should not change anything with the 
  # present formulation of the model, but if I want to model omega as a function of neighbors the previous 
  # year (e.g. rescue effect) I might need to do it.
  
  x[x > thr] <- 0
  d <- as.numeric(x[x > 0])
  tmp <- apply(x, 1, FUN = function(x){which(x > 0)})
  i <- unlist(tmp)
  p <- cumsum(c(0, reduce(map(tmp, length), c)))
  
  list(d = d,
       i = i,
       p = p)
}

predict_over <- function(x, alpha, beta, link = invlogit){
  cbind(x, t(apply(invlogit(cbind(alpha, beta) %*% t(cbind(1, x))), 2, my_quantile)))
}

predict_over_all <- function(X, alpha, beta, covs){
  if(ncol(X) > 1){
    X <- X[, covs]
    out <- map(1:ncol(X), ~ data.frame(cov = covs[.x],
                                       predict_over(seq(min(X[, .x]), max(X[, .x]),  0.1),
                                                    alpha, beta[, .x]))) %>%
      reduce(rbind)
  } else { 
    out <- data.frame(cov = covs,
                      predict_over(seq(min(X), max(X),  0.1),alpha, beta))
  }
  out
}

my_quantile <- function(x){
  q <- quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)
  names(q) <- c("inf", "med", "sup")
  q
}

truncate_to_quantile <- function(x, q){
  ifelse(x > quantile(x, q), quantile(x, q), x)
}
