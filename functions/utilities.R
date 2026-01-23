logit <- function(p) log(p/(1-p))
invlogit <- function(x) 1/(1 + exp(-x))

getSparse <- function(x, thr){
  # NB: I can replace 0 in the following lines by Inf and (> by !=). It should not change anything with the 
  # present formulation of the model, but if I want to model omega as a function of neighbors the previous 
  # year (e.g. recue effect) I might need to do it.
  
  x[x > thr] <- 0
  d <- as.numeric(x[x > 0])
  tmp <- apply(x, 1, FUN = function(x){which(x > 0)})
  i <- unlist(tmp)
  p <- cumsum(c(0, reduce(map(tmp, length), c)))
  
  list(d = d,
       i = i,
       p = p)
}
