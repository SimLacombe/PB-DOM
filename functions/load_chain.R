load_chain <- function(path, T, N, z_mean = TRUE){
  out <- map(list.files(path, full.names = TRUE, pattern = "MCMC"), readRDS) %>% 
    reduce(rbind)
  zm <- NULL
  if(z_mean){
    zm <- apply(out[, grepl("z\\[", colnames(out))], 2, mean) %>%
      matrix(nrow = T, ncol = N)
    out <- out[, !grepl("z\\[", colnames(out))]
  }
  
  list(pars = out, zm = zm)
}
