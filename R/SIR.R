# Sampling-importance resampling ------------------------------------------

# Sampling-importance resampling approximation to the posterior distribution
sir = function(lp, priors, parmodes, algorithm = "Nelder", nStart = 20,
               n = 1e3, nratio = 10, df = 3, gr = NULL, ...) {
  # Retrieve MAP and approximation to variance-covariance matrix
  map_vcov = getMAP(lp, priors, parmodes, algorithm = algorithm, gr = gr,
                    nStart = nStart, ...)
  # Generate large sample using mvt
  nlarge = n*nratio
  theta = mvnfast::rmvt(n = nlarge, mu = map_vcov[[1]], sigma = map_vcov[[2]], df = df,
                        ncores = future::availableCores())
  colnames(theta) = names(pars)
  # Calcualte log density for proposal distribution (lp) and target/posterior distribution (lf)
  lprop = mvnfast::dmvt(theta, mu = map_vcov[[1]], sigma = map_vcov[[2]], df = df,
                        log = TRUE, ncores  = future::availableCores())
  ltarget = furrr::future_map_dbl(1:nlarge, ~lp(theta[.x,]))
  # Resample without replacement a smaller sample proportionally to ratio of densities
  wt = exp(ltarget - lprop)
  indices = sample(1:nlarge, size = n, prob = wt, replace = FALSE)
  theta[indices,]
}

