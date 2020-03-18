# Quadratic approximation to posterior distribution -----------------------

# Generate random value from the prior distribution (depends on the parameter
# mode)
rprior = function(par, mode) {
 if(mode == "-") {
   rnorm(1, mean = par[1], sd = par[2]/par[1])
 } else if(mode == "0-1") {
   distros::rbeta2(1, mu = par[1], cv = par[2])
 } else {
   distros::rlnorm2(1, mu = par[1], cv = par[2]/par[1])
 }
}

# Generate random values from the prior distributions and run the optimization
# from that starting point
runOptim = function(fun, priors, parmodes) {
  pars = map2_dbl(priors, parmodes, ~rprior(.x, .y))
  pars = unconstrain1D(pars, parmodes)
  fun(pars)
}

# Retrieve MAP and vcov approximation based on curvature of log posterior
getMAP = function(lp, priors, parmodes, algorithm = "Nelder", itermax = 1e5,
                  nStart = 10, gr = NULL, ...) {

  # Fitting fun
  fitfun = function(x) {
    optim(par = x, fn = lp, gr = gr, ...,
          method = algorithm, hessian = FALSE,
          control = list(fnscale = -1, parscale = abs(x), maxit = itermax))
  }

  # Start from different random starting points (this should run in parallel!)
  fits = map(1:nStart, ~runOptim(fitfun, priors, parmodes))

  # Choose the best fit case
  values = map_dbl(fits, ~.x$value)
  bestfit = fits[[which.max(values)]]

  # Extract MAP
  MAP = bestfit$par
  names(MAP) = names(priors)

  # Compute Hessian as Jacobian from the gradient gunction
  hess = try(numDeriv::jacobian(func = gr, x = MAP, method.args = list(eps=1e-6,
            d=1e-6, zero.tol=sqrt(.Machine$double.eps), r=6, v=2), ...), silent = TRUE)
  if (inherits(hess, "try-error")) stop("Unable to compute Hessian")
  if(!isSymmetric(hess, tol = 1e-7)) hess = 0.5 * (t(hess) + hess)

  # Return the MAP and variance-covariance matrix
  vcov = solve(-hess)
  list(MAP = MAP, vcov = vcov)
}


# Quadratic approximation to the posterior distribution
quad = function(lp, priors, parmodes, algorithm = "Nelder", nStart = 20, n = 1e3, gr = NULL, ...) {
  # Retrieve MAP and approximation to variance-covariance matrix
  map_vcov = getMAP(lp, priors, parmodes, algorithm = algorithm, gr = gr,
                    nStart = nStart, ...)
  # Generate sample from multivariate normal distribution
  cholSigma = try(chol(map_vcov[[2]]), silent = TRUE)
  if(inherits(cholSigma, "try-error")) {
    cholSigma = suppressWarnings(try(chol(map_vcov[[2]], pivot = TRUE), silent = TRUE))
  }
  sample = try(mvnfast::rmvn(n = n, mu = map_vcov[[1]], sigma = cholSigma,
                ncores = future::availableCores(), isChol = TRUE), silent = TRUE)
  if(inherits(sample, "try-error")) {
    warning("Could not calculate covariance matrix - removing correlations among parameters")
    vcov = abs(map_vcov[[2]])
    vcov[lower.tri(vcov)] = 0
    vcov[upper.tri(vcov)] = 0
    sample = mvnfast::rmvn(n = n, mu = map_vcov[[1]], sigma = vcov,
                           ncores = future::availableCores())
  }
  sample
}
