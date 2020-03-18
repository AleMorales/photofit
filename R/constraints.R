
# Parameters of a prior distribution (Normal in unconstrained space)
#' @export
prior = function(mu = 0, sd = 1) {
  c(mu = mu, sd = sd)
}

# Functions to help with constraints on parameter space (move from constrained
# to unconstrained)
logistic = function(x) 1/(exp(-x) + 1)
logit = function(x) log(x/(1 - x))
identity = function(x) x

unconstrain = c("-" = identity, "0-1" = logit, ">0" = log)
constrain = c("-" = identity, "0-1" = logistic, ">0" = exp)

constrain1D = function(pars, parmodes) {
  for(name in names(pars))
    pars[[name]] = constrain[[parmodes[name]]](pars[[name]])
  pars
}

constrain2D = function(pars, parmodes) {
  for(name in colnames(pars)) {
    pars[,name] = constrain[[parmodes[name]]](pars[,name])
  }
  pars
}

unconstrain1D = function(pars, parmodes) {
  for(name in names(pars))
    pars[[name]] = unconstrain[[parmodes[name]]](pars[[name]])
  pars
}


# Need to avoid values that cause problems with parscale (0) or with some
# transformations
checkPars = function(pars, parmodes) {
  if(any(pars == 0)) stop("Initial values cannot be 0")
  if(any(pars == 1)) stop("Initial values cannot be 1")
}


# Prior distributions depending on the type of constraint

