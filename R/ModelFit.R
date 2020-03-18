
# Fit a non-linear model to data using Bayesian approach ------------------
fitModel = function(model, pars, method = "quad", algorithm = "BFGS", nStart = 20, ...) {

  # Prepare parameters and functions
  checkPars(pars, model$parmodes)
  parsu = unconstrain1D(pars, model$parmodes)
  fn = function(x) lp(x, model)
  gr = function(x) gr_lp(x, model)

  # Sample the posterior distribution using different numerical methods
  if(method == "quad") {
    theta = quad(lp = fn, priors = model$priors, parmodes = model$parmodes,
                 gr = gr, algorithm = algorithm, nStart = nStart, ...)
    colnames(theta) = names(pars)
    theta = constrain2D(theta, model$parmodes)
  } else if (method == "sir") {
    theta = sir(lp = fn, priors = model$priors, parmodes = model$parmodes,
                gr = gr, algorithm = algorithm, nStart = nStart, ...)
    colnames(theta) = names(pars)
    theta = constrain2D(theta, model$parmodes)
  } else {
    stop('Method must be one of the following: "quad" or "sir"')
  }

  # Update the model with the sample from the posterior and return to the user
  model$posterior = theta
  return(model)
}

# Wrappers to construct log-posterior and its gradient --------------------

# Calculate log-likelihood...
ll = function(pars, fit) {
  data = fit$data
  Apred = fit$model(c(pars, fit$fixed), data)
  LL = sum(dnorm(x = data$A, mean = Apred, sd = pars[length(pars)], log = TRUE))
}

# ...and its gradient
gr_ll = function(pars, fit) {
  np = length(pars)
  parnames = names(pars)[-np]
  data = fit$data
  A = data$A
  sigma = pars[np]
  Apred = fit$model(c(pars, fit$fixed), data)
  gr_Apred = fit$gradient(c(pars, fit$fixed), data, parnames)
  gr_LL = numeric(np)
  for(i in 1:(np - 1)) {
    gr_LL[i] = sum((A - Apred)*gr_Apred[,i])/sigma^2
  }
  gr_LL[np] = -length(A)/sigma + sum((A - Apred)^2)/sigma^3
  gr_LL
}


# Calculate log posterior with Jacobian corrections...
lp = function(pars, fit) {
  pars = constrain1D(pars, fit$parmodes)
  LP = ll(pars, fit)
  for(i in 1:length(pars)) {
    name = names(pars)[i]
    prior = fit$priors[[name]]
    LP = LP + dnorm(pars[i], prior[1], prior[2], log = TRUE)
    if(fit$parmodes[name] == ">0")
      LP = LP + log(pars[i])
    else if(fit$parmodes[name] == "0-1")
      LP = LP + log(pars[i]) + log(1 - pars[i])
  }
  LP
}

# ...and its gradient
gr_lp = function(upars, fit) {
  pars = constrain1D(upars, fit$parmodes)
  np = length(pars)
  gr_lp = gr_ll(pars, fit)
  for(i in 1:np) {
    name = names(pars)[i]
    prior = fit$priors[[name]]
    gr_lp[i] = gr_lp[i] - (pars[i] - prior[1])/prior[2]^2
    if(fit$parmodes[name] == ">0") {
      gr_lp[i] = (gr_lp[i] + 1/pars[i])*exp(upars[i])
    } else if(fit$parmodes[name] == "0-1") {
      gr_lp[i] = (gr_lp[i] + 1/pars[i] - 1/(1 - pars[i]))*exp(-upars[i])/(1 + exp(-upars[i]))^2
    }
  }
  gr_lp
}


