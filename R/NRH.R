
# NRH Model -------------------------------------------------------------------

NRH = function(pars, data) {
  alpha = pars[["alpha"]]
  theta = pars[["theta"]]
  Amax = pars[["Amax"]]
  Rd = pars[["Rd"]]
  PAR = data$PAR
  (alpha*PAR + Amax - sqrt((alpha*PAR + Amax)^2 - 4*alpha*theta*Amax*PAR))/(2*theta) - Rd
}

gr_NRH = function(pars, data, parnames) {
  alpha = pars[["alpha"]]
  theta = pars[["theta"]]
  Amax = pars[["Amax"]]
  Rd = pars[["Rd"]]
  PAR = data$PAR
  dA_dalpha = (1/2)*(PAR + (-1/2)*(2*PAR*(Amax + PAR*alpha) - 4*PAR*Amax*theta)/sqrt(-4*PAR*Amax*theta*alpha + (Amax + PAR*alpha)^2))/theta
  dA_dtheta = (-1/2)*(Amax + PAR*alpha - sqrt(-4*PAR*Amax*theta*alpha + (Amax + PAR*alpha)^2))/theta^2 + PAR*Amax*alpha/(theta*sqrt(-4*PAR*Amax*theta*alpha + (Amax + PAR*alpha)^2))
  dA_dAmax = (1/2)*(1 + (-1/2)*(-4*PAR*theta*alpha + 2*(Amax + PAR*alpha))/sqrt(-4*PAR*Amax*theta*alpha + (Amax + PAR*alpha)^2))/theta
  dA_dRd = -1
  cbind(alpha = dA_dalpha, theta = dA_dtheta, Amax = dA_dAmax, Rd = dA_dRd)[,parnames]
}

# Generate list of parameters for prior distributions

#' @export
priors_NRH = function(alpha = prior(0.05, 0.02),
                      theta = prior(0.7, 1),
                      Amax  = prior(15, 15),
                      Rd   = prior(1, 1),
                      sigma = prior(1, 1)) {
  list(alpha = alpha, theta = theta, Amax = Amax, Rd = Rd, sigma = sigma)
}


#' @export
initial_NRH = function(alpha = 0.05, theta = 0.7, Amax = 15,
                       Rd = 0.75, sigma = 0.5) {
  c(alpha = alpha, theta = theta, Amax = Amax, Rd = Rd, sigma = sigma)
}

# Constraint on each parameter
parmodesNRH = c(alpha = "0-1", theta = "0-1", Amax = ">0", Rd = ">0", sigma = ">0")


#' @export
fitNRH = function(data, priors = priors_NRH(), pars = initial_NRH(),
                  fixed = c(), method = "quad", algorithm = "BFGS", ...) {
  # Fixed parameters not to be fitted to the data
  fixed_pars = list()
  for(name in fixed) {
    fixed_pars[[name]] = pars[[name]]
    pars = pars[-which(names(pars) == name)]
  }
  modelNRH = structure(list(data = as.data.frame(data), priors = priors, model  = NRH,
                            gradient = gr_NRH, parmodes = parmodesNRH,
                            fixed = fixed_pars, predictor = "PAR"), class = "NRH")
  fit = fitModel(modelNRH, pars, method = method, algorithm = algorithm, ...)
  class(fit) = c("NRH", "photofit")
  return(fit)
}


