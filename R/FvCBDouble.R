# Simple ACi and LRC model (no gm) ---------------------------------------------

ACiLRC_limitations  = function(pars, data) {
  # Parameters
  Vcmax = pars[["Vcmax"]]
  KmC   = pars[["KmC"]]
  KmO   = pars[["KmO"]]
  CiStar = pars[["CiStar"]]
  alpha = pars[["alpha"]]
  theta = pars[["theta"]]
  Jmax  = pars[["Jmax"]]
  TPU   = pars[["TPU"]]
  Rd    = pars[["Rd"]]
  # Environment
  Ci    = data$Ci
  PAR   = data$PAR
  O2    = 210
  # Limiting factors to CO2 assimilation
  Ac    = Vcmax*(Ci - CiStar)/(Ci + KmC*(1 + O2/KmO))
  J     = (alpha*PAR + Jmax - sqrt((alpha*PAR + Jmax)^2 -
                                     4*alpha*theta*Jmax*PAR))/(2*theta)
  Aj    = J*(Ci - CiStar)/(4*Ci + 8*CiStar)
  At    = 3*TPU
  cbind(Ac, Aj, At) - Rd
}

ACiLRC = function(pars, data) {
  lims = ACiLRC_limitations(pars, data)
  pmin(pmin(lims[,1], lims[,2]), lims[,3])
}

gr_ACiLRC = function(pars, data, parnames) {
  # Parameters
  Vcmax = pars[["Vcmax"]]
  KmC   = pars[["KmC"]]
  KmO   = pars[["KmO"]]
  CiStar = pars[["CiStar"]]
  alpha = pars[["alpha"]]
  theta = pars[["theta"]]
  Jmax  = pars[["Jmax"]]
  TPU   = pars[["TPU"]]
  Rd    = pars[["Rd"]]
  # Environment
  Ci    = data$Ci
  PAR   = data$PAR
  O2 = 210
  n = nrow(data)
  # Limiting factors to CO2 assimilation
  Ac = Vcmax*(Ci - CiStar)/(Ci + KmC*(1 + 210/KmO))
  J     = (alpha*PAR + Jmax - sqrt((alpha*PAR + Jmax)^2 -
                                     4*alpha*theta*Jmax*PAR))/(2*theta)
  Aj    = J*(Ci - CiStar)/(4*Ci + 8*CiStar)
  At = rep(3*TPU, n)

  # Figure out which are the limiting steps
  limiting = purrr:::map_dbl(1:n, ~which.min(c(Ac[.x], Aj[.x], At[.x])))

  # Derivatives for Ac-limited photosynthesis
  .e1 = 1 + 210/KmO
  .e2 = Ci + KmC * .e1
  .e3 = Ci - CiStar
  .e4 = .e2^2
  dA_dVcmax = ifelse(limiting == 1, .e3/.e2, 0)
  dA_dKmC   = ifelse(limiting == 1, -(Vcmax*.e1 * .e3/.e4), 0)
  dA_dKmO   = ifelse(limiting == 1, 210*(KmC * Vcmax * .e3/(KmO^2*.e4)), 0)
  dA_dCiStar_Vcmax = -(Vcmax/.e2)

  # Derivatives for Aj-limited photosynthesis
  .e1 <- alpha * PAR
  .e2 <- .e1 + Jmax
  .e4 <- alpha * Jmax * PAR
  .e8 <- 4 * Ci + 8 * CiStar
  .e9 <- sqrt(.e2^2 - 4 * (.e4 * theta))
  .e10 <- Ci - CiStar
  .e11 <- theta * .e8
  .e12 <- 2 * .e2
  .e13 <- 2 * .e11
  .e14 <- .e2 - .e9
  dA_dJmax = ifelse(limiting == 2, (1 - 0.5*((.e12 - 4*(.e1*theta))/.e9))*.e10/.e13, 0)
  dA_dalpha = ifelse(limiting == 2, PAR*(1 - 0.5*((.e12 - 4*(Jmax * theta))/.e9))*.e10/.e13, 0)
  dA_dtheta = ifelse(limiting == 2, (.e4/(theta*.e9) - 2*(.e14/(2*theta)^2))*.e10/.e8, 0)
  dA_dCiStar_J = -((0.5 + 4*(.e10/.e8))*.e14/.e11)

  # Combine and others
  dA_dCiStar = ifelse(limiting == 1, dA_dCiStar_Vcmax,
                     ifelse(limiting == 2, dA_dCiStar_J, 0))
  dA_dRd    = rep(-1, n)
  dA_dTPU   = ifelse(limiting == 3, 3, 0)

  cbind(Vcmax = dA_dVcmax, KmC = dA_dKmC, KmO = dA_dKmO, CiStar = dA_dCiStar,
        Jmax = dA_dJmax, alpha = dA_dalpha, theta = dA_dtheta, TPU = dA_dTPU,
        Rd = dA_dRd)[,parnames]
}


# Generate list of parameters for prior distributions

#' @export
priors_ACiLRC = function(Vcmax = prior(100, 100),
                      CiStar = prior(15, 15),
                      KmC = prior(230, 50),
                      KmO = prior(195, 40),
                      Jmax = prior(150,150),
                      alpha = prior(0.3,0.1),
                      theta = prior(0.7,1),
                      TPU = prior(10, 10),
                      Rd = prior(1,1),
                      sigma = prior(1,1)) {
  list(Vcmax = Vcmax, CiStar = CiStar, KmC = KmC, KmO = KmO, Jmax = Jmax,
       alpha = alpha, theta = theta, TPU = TPU, Rd = Rd, sigma = sigma)
}


#' @export
initial_ACiLRC = function(Vcmax = 100, CiStar = 15, KmC = 230, KmO = 195, Jmax = 150,
                          alpha = 0.30, theta = 0.7, TPU = 10, Rd = 0.75, sigma = 0.9) {
  c(Vcmax = Vcmax, CiStar = CiStar, KmC = KmC, KmO = KmO, Jmax = Jmax,
    alpha = alpha, theta = theta, TPU = TPU, Rd = Rd, sigma = sigma)
}


# Constraint on each parameter
parmodesACiLRC = c(Vcmax = ">0", CiStar = ">0", KmC = ">0", KmO = ">0", Jmax = ">0",
                   alpha = "0-1", theta = "0-1", TPU = ">0", Rd = ">0", sigma = ">0")


#' @export
fitACiLRC = function(data, priors = priors_ACiLRC(), pars = initial_ACiLRC(),
                     fixed = c(), curves = NULL, method = "quad",
                     algorithm = "BFGS", ...) {
  if(is.null(curves)) stop("Please use the argument 'curves' to define to which
                         curve each data point belongs to.")
  # Fixed parameters not to be fitted to the data
  fixed_pars = list()
  for(name in fixed) {
    fixed_pars[[name]] = pars[[name]]
    pars = pars[-which(names(pars) == name)]
  }
  model = structure(list(data = as.data.frame(data), priors = priors, model  = ACiLRC,
                         limitations = ACiLRC_limitations, gradient = gr_ACiLRC,
                         parmodes = parmodesACiLRC,
                         fixed = fixed_pars, curves = curves,
                         predictors = c("Ci","PAR")), class = "ACiLRC")
  fit = fitModel(model, pars, method = method, algorithm = algorithm, ...)
  class(fit) = c("ACiLRC", "photofit")
  return(fit)
}

