% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FvCBAci.R
\name{fitACi}
\alias{fitACi}
\title{Fit A-Ci curve}
\usage{
fitACi(
  data,
  priors = priors_ACi(),
  pars = initial_ACi(),
  fixed = c("KmC", "KmO"),
  method = "quad",
  algorithm = "BFGS",
  ...
)
}
\arguments{
\item{data}{Measurements as a \code{data.frame} with two columns called \code{Ci} and \code{A}
representing intercellular CO₂ molar fraction and net CO₂ assimilation, respectively.}

\item{priors}{List of parameters for the prior distributions as generated by
\code{\link{priors_ACi}}.}

\item{pars}{Vector of initial values of parameters as generated by
\code{\link{initial_ACi}}.}

\item{fixed}{Vector with names of parameters that will not be fitted to the
data but fixed to their initial values. Same names as in \code{\link{priors_ACi}}
and \code{\link{initial_ACi}} should be used.}

\item{method}{Method to approximate the posterior distribution of parameters. See Details.}

\item{algorithm}{Name of non-linear optimization algorithm used to calculate
maximum a posteriori. See Details.}
}
\description{
\code{fitACi} will fit a steady-state model of photosynthesis to a response curve
measuring CO₂ assimilation and intercellular CO₂ molar fraction.
}
\details{
A classic FvCB model that implicitly assumes infinite mesophyll conductance is used.
Net CO₂ assimilation for every measurement point is calculated as the minimum of three
potential rates limited by Rubisco kinetics (\code{Ac}), electron transport (\code{Aj})
and triose phosphate utilisation (\code{At}):

\code{Ac = Vcmax*(Ci - CiStar)/(Ci + KmC*(1 + 210/KmO))}

\code{Aj = J*(Ci - CiStar)/(4*Ci + 8*CiStar)}

\code{At = 3*TPU}

The parameters that can be fitted are:

\itemize{
 \item \code{Vcmax} Maximum rate of carboxylation (μmol m⁻² s⁻¹).
 \item \code{CiStar} Intercellular CO₂ compensation point (μmol mol⁻¹).
 \item \code{Kmc} Michaelis-Menten constant of carboxylation with respect to CO₂ (μmol mol⁻¹).
 \item \code{Kmo} Michaelis-Menten constant of carboxylation with respect to O₂ (mmol mol⁻¹).
 \item \code{J} Maximum rate of electron transport under measurement conditions (μmol m⁻² s⁻¹).
 \item \code{TPU} Maximu rate of triose phosphate utilisation (μmol m⁻² s⁻¹).
 \item \code{Rd} Mitochondrial respiration in the light (μmol m⁻² s⁻¹).
 \item \code{sigma} Standard deviation of measurement error (μmol m⁻² s⁻¹).
}

The default initial values are generated by \code{\link{initial_ACi}} whereas
the mean and standard deviations of Gaussian prior distributions are generated
by \code{\link{priors_ACi}}.

Values allowed for the \code{method} argument are \code{"quad"} (quadratic
approximation) and \code{"sir"} (sampling importance resampling).

All methods to approximate the posterior start with a calculated of the
maximum a posteriori. The argument \code{algorithm} corresponds to the \code{method}
argument of the function \code{\link{optimx}} (see help for this function for
a list of possible algorithms).
}
\examples{
data = data.frame(Ci = c(245.39, 164.70, 87.59, 50.22, 332.55, 333.15, 421.98,
                  515.53, 612.01, 707.88, 904.34, 1099.13, 1295.58, 1485.47),
           A = c(8.90, 5.72, 1.94, -0.18, 11.38, 11.49, 13.49, 14.66, 15.10,
                 15.61, 15.49, 15.43, 14.97, 15.51))

# Quadratic approximation
fitQuad = fitACi(data = data, fixed = c("KmC", "KmO"))
summary(fitQuad)
plot(fitQuad)
plot(fitQuad, type = "marginal")
plot(fitQuad, type = "pairs")

# Quadratic approximation - Avoiding TPU
fitQuad2 = fitACi(data = data, pars = initial_ACi(TPU = 100),
                 fixed = c("KmC", "KmO", "TPU"))
summary(fitQuad2)
plot(fitQuad2)
plot(fitQuad2, type = "pairs")

# SIR approximation
fitSIR = fitACi(data = data, pars = initial_ACi(TPU = 100),
                 fixed = c("KmC", "KmO", "TPU"))
summary(fitSIR)
plot(fitSIR)
plot(fitSIR, type = "marginal")
plot(fitSIR, type = "pairs")
}
\seealso{
initial_ACi, priors_ACi
}
