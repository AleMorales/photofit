

#' Summary of a fitted model
#'
#' @param fit Object of class \code{photofit}
#' @param cor Method to calculate correlation (see \code{\link{cor}}). If
#'            \code{NULL}, correlations will not be reported.
#' @param include_error Whether to include the experimental error as parameter
#'
#' @details Compute the mean, median and typical quantiles from each posterior
#' marginal distribution as well as the correlation among these parameters.
#'
#' The methods for calculated the correlations are \code{"pearson"},
#' \code{"spearman"} and \code{"kendall"}.
#'
#' @return Summary of posterior distribution
#' @export
summary.photofit = function(fit, cor = NULL, include_error = TRUE) {
  if(include_error)
    sample = fit$posterior
  else
    sample = fit$posterior[,-which(colnames(fit$posterior) == "sigma")]
  mu = colMeans(sample)
  sd = apply(sample, 2, sd)
  med = apply(sample, 2, median)
  mad = apply(sample, 2, mad)
  lb = apply(sample, 2, quantile, probs = 0.055)
  ub = apply(sample, 2, quantile, probs = 0.945)
  marginal = cbind(mu, sd, med, mad, lb, ub)
  rownames(marginal) = colnames(sample)
  colnames(marginal) = c("mean", "sd", "median", "mad", "5.5%", "94.5%")
  cat("Marginal posterior estimates:\n\n")
  print(marginal)
  if(!is.null(cor)) {
    cors = cor(sample, method = cor)
    cat("\nCorrelation matrix:\n\n")
    print(cors)
    out = list(marginal = marginal, correlation = cors)
  } else {
    out = list(marginal = marginal)
  }
  invisible(out)
}

#' Sample from posterior prediction
#'
#' @param fit Object of class \code{photofit}
#' @param newdata New data for which predictions are requested
#' @param nsamples Number of samples to obtain from the posterior prediction
#' @param include_error Indicate whether predictions should include observation error or not.
#'              Defaults to \code{TRUE}.
#'
#' @return Matrix with sample from predictive posterior with \code{nsamples} rows.
#' @export
predict.photofit = function(fit, newdata = fit$data, nsamples = nrow(fit$posterior),
                            include_error = FALSE) {
  # Get sample from posterior
  if(nsamples == nrow(fit$posterior)) {
    sample = fit$posterior %>% unclass %>% as.data.frame %>% as.list
  } else {
    sample = fit$posterior[sample(1:nrow(fit$posterior), nsamples,
                              replace = TRUE),] %>% unclass  %>% as.data.frame %>% as.list
  }
  # Calculate model predictions for each row in newdata and sample
  out = matrix(NA, ncol = nrow(newdata), nrow = length(sample[[1]]))
  for(i in 1:ncol(out)) {
    A = fit$model(c(sample, fit$fixed), newdata[i,])
    if(include_error) {
      A = rnorm(length(A), mean = A, sd = sample[["sigma"]])
    }
    out[,i] = A
  }
  return(out)
}



#' Sample from posterior prediction of limiting factors
#'
#' @param fit Object of class \code{photofit}
#' @param newdata New data for which predictions are requested
#' @param nsamples Number of samples to obtain from the posterior prediction
#' @param include_error Indicate whether predictions should include observation error or not.
#'              Defaults to \code{TRUE}.
#'
#' @return Matrix with sample from predictive posterior with \code{nsamples} rows.
#' @export
predictLimitations = function(fit, newdata = fit$data, nsamples = nrow(fit$posterior),
                            include_error = FALSE) {
  # Get sample from posterior
  if(nsamples == nrow(fit$posterior)) {
    sample = fit$posterior %>% unclass %>% as.data.frame %>% as.list
  } else {
    sample = fit$posterior[sample(1:nrow(fit$posterior), nsamples,
              replace = TRUE),] %>% unclass  %>% as.data.frame %>% as.list
  }
  # Calculate model predictions for each row in newdata and sample
  out = array(NA, dim = c(nrow = length(sample[[1]]), nrow(newdata), 3))
  for(i in 1:ncol(out)) {
    A = fit$limitations(c(sample, fit$fixed), newdata[i,])
    if(include_error) {
      dims = dim(A)
      A = rnorm(length(A), mean = A, sd = sample[["sigma"]])
      dim(A) = dims
    }
    for(j in 1:3) {
      out[,i,j] = A[,j]
    }
  }
  return(out)
}



#' Plot results of a model fit
#'
#' @param fit Object that results fof class \code{photofit}
#' @param type String indicated the type of figure to plot. See "details".
#' @param prob Probability use to calculate prediction intervals for the plot of type "response".
#' @param ... Other arguments to the function "predict" for plots of type "response"
#'
#' @details Three types of plots are supported, depending on the string assigned to the argument
#' \code{type}:
#'
#' 1. \code{type = "response"}: This will plot the measured response curve (as points),
#' the median prediction and the quantile intervals for the mean and for future observations
#' (i.e. taking into account measurement error) for the given probability (argument \code{prob}).
#'
#'
#' 2. \code{type = "marginal"}" This will be display density plots for every marginal distribution
#' in the posterior.
#'
#' 3. \code{type = "pairs"}: This will display contour plot of probability density for each
#' pair of parameters, the marginal distributions (as in the previous type) and the correlation
#' between every pair of parameters (using Pearson's correlation coefficient).
#'
#' @export
#' @examples
#' data = data.frame(PAR = c(25, 50, 100, 150, 250, 350, 501, 503, 750, 1000),
#'                   A = c(0.51, 2.04, 4.41, 6.08, 8.03, 9.36, 9.85, 10.51, 11.14, 11.76))
#' fit = fitNRH(data)
#' summary(fit)
#' plot(fit)
#' plot(fit, type = "limitations")
#' plot(fit, type = "marginal")
#' plot(fit, type = "pairs")
plot.photofit = function(fit, type = "response", prob = 0.89, include_error = FALSE,...) {
  # Choose whether to plot the estimated experimental error or not
  if(include_error) {
    sample = fit$posterior
  } else {
    sample = fit$posterior[,-which(colnames(fit$posterior) == "sigma")]
  }
  # Switch to the correct type of plot
  curpar = par(no.readonly = TRUE)
  if(type == "response") {
    plotResponse(fit, prob)
  } else if (type == "limitations") {
    plotLimitations(fit, prob, include_error = include_error)
  } else if (type == "marginal") {
    par(mfrow = c(2,2), mar = c(5,3.5,0.5,0.5), las = 1)
    for(name in colnames(sample)) {
      plot(density(sample[,name], bw = "SJ"), xlab = name, ylab = "", main = NA)
    }
  } else if (type == "pairs") {
    p = ggpairs(as.data.frame(sample),
                lower = list(continuous = "density"),
                upper = list(continuous = wrap("cor",method = "pearson")),
                diag = list(continuous = wrap("densityDiag", bw = "SJ")),
    )
    print(p)
  } else {
    stop('The argument "type" must be "response", "limitations, "marginal"
         or "pairs".')
  }
  par(curpar)
}

