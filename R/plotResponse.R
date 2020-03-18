
# Generic method that varies depending on the class of object

plotResponse = function (fit, ...) {
  UseMethod("plotResponse", fit)
}

# Generate beautiful labels for the axes of the figures
prettify = function(x) {
  if(x == "PAR") {
    expression(PAR~(mu*mol~m^{-2}~s^{-1}))
  } else if (x == "Ci") {
    expression(Ci~(mu*mol~mol^{-1}))
  } else {
    x
  }
}


# Default response plot (one response curve)
plotResponse.photofit = function(fit, prob) {
  par(mar = c(5,5,0.5,0.5))
  # Construct continuous range of values within the range of the predictor
  rangeX = range(fit$data[[fit$predictor]])
  x = seq(rangeX[1], rangeX[2], length.out = 20)
  newdata = data.frame(x, A = NA)
  names(newdata) = c(fit$predictor,"A")
  # Make predictions for the range of X with and without propagation experimental error
  preMean = predict(fit, newdata = newdata, include_error = FALSE)
  preObs = predict(fit, newdata = newdata, include_error = TRUE)
  # Compute quantiles for a given probability from each prediction
  preMedian = apply(preMean, 2, median)
  preLBMean = apply(preMean, 2, quantile, probs = (1 - prob)/2)
  preUBMean = apply(preMean, 2, quantile, probs = 1 - (1 - prob)/2)
  preLBObs = apply(preObs, 2, quantile, probs = (1 - prob)/2)
  preUBObs = apply(preObs, 2, quantile, probs = 1 - (1 - prob)/2)
  # Plot the results
  orig_x = fit$data[[fit$predictor]]
  plot(orig_x, fit$data$A, pch = 16, las = 1,
       ylab = expression(A~(mu*mol~m^{-2}~s^{-1})),
       xlab = prettify(fit$predictor))
  lines(x, preMedian, lwd = 1.5)
  polygon(c(x, rev(x)), c(preUBMean, rev(preLBMean)), col = rgb(0,0,0,0.45),
          border = NA)
  polygon(c(x, rev(x)), c(preUBObs, rev(preLBObs)),
          col = rgb(0,0,0,0.15), border = NA)
}


# Response plot when both LRC and ACi are fitted simultaneously
#' @export
plotResponse.ACiLRC = function(fit, prob) {
  par(mfrow = c(1,2), mar = c(5,5,0.5,0.5))
  for(i in 1:2) {
    curveType = fit$predictor[i]
    # Construct continuous range of values within the range of the predictor
    curve_data = fit$data[which(fit$curves == curveType),]
    orig_x = curve_data[,curveType]
    rangeX = range(orig_x)
    x = seq(rangeX[1], rangeX[2], length.out = 50)
    # THIS PART DOES NOT GENERALIZE TO THE CONCEPT OF "N CURVES"
    if(curveType == "PAR")
      newdata = data.frame(PAR = x, A = NA, Ci = mean(curve_data$Ci))
    else
      newdata = data.frame(Ci = x, A = NA, PAR = mean(curve_data$PAR))
    # Make predictions for the range of X with and without propagation experimental error
    preMean = predict(fit, newdata = newdata, include_error = FALSE)
    preObs = predict(fit, newdata = newdata, include_error = TRUE)
    # Compute quantiles for a given probability from each prediction
    preMedian = apply(preMean, 2, median)
    preLBMean = apply(preMean, 2, quantile, probs = (1 - prob)/2)
    preUBMean = apply(preMean, 2, quantile, probs = 1 - (1 - prob)/2)
    preLBObs = apply(preObs, 2, quantile, probs = (1 - prob)/2)
    preUBObs = apply(preObs, 2, quantile, probs = 1 - (1 - prob)/2)
    # Plot the results
    plot(orig_x, curve_data$A, pch = 16, las = 1,
         ylab = expression(A~(mu*mol~m^{-2}~s^{-1})),
         xlab = prettify(curveType))
    lines(x, preMedian, lwd = 1.5)
    polygon(c(x, rev(x)), c(preUBMean, rev(preLBMean)), col = rgb(0,0,0,0.45),
            border = NA)
    polygon(c(x, rev(x)), c(preUBObs, rev(preLBObs)),
            col = rgb(0,0,0,0.15), border = NA)
  }
}
