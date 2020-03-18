
# Generic method that varies depending on the class of object
plotLimitations = function (fit, ...) {
  UseMethod("plotLimitations", fit)
}

# Generate beautiful labels for the x axes of the response curves
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
plotLimitations.photofit = function(fit, prob, include_error = FALSE) {
  par(mar = c(5,5,0.5,0.5))
  # Construct continuous range of values within the range of the predictor
  rangeX = range(fit$data[[fit$predictor]])
  x = seq(rangeX[1], rangeX[2], length.out = 50)
  newdata = data.frame(x, A = NA)
  names(newdata) = c(fit$predictor,"A")
  # Make predictions for the range of X with and without propagation experimental error
  preMean = predictLimitations(fit, newdata = newdata, include_error = FALSE)
  if(include_error)
    preObs  = predictLimitations(fit, newdata = newdata, include_error = TRUE)
  # Compute quantiles for a given probability from each prediction
  preMedian = apply(preMean, 3, function(x) apply(x, 2, median))
  preLBMean = apply(preMean, 3, function(x) apply(x, 2, quantile, probs = (1 - prob)/2))
  preUBMean = apply(preMean, 3, function(x) apply(x, 2, quantile, probs = 1 - (1 - prob)/2))
  if(include_error) {
    preLBObs  = apply(preObs, 3, function(x) apply(x, 2, quantile, probs = (1 - prob)/2))
    preUBObs  = apply(preObs, 3, function(x) apply(x, 2, quantile, probs = 1 - (1 - prob)/2))
  }
  # Plot the results
  orig_x = fit$data[[fit$predictor]]
  plot(orig_x, fit$data$A, pch = 16, las = 1,
       ylab = expression(A~(mu*mol~m^{-2}~s^{-1})),
       xlab = prettify(fit$predictor))
  for(i in 1:3) {
    lines(x, preMedian[,i], lwd = 2, col = i + 1)
    polygon(c(x, rev(x)), c(preUBMean[,i], rev(preLBMean[,i])),
            col = rgb(c(1,0,0)[i],c(0,1,0)[i],c(0,0,1)[i],0.45),
            border = NA)
    if(include_error)
        polygon(c(x, rev(x)), c(preUBObs[,i], rev(preLBObs[,i])),
              col = rgb(c(1,0,0)[i],c(0,1,0)[i],c(0,0,1)[i],0.2), border = NA)
  }
  legend("bottomright", legend = c("Ac", "Aj", "At"), col = 2:4, fill = 2:4)
}


# Response plot when both LRC and ACi are fitted simultaneously
#' @export
plotLimitations.ACiLRC = function(fit, prob, include_error = FALSE) {
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
    preMean = predictLimitations(fit, newdata = newdata, include_error = FALSE)
    if(include_error)
      preObs  = predictLimitations(fit, newdata = newdata, include_error = TRUE)
    # Compute quantiles for a given probability from each prediction
    preMedian = apply(preMean, 3, function(x) apply(x, 2, median))
    preLBMean = apply(preMean, 3, function(x) apply(x, 2, quantile, probs = (1 - prob)/2))
    preUBMean = apply(preMean, 3, function(x) apply(x, 2, quantile, probs = 1 - (1 - prob)/2))
    if(include_error) {
      preLBObs  = apply(preObs, 3, function(x) apply(x, 2, quantile, probs = (1 - prob)/2))
      preUBObs  = apply(preObs, 3, function(x) apply(x, 2, quantile, probs = 1 - (1 - prob)/2))
    }
    # Plot the results
    plot(orig_x, curve_data$A, pch = 16, las = 1,
         ylab = expression(A~(mu*mol~m^{-2}~s^{-1})),
         xlab = prettify(curveType))
    for(j in 1:3) {
      lines(x, preMedian[,j], lwd = 2, col = j + 1)
      polygon(c(x, rev(x)), c(preUBMean[,j], rev(preLBMean[,j])),
              col = rgb(c(1,0,0)[j],c(0,1,0)[j],c(0,0,1)[j],0.45),
              border = NA)
      if(include_error)
        polygon(c(x, rev(x)), c(preUBObs[,j], rev(preLBObs[,j])),
                col = rgb(c(1,0,0)[j],c(0,1,0)[j],c(0,0,1)[j],0.2), border = NA)
    }
    legend("bottomright", legend = c("Ac", "Aj", "At"), col = 2:4, fill = 2:4)
  }
}
