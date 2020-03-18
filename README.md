# photofit

**This package is development and the API is expected to change**

## General usage and assumptions

The package `photofit` implements functions to fit specific models of photosynthesis
to specific types of measurement, providing a separate function for each 
combination of model and type of measurements. The measurements consist of one
or more steady-state response curves of net CO₂ assimilation to environmental variables.

For each model, any subset of the parameters (including all parameters) can be
estimated from the data using a Bayesian approach. A Normal prior distribution
parameterized by mean and standard deviation is assumed for each parameter. The 
measurement errors are also assumed to follow a Normal distribution and to be
independent for each measurement points. Auxilliary functions are available with
default values for the prior distributions as well as initial values required by
the optimization algorithms.

The result of fitting the model to the data is a posterior distribution 
representing what we have learnt about the parameters after incorporating the 
data (this is a function of both the data and the prior knowledge incorported into 
the prior distribution). As the posterior distribution does not have a closed
form, the fitting function will generate a sample from this distribution using
a diversity of samplers. Currently, these samplers are:

- Quadratic approximation (`quad`): The maximum of the posterior distribution is estimated
with non-linear optimization in unconstrained scale (i.e using log and logit
transformations on constrained parameters) and a multivariate Normal 
distribution is fitted to the unconstrained posterior by inverting the 
Hessian matrix evaluated at the maximum. From this fitted distribution a random 
sample of values is generated in the original scale. The underlying assumption
is that the posterior distribution can be represented by a copula with Normal,
LogNormal and LogitNormal marginal distributions (depending on constraints).

- Sampling importance resampling (`sir`): A first quadratic approximation to the
posterior is obtained. This approximation is used as proposal distribution and sampling
importance resampling is applied to obtain a more accurate sample from the posterior
distribution without making distributional assumptions.

Several functions are provided to analyse the output of fitting a model:

- `summary` - Calculate summaries of the posterior distribution. Check `?summary.photofit` for details.

- `predict` - Make predictions with the model for new environmental factors. Check `?predict.photofit` for details.

- `predictLimitations` - Like `predict` but it returns the different potentially limiting factors. Check `?predictLimitations` for details.

- `plot` - Visualize the results of fitting the model to the data. The argument `type` determines the type of plot to be created (check `?plot.photofit` for details):  
  - `type = "response"` - The default value, it will plot the original measured curve with the model predictions overlaid over the data.  
  - `type = "limitations"` - Like for `"response"` but each potentially limiting factor is plotted.  
  - `type = "marginal"` - Plot the posterior marginal distribution of each parameter using density plots.  
  - `type = "pairs"` - Matrix of two dimensional density plots for every pair of parameters, marginal density plot for every parameter and correlation coefficients for every pair of parameters.  

## Installation

This package is still not available on CRAN, so it needs to be installed using the 
`devtools` package as follows:

```r
devtools::install_github("AleMorales/photofit")
```

To update to the latest version use:

```r
devtools::update_packages()
```

## Models and response curves currently implemented

- `fitACi` - Fit FvCB model to A-Ci response curve.

- `fitAciLRC` - Fit FvCB model simultaneously to A-Ci and A-PAR curves.

- `fitNRH` - Fit non-rectangular hyperbola to A-PAR curve.


## Example: A-Ci curve

Let's imagine an experiment that generates the following response curve between 
CO~2~ assimilation (A) and intercellular CO~2~ concentration (Ci).

```r
data = data.frame(
  Ci = c(245, 165, 85, 50, 335, 335, 420, 515, 610, 700, 900, 1000, 1300, 1500),
  A = c(8.9, 5.7, 1.9, -0.2, 11.4, 11.5, 13.5, 14.7, 15.1, 15.6, 15.5,  15.4, 15.5, 15.5)
)
```

We want to fit a classic biochemical model of photosynthesis to this curve. This
is achieved with the function `fitAci`. The function contains sensible defaults,
so we can have a default fit by just passing the data:

```r
fit = fitACi(data = data)
```

We can visually check how well the model fits the data with a `plot`:

```r
plot(fit)
```

As a Bayesian approach is used, the prediction is not a single curve, as each 
combination of parameters from the posterior distribution will result in a 
different curve. Therefore, the predictive distribution is represented by a 
curve representing the median of the distribution and two ribbons representing
probability intervals with and without including measurement error (light and
dark ribbons, respectively). We can turn off the ribbons by setting `prob = 0` 
(but remember that this is just the median of a distribution):

```r
plot(fit, prob = 0)
```

We can visualize the different potentially limiting factors by using 
`type = "limitations"`

```r
plot(fit, prob = 0, type = "limitations")
```

Note that TPU does not limit the curve at any point, so there is actually no 
information in the data to update the prior for this parameter. This is quite 
common in 


If we want to get information on the parameters we can use the `summary` function. 
This reports different summary statistics of the marginal posterior distributions:

```r
summary(fit)
```

We can include the correlation matrix of the posterior (in the original scale) by
passing the type of correlation coefficient to calculate (e.g. `cor = "pearson"`):

```r
summary(fit, cor = "pearson")
```

A visual representation of same information can be obtained by plotting with `type = "pairs"`:

```r
plot(fit, type = "pairs")
```

## Future improvements

- The algorithms from BayesianTools will be added to increase the
list of samplers available, including several Markov Chain Monte Carlo and
Sequential Monte Carlo algorithms.

- Newer models and response curves (including combining multiple response curves) 
will be added/

