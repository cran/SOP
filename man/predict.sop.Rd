\name{predict.sop}
\alias{predict.sop}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Prediction from a fitted SOP model}

\description{The function takes a fitted \code{sop} object and produces predictions for the original data if the argument \code{newdata} is not set or predictions for new data if \code{newdata} is specified. Predictions can be accompanied by standard errors, based on the Bayesian posterior distribution of the model coefficients.
}

\usage{
\method{predict}{sop}(object, newdata, type = c("response", "link", "terms"), 
      se.fit = FALSE, ...)
}
\arguments{ 
 \item{object}{a fitted \code{sop} object as produced by \code{sop()}.}
 \item{newdata}{a data frame containing the values of the model covariates at which predictions are required. If this is not provided then predictions corresponding to the original data are returned. If the data frame \code{newdata} is provided then it should contain all the variables needed for prediction: a warning is generated if not. If \code{newdata} contains a variable \code{offset}, it is included into the predictions when \code{type = "link"} and \code{type = "response"}.}
\item{type}{When this has the value \code{"link"} the linear predictor fitted values or predictions (possibly with associated standard errors) 
are returned. 
When \code{type = "terms"} each component of the linear predictor is returned separately (possibly with approximate standard errors): this includes parametric model components, followed by each smooth component, but excludes any offset and any intercept.
When \code{type = "response"} (default) fitted values or predictions on the scale of the response are returned (possibly with approximate standard errors).
}
\item{se.fit}{when this is TRUE (not default) standard error estimates are returned for each prediction.}
\item{...}{ other arguments. Not yet implemented.}

}

\value{ 
A vector/matrix (or list, with elements \code{fit} and \code{se.fit}, is \code{se  = TRUE}) equal to:

\item{"link"}{a vector of linear predictor values.}

\item{"response"}{a vector of linear predictor values on the scale of the response.}

\item{"terms"}{a matrix with a column per term, and may have an attribute "constant".}
}

%\references{
%Rodriguez-Alvarez, M.X., Lee, D. J., Kneib, T., Durban, M., and Eilers, P. (2015). Fast smoothing parameter separation in multidimensional generalized P-splines: the SAP algorithm. \emph{Statistics and Computing}, \bold{25 (5)}, 941--957.

%Rodriguez-Alvarez, M.X., Durban, M., Lee, D. J. and Eilers, P. (2019). On the estimation of variance parameters in non-standard generalised linear mixed models: application to penalised smoothing. \emph{Statistics and Computing}, \bold{29 (3)}, 483--500.

%}
%\author{ Maria Xose Rodriguez-Alvarez  \email{mxrodriguez@bcamath.org},  Manuel Oviedo de la Fuente.}
\seealso{  \code{\link{sop}},  \code{\link{plot.sop}}}

\examples{
library(SOP)
## Example training/set
# Simulate the data
set.seed(123)
n <- 1000
sigma <- 0.5
x <- runif(n)
f0 <- function(x)2*sin(pi*x)
f <- f0(x)
y <- f + rnorm(n, 0, sigma)
da <- data.frame(x = x, y = y)# all data
rand <-  sample(2, 610, replace=TRUE, prob=c(0.6,0.4))
traindata <- da[rand==1,] # training data
valdata <- da[rand==2,] # validation data  
plot(y ~ x, data = traindata, pch = 20, col = gray(.7))
points(y ~ x, data = valdata, pch = 20, col = gray(.2))

# Fit the model in the training data
m0 <- sop(formula = y ~ f(x, nseg = 10), data = traindata)
lines(fitted(m0)[order(traindata$x)]~traindata$x[order(traindata$x)], 
       col="red", lwd=2)

# Predict and plot in the data used for the fit
po <- predict(m0)
plot(y ~ x, data = traindata, pch = 20, col = gray(.7))
lines(po[order(traindata$x)] ~ traindata$x[order(traindata$x)], 
      col="red", lwd=2)

# Predict and plot in new data
pn <- predict(m0, newdata = valdata)
plot(y ~ x, data = traindata, pch = 20, col = gray(.7))
lines(pn[order(valdata$x)] ~ valdata$x[order(valdata$x)], col = "yellow", lwd = 2)

# Example Gamma distribution
# Simulate the data
set.seed(123)
n <- 1000
alpha <- 0.75
x0 <- runif(n)
x1 <- x0*alpha + (1-alpha)*runif(n)
x2 <- runif(n)
x3 <- x2*alpha + (1-alpha)*runif(n)
x4 <- runif(n)
x5 <- runif(n)

f0 <- function(x)2*sin(pi*x)
f1 <- function(x)exp(2*x)
f2 <- function(x) 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10

f <- f0(x0) + f1(x1) + f2(x2)
y <- rgamma(f,exp(f/4),scale=1.2)

df <- data.frame(y = y, x0 = x0, x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)

# Fit the model
m1 <- sop(formula = y ~ f(x0, nseg = 17) + f(x1, nseg = 17) + 
      f(x2, nseg = 17) + f(x3, nseg = 17) + 
      f(x4, nseg = 17) + f(x5, nseg = 17), 
      family = Gamma(link = log), data = df)
summary(m1)

# Predict in a new dataframe
x <- seq(max(c(min(x1),min(x3))), min(c(max(x1),max(x3))), l = 100)
df.p <- data.frame(x0 = x, x1 = x, x2 = x, x3 = x, x4 = x, x5 = x)
p <- predict(m1, type = "terms", newdata = df.p)
colnames(p)

# Plot the different smooth terms
op <- par(mfrow = c(2,3))
plot(m1, select = 1)
lines(x, p[,1], col = "red")
plot(m1, select = 2)
lines(x, p[,2], col = "red")
plot(m1, select = 3)
lines(x, p[,3], col = "red")
plot(m1, select = 4)
lines(x, p[,4], col = "red")
plot(m1, select = 5)
lines(x, p[,5], col = "red")
plot(m1, select = 6)
lines(x, p[,6], col = "red")
par(op)
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..

