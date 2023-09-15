% File SOP/man/sop.control.Rd
% Part of the R package, http://www.R-project.org
% Copyright 1995-2007 R Core Development Team
% Distributed under GPL 2 or later

\name{sop.control}
\alias{sop.control}
\title{Function for controlling SOP fitting}
\description{
The function controls some of the fitting parameters of \code{\link{sop}}. Typically only used when calling \code{sop()}.
}

\usage{
sop.control(maxit = 200, epsilon = 1e-6, trace = FALSE)
}
\arguments{
  \item{maxit}{numerical value indicating the maximum number of iterations. Default set to 200 (see Details).}
  \item{epsilon}{numerical value indicating the tolerance for the convergence criterion. Default set to 1e-6 (see Details).}
  \item{trace}{logical indicating if output should be produced for each iteration.}
}
\details{
For Gaussian response variables, the implemented algorithm is an iterative procedure, with the fixed and random effects as well as the variance components being updated at each iteration. To check the convergence of this iterative procedure, the (REML) deviance is monitored.
For non-Gaussian response variables, estimation is based on Penalized Quasi-likelihood (PQL) methods. Here, the algorithm is a two-loop algorithm: the outer loop corresponds to the Fisher-Scoring algorithm (monitored on the basis of the change in the linear predictor between consecutive iterations), and the inner loop corresponds to that described for the Gaussian case. 
}
\value{
A list with the arguments as components.
}
\references{
Rodriguez-Alvarez, M.X., Lee, D. J., Kneib, T., Durban, M., and Eilers, P. (2015). Fast smoothing parameter separation in multidimensional generalized P-splines: the SAP algorithm. \emph{Statistics and Computing}, \bold{25 (5)}, 941--957.

Rodriguez-Alvarez, M.X., Durban, M., Lee, D. J. and Eilers, P. (2019). On the estimation of variance parameters in non-standard generalised linear mixed models: application to penalised smoothing. \emph{Statistics and Computing}, \bold{29 (3)}, 483--500.

}
\seealso{\code{\link{sop}}.}
\examples{
library(SOP)
# Simulate the data
set.seed(123)
n <- 1000
sigma <- 0.5
x <- runif(n)
f0 <- function(x) 2*sin(pi*x)
f <- f0(x)
y <- f + rnorm(n, 0, sigma)
dat <- data.frame(x = x, y = y)

# Fit the model
m0 <- sop(formula = y ~ f(x, nseg = 10), data = dat, control = list(trace = FALSE))
summary(m0)
}
\keyword{regression}
