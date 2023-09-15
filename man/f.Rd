\name{f}
\alias{f}
\title{Defining smooth terms in SOP formulae}
\description{Auxiliary function used to define smooth terms within \code{sop()} model formulae. The function does not evaluate the smooth - it exists purely to help set up a model using P-spline based smoothers.
}

\usage{f(..., nseg = 10, pord = 2 , degree = 3)}

\arguments{
\item{...}{a list of up to three variables to construct the smooth term.}
\item{nseg}{the number of segments for the (marginal) B-spline bases used to represent the smooth term. Numerical vector of length equal to the number of covariates. Atomic values are also valid, being recycled. The default value is 10.}
\item{pord}{penalty order. Numerical vector of length equal to the number of covariates. Atomic values are also valid, being recycled. The default value is 2 (second-order penalty).}
\item{degree}{the order of the polynomial for the (marginal) B-spline bases for this term. Numerical vector of length equal to the number of covariates. Atomic values are also valid, being recycled. The default value is 3 (cubic B-splines).}
}

\details{
The functions \code{f()} is designed to represent either a one dimensional smooth functions for main effects of an continuous explanatory variable or two or three dimensional smooth functions representing two way and three way interactions of continuous variables. By default, the values of the arguments  \code{nseg}, \code{pord} and \code{degree}  are repeated to the length of the explanatory covariates. The two and three dimensional smooth terms are constructed using the tensor-product of marginal (one-dimensional) B-spline bases and anisotropic penalties are considered.

%The contributions of the main effects and interactions can be seperated. For example the forumae \code{f(x1,x2)} and \code{f(x1)+f(x2)+f(x1,x2)} will fit essencialy the same model but in the later the contribution of the main effects and the interactions will beeen shown seperatly, i.e. the different effective degrees.   

}

\value{ 
The function is interpreted in the formula of a \code{sop} model and creates the right framework for fitting the smoother. List containing the following elements:

\item{vars}{names of the covariates involved in the smooth term.}

\item{nseg}{the number of segments for the (marginal) B-spline basis for each covariate.}

\item{pord}{the penalty order (numerical vector of length equal to the number of covariates).}

\item{degree}{the order of the polynomial for the (marginal) B-Spline bases for this term (numerical vector of length equal to the number of covariates). }

\item{dim}{The dimension of the smoother - i.e. the number of covariates that it is a function of.}

\item{label}{labels terms.}
}

%\author{ Maria Xose Rodriguez-Alvarez  \email{mxrodriguez@bcamath.org}, Manuel Oviedo de la Fuente.}
\seealso{\code{\link{ad}}, \code{\link{rae}}, \code{\link{sop}}}
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
m0 <- sop(formula = y ~ f(x, nseg = 10), data = dat)
summary(m0)

# Plot results
plot(y ~ x, data = dat)
ox <- order(dat$x)
lines(f[ox] ~ dat$x[ox], lwd = 2)
lines(fitted(m0)[ox] ~ dat$x[ox], col = "red", lwd = 2)
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..
