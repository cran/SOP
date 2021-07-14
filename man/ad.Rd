\name{ad}
\alias{ad}
\title{Adaptive smooth terms in a SOP model formula}

\description{Auxiliary function used to define adaptive smooth terms in a \code{sop} model formula. The function does not do any fitting but is used as part of a \code{sop()} model formula.}

\usage{
ad(..., nseg = 10, pord = 2, degree = 3, nseg.sp = 5, degree.sp = 3)
}

\arguments{
  \item{\dots}{the x-variable (continuous) to be used for adaptive smoothing. Currently, only one dimensional adaptive smoothers are allowed.}
  \item{nseg}{the number of segments for the B-spline basis used to represent the smooth term. The default value is 10.}
  \item{pord}{penalty order. The defauly value is 2 (second order penalty).}
  \item{degree}{the order of the polynomial for the B-spline basis for this term. The default value is 3 (cubic B-splines).}
  \item{nseg.sp}{the number of segments for the B-spline basis used to `smooth' the smoothing parameters. The default value is 5.}
  \item{degree.sp}{the order of the polynomial for the B-spline basis used for smoothing of the smoothing parameters. The default value is 3 (cubic B-splines).}
}
\details{
The function \code{ad()} can be use to fit an adaptive smooth function of \code{x}.
An `adaptive'  smooth function is one in which the smoothing parameter is allowed to vary over the range of the explanatory variable \code{x}. Details can be found in Rodriguez-Alvarez \emph{at. al} (2019).  
}
\value{The function is interpreted in the formula of a \code{sop} model and creates the right framework for fitting the adaptive smoother. List containing the following objects:
\item{vars}{name of the covariates involved in the adaptive smooth term.}

\item{nseg}{the number of segments for the B-spline basis.}

\item{pord}{the penalty order.}

\item{degree}{the order of the polynomial for the B-Spline basis for this term. }

\item{nseg.sp}{the number of segments for the B-spline basis used to `smooth' the smoothing parameters.}

\item{degree.sp}{the order of the polynomial for the B-spline basis used for smoothing of the smoothing parameters.}

\item{dim}{The dimension of the smoother - i.e. the number of covariates that it is a function of.}

\item{label}{labels terms.}    
}
\references{
%Rodriguez-Alvarez, M.X., Lee, D. J., Kneib, T., Durban, M., and Eilers, P. (2015). Fast smoothing parameter separation in multidimensional generalized P-splines: the SAP algorithm. \emph{Statistics and Computing}, \bold{25 (5)}, 941--957.

Rodriguez-Alvarez, M.X., Durban, M., Lee, D. J. and Eilers, P. (2019). On the estimation of variance parameters in non-standard generalised linear mixed models: application to penalised smoothing. \emph{Statistics and Computing}, \bold{29 (3)}, 483--500.
}
%\author{ Maria Xose Rodriguez-Alvarez  \email{mxrodriguez@bcamath.org}, Manuel Oviedo de la Fuente.}
\seealso{\code{\link{f}}, \code{\link{rae}}, \code{\link{sop}}}
\examples{
library(SOP)
# Simulate the data
set.seed(123)
n <- 1000
x <- runif(n, 0.0001, 1)
doppler.function <- function(x) sin(4 / (x + 0.1)) + 1.5
mu <- doppler.function(x)
sigma <- 0.2
y <- mu + sigma*rnorm(n)
dat <- data.frame(x = x, y = y)

# Fit the models
# With addaptive smoothing
m0 <- sop(formula = y ~ ad(x, nseg = 197, nseg.sp = 17), data = dat, 
         control =  list(trace = FALSE, epsilon = 1e-03))

# Without addaptive smoothing
m1 <- sop(formula = y ~ f(x, nseg = 197), data = dat, 
         control =  list(trace = FALSE, epsilon = 1e-03))         

# Plot results
plot(y ~ x, data = dat)
ox <- order(dat$x)
lines(fitted(m0)[ox] ~ dat$x[ox], col = 2, lwd = 2)
lines(fitted(m1)[ox] ~ dat$x[ox], col = 4, lwd = 2)
legend("topright", c("Theoretical", "Adaptive", "Non Adaptive"), 
    col = c(1,2,4), lty = 1, lwd = 2, bty = "n")	
}
\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more .