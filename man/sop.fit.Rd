\name{sop.fit}
\alias{sop.fit}
\title{
Fitting generalised linear mixed models with overlapping precision matrices.
}
\description{
This is an internal function of package SOP. It is used to fit SOP models by specifying the design matrices for the fixed and random effects as well as the precision matrices for each variance component in the model.
}
\usage{
sop.fit(y, X, Z, weights = NULL, G = NULL, vcstart = NULL, 
  etastart = NULL, mustart = NULL, offset = NULL, 
  family = gaussian(), control = sop.control())     
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y}{vector of observations of length n.}
\item{X}{design matrix for the fixed effects (dimension n x p).}
\item{Z}{design matrix for the random effects (of dimension n x q).}
\item{weights}{an optional vector of 'prior weights' to be used in the fitting process. If NULL (default), the weights are considered to be one.}
\item{G}{a list with the diagonal elements of the precision matrices for each variance component in the model. Each element of the list is a vector of the same length as the number of columns in Z (i.e. q). The vector can be padded out with zeroes to indicate random coefficient not `affected' by the variance component (see details).}
\item{vcstart}{optional numeric vector. Initial values for the variance components (including the error variance as the first element of the vector). If \code{NULL}, all variance components are initialised to one.}
\item{etastart}{initial values for the linear predictor.}
\item{mustart}{initial values for the expected response.}
\item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of observations.}            
\item{family}{object of class \code{\link{family}} specifying the distribution and link function.} 
\item{control}{a list of control values to replace the default values returned by the function \code{\link{sop.control}.}}          
}
\details{
\code{sop.fit} is the workhorse function: it is typically not normally called directly but can be more efficient where the response vector 'y', design matrixs 'X' and 'Z', and precision matrices 'G' have already been calculated. Currently, the funcion only allows for diagonal precision matrices (possibly overlappping).
}
\value{
  A list containing the following objects: 
  \item{b.fixed}{the estimated fixed effect coefficients.}
  \item{b.random}{the predicted random effect coefficients.} 
  \item{residuals}{the (deviance) residuals.}
  \item{fitted.values}{the fitted values.}
  \item{linear.predictor}{the values of the linear predictor.}
  \item{X}{the fixed effect design matrix.}
  \item{Z}{the random effect design matrix.}
  \item{y}{the response.}
  \item{weights}{the prior weights.}
  \item{family}{the distribution family.}
  \item{out}{a list with i) \code{tol.ol} tolerance parameter (outer loop); ii) \code{it.ol} number of iteration (outer loop); iii) \code{tol.il} tolerance parameter (inner loop); \code{it.il} (number of iteration (inner loop)), iv) \code{vc} variance components estimates, v) \code{edf} effective degrees of freedom.}
  \item{deviance}{the deviance.}
  \item{null.deviance}{the null deviance.}
  \item{Vp}{Bayesian posterior covariance matrix for the coefficients.}
}

\references{
Rodriguez-Alvarez, M.X., Lee, D. J., Kneib, T., Durban, M., and Eilers, P. (2015). Fast smoothing parameter separation in multidimensional generalized P-splines: the SAP algorithm. \emph{Statistics and Computing}, \bold{25 (5)}, 941--957.

Rodriguez-Alvarez, M.X., Durban, M., Lee, D. J. and Eilers, P. (2019). On the estimation of variance parameters in non-standard generalised linear mixed models: application to penalised smoothing. \emph{Statistics and Computing}, \bold{29 (3)}, 483--500.
}

%\author{Maria Xose Rodriguez-Alvarez  \email{mxrodriguez@bcamath.org},  Manuel Oviedo de la Fuente.}

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
library(SOP)
# Simulate the data
set.seed(123)
n <- 1000
sigma <- 0.1
x1 <- runif(n, -1, 1)
x2 <- runif(n, -1, 1)
f0 <- function(x1, x2) cos(2*pi*sqrt((x1 - 0.5)^2 + (x2 - 0.5)^2))
f <- f0(x1, x2)
y <- f + rnorm(n, 0, sigma)
dat <- data.frame(x1 = x1, x2 = x2, y = y)

# Save but not fit the model
m0_nfit <- sop(formula = y ~ f(x1, x2, nseg = 10), data = dat, 
	    	fit = FALSE)

# Now fit using sop.fit()
m0 <- sop.fit(X = m0_nfit$X, Z = m0_nfit$Z, G = m0_nfit$G, 
	y = m0_nfit$y, weights = m0_nfit$weights, 
	   control = list(trace = FALSE))

names(m0)
}
% R documentation directory.
\keyword{regression}%
