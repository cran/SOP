\name{sop}
\alias{sop}
\title{
Estimation of generalised additive P-spline regression models with overlapping penalties.
}
\description{
The function \code{sop()} fits generalised additive regression models. For the smooth terms, it uses P-splines (Eilers and Marx, 1996) and it can cope with one, two and three dimensional smooth terms. The innovation of the function is that smoothing/variance parameters are estimated on the basis of the SOP method; see Rodriguez-Alvarez \emph{et al.} (2015) and Rodriguez-Alvarez \emph{et al.} (2019) for details. This speeds up the fit. 
}
\usage{
sop(formula, data = list(),  family = gaussian(), weights = NULL, offset = NULL, 
    control = sop.control(), fit = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{formula}{a sop formula. This is exactly like the formula for a GLM except that (1) P-splines in one, two and three dimensions (\code{\link{f}}), (2)spatially adaptive P-splines in 1 dimension (\code{\link{ad}}); and (3) random effects (\code{\link{rae}}) can be added to the right hand side of the formula.}
\item{data}{a data frame containing the model response variable and covariates required by the formula.} 
\item{family}{object of class \code{\link{family}} specifying the distribution and link function.} 
\item{weights}{prior weights on the contribution of the data to the log likelihood. Note that a weight of 2, for example, is equivalent to having made exactly the same observation twice. If you want to reweight the contributions of each datum without changing the overall magnitude of the log likelihood, then you should normalize the weights e.g. \code{weights <- weights/mean(weights)}). If NULL (default), the weights are considered to be one.}
\item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of observations. }                  
%\item{etastart}{Initial values for the linear predictor}
%\item{mustart}{Initial values for the expected response.}
\item{control}{a list of control values to replace the default values returned by the function \code{\link{sop.control}.}}
\item{fit}{logical. If TRUE, the model is fitted.}
}
\details{
The \code{sop()} can be used to fit generalised additive models. It works similarly to the function  \code{gam()}  of the package \pkg{mgcv}. The function \code{sop()} uses P-splines (Eilers and Marx, 1996), one of the option on \code{gam()}. Estimation is based on the equivalence between P-splines and linear mixed models, and variance/smoothing parameters are estimated based on restricted maximum likelihood (REML) using the separation of overlapping precision matrices (SOP) method described in Rodriguez-Alvarez \emph{et al.} (2015) and Rodriguez-Alvarez \emph{et al.} (2019).
The function \code{sop()} can be seen as a faster alternative to \code{gam()} for some data sets.  
}
\value{
  An object of class \sQuote{sop}. 
  It is a list containing the following objects: 
  \item{b.fixed}{the estimated fixed effect coefficients (present if \code{fit = TRUE}).}
  \item{b.random}{the predicted random effect coefficients (present if \code{fit = TRUE}).}   
  \item{fitted.values}{the fitted values (present if \code{fit = TRUE}).}
  \item{linear.predictor}{the values of the linear predictor (present if \code{fit = TRUE}).}
  \item{residuals}{the (deviance) residuals (present if \code{fit = TRUE}).}
  \item{X}{the fixed effect design matrix.}
  \item{Z}{the random effect design matrix.}
  \item{G}{a list containing information about the precision/penalty matrices (one for each smoothing/variance parameter in the model).}
  \item{y}{the response}
  \item{weights}{the prior weights.}
  \item{family}{the distribution family.}
  \item{out}{a list with i) \code{tol.ol} tolerance parameter (outer loop); ii) \code{it.ol} number of iteration (outer loop); iii) \code{tol.il} tolerance parameter (inner loop); \code{it.il} (number of iteration (inner loop)), iv) \code{vc} variance components estimates, v) \code{edf} effective degrees of freedom (present if \code{fit = TRUE}).}
   \item{deviance}{the deviance (present if \code{fit = TRUE}).}
   \item{null.deviance}{the null deviance (present if \code{fit = TRUE}).}
   \item{Vp}{Bayesian posterior covariance matrix for the coefficients (present if \code{fit = TRUE}).}
   \item{call}{the function call.}
   \item{data}{the data.}
   \item{formula}{the model formula.}
   \item{lin}{a list containing information about the parametric/linear.}
   \item{random}{a list containing information about the random effects.}
   \item{f}{a list containing information about the smoothers.}       
   \item{na.action}{vector with the observations (position) deleted due to missingness.}
   \item{names.terms}{the terms used in the formula.}
   \item{model.terms}{the explanatory variables.}   
   \item{nterms}{the number of linear, random and smooth terms in the formula.}
}

\references{
Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with B-splines and penalties. \emph{Statistical Science}, \bold{11 (2)}, 89--121.

Rodriguez-Alvarez, M.X., Lee, D. J., Kneib, T., Durban, M., and Eilers, P. (2015). Fast smoothing parameter separation in multidimensional generalized P-splines: the SAP algorithm. \emph{Statistics and Computing}, \bold{25 (5)}, 941--957.

Rodriguez-Alvarez, M.X., Durban, M., Lee, D. J. and Eilers, P. (2019). On the estimation of variance parameters in non-standard generalised linear mixed models: application to penalised smoothing. \emph{Statistics and Computing}, \bold{29 (3)}, 483--500.
}

%\author{Maria Xose Rodriguez-Alvarez  \email{mxrodriguez@bcamath.org},  Manuel Oviedo de la Fuente.}

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}

\examples{
library(SOP)
## An example of use of SOP package with tensor product B-splines in 2D
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

# Theoretical surface
np <- 50
x1p <- seq(-1, 1, length = np)  
x2p <- seq(-1, 1, length = np)
fp <- cos(2 * pi * sqrt(outer((x1p - 0.5) ^ 2, (x2p - 0.5) ^ 2, '+')))

image(x1p, x2p, matrix(fp, np, np), main = 'f(x1,x2) - Theor', 
     col = topo.colors(100))

# Fit the model
m0 <- sop(formula = y ~ f(x1, x2, nseg = 10), data = dat, 
        control = list(trace = FALSE))

summary(m0)
plot(m0, col = topo.colors(100))

## An example of use of SOP package with several smooth terms and Gamma distribution
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
m1 <- sop(formula = y ~ f(x0, nseg = 17) + 
                             f(x1, nseg = 17) + 
                             f(x2, nseg = 17) + 
                             f(x3, nseg = 17) + 
                             f(x4, nseg = 17) + 
                             f(x5, nseg = 17), family = Gamma(link = log), data = df)
summary(m1)
plot(m1)

## An example of use of SOP package for the analysis of field trials experiments.
## Taken from the SpATS package.
require(SpATS)
data(wheatdata)

# Create factor variable for row and columns
wheatdata$R <- as.factor(wheatdata$row)
wheatdata$C <- as.factor(wheatdata$col)

# package SOP
m2 <- sop(formula = yield ~ colcode + rowcode + 
  f(col, row, nseg = c(10, 10)) + 
  rae(geno) + rae(R) + rae(C), data = wheatdata)
summary(m2)
plot(m2, col = topo.colors(100), pages = 1)

# Package SpATS: more adequate for this analysis.
# SpATS has been explicitly developed for the analysis field trials experiments. 
m3 <- SpATS(response = "yield", 
    spatial = ~ SAP(col, row, nseg = c(10,10), degree = 3, pord = 2, center = TRUE), 
    genotype = "geno",
    genotype.as.random = TRUE,
    fixed = ~ colcode + rowcode, random = ~ R + C, data = wheatdata, 
    control =  list(tolerance = 1e-06))
summary(m3)
plot(m3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}%
