\name{summary.sop}
\alias{summary.sop}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary method for a fitted SOP model.
}
\description{
Summary method for a fitted SOP model.
}
\usage{
\method{summary}{sop}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{sop} as produced by \code{sop()}.}
  \item{\dots}{further arguments passed to or from other methods. Not yet implemented.}
}
%\details{
%\code{print.summary.lm} tries to be smart about formatting the coefficients.
%The coefficients component of the result gives the estimated coefficients 
%}
\value{
The function \code{summary.sop} computes and returns a list of summary statistics of the fitted model given in \code{object}, using the components (list elements) \code{"call"} and \code{"terms"} from its argument, plus

\item{call}{the matched call.}

%\item{residuals}{the \emph{deviance} residuals.}

\item{b.random}{a vector with the predicted random effects coefficients.}

\item{b.fixed}{a vector with the estimated fixed effects coefficients.}

\item{r.sq.adj}{the (adjusted) \eqn{R^2}, i.e., \sQuote{fraction of variance explained by the model}, \deqn{R^2 = 1 - \frac{\sum_i{R_i^2}/(n-df)}{\sum_i(y_i- y^*)^2/(n-1)},}{R^2 = 1 - Sum(R[i]^2)/(n-df) / Sum((y[i]- y*)^2)/(n-1),} where \eqn{R_i = w_i(y_i - \mu_i)} and \eqn{y^*}{y*} is the (weighted) mean of \eqn{y_i}{y[i]}.}

\item{deviance}{the deviance.}

\item{null.deviance}{the null deviance.}

\item{dev.expl}{proportion of the null deviance explained by the model.}

\item{n}{number of data.}

\item{iter}{number of iterations.}

\item{residual.df}{residual degrees of freedom.}

\item{edf}{a vector with the estimated degrees of freedom for the (smooth and random) model terms.}

\item{formula}{the model formula.}

\item{family}{the family used.}

 \item{na.action}{vector with the observations (position) deleted due to missingness.}

%a \eqn{p \times 4}{p x 4} matrix with columns for the estimated coefficient, its standard error, t-statistic and corresponding (two-sided) p-value.  Aliased coefficients are omitted.}

%\item{aliased}{named logical vector showing if the original coefficients are aliased.}

%\item{sigma}{the square root of the estimated variance of the random error \deqn{\hat\sigma^2 = \frac{1}{n-p}\sum_i{w_i R_i^2},}{\sigma^2 = 1/(n-p) Sum(w[i] R[i]^2),} where \eqn{R_i}{R[i]} is the \eqn{i}-th residual, \code{residuals[i]}.}
%\item{df}{degrees of freedom, a 3-vector \eqn{(p, n-p, p*)}, the first being the number of non-aliased coefficients, the last being the total number of coefficients.}

%\item{adj.r.squared}{the above \eqn{R^2} statistic \sQuote{\emph{adjusted}}, penalizing for higher \eqn{p}.}
}
\seealso{
\code{\link{sop}}, \code{\link{summary}}
}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
