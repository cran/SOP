\name{print.sop}
\alias{print.sop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print method for sop objects
}
\description{
Print method for sop objects
}
\usage{
\method{print}{sop}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{x}{an object of class \code{sop} as produced by \code{sop()}}
  	\item{\dots}{further arguments passed to or from other methods. Not yet implemented.}
}
\value{
Prints some summary statistics of the fitted model.
}
\seealso{
\code{\link{sop}}
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
m0
}

\keyword{smooth}
