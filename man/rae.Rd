
\name{rae}
\alias{rae}
\title{Defining random effects in SOP formula}
\description{Auxiliary function used to define random effects terms in a \code{sop} model formula.}

\usage{rae(x)}

\arguments{
\item{x}{the x-variable (factor) that defines the random effects term.}
}

\details{
The functions is designed to represent random effects in SOP formulae. 
}

\value{The function is interpreted in the formula of a \code{sop} model and creates the right framework for fitting the random effect. List containing the following elements:
  \item{x}{name of the covariate involved.}
}

%\author{ Maria Xose Rodriguez-Alvarez  \email{mxrodriguez@bcamath.org}, Manuel Oviedo de la Fuente.}

\seealso{ \code{\link{f}}, \code{\link{ad}}, \code{\link{sop}}}

\examples{
library(SOP)
require(SpATS)
## An example of use of SOP package for the analysis of field trials experiments.
## Taken from the SpATS package.
data(wheatdata)

# Create factor variable for row and columns
wheatdata$R <- as.factor(wheatdata$row)
wheatdata$C <- as.factor(wheatdata$col)

# package SOP
m0 <- sop(formula = yield ~ colcode + rowcode + 
  f(col, row, nseg = c(10, 10)) + 
  rae(geno) + rae(R) + rae(C), data = wheatdata)
summary(m0)
plot(m0, col = topo.colors(100))
}

\keyword{models} \keyword{smooth} \keyword{regression}%-- one or more ..
