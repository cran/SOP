\name{plot.sop}
\alias{plot.sop}
\title{Default SOP plotting}
\description{ Takes a fitted \code{sop} object produced by \code{sop()} and plots the component smooth functions that make it up, on the scale of the linear predictor. 
% Optionally produces term plots for parametric model components as well.
}
\usage{
\method{plot}{sop}(x, rug = TRUE, pages = 0, select = NULL, grid, ...) 
}
\arguments{ 
  
\item{x}{ a fitted \code{sop} object as produced by \code{sop()}.}

\item{rug}{when TRUE (default) then the covariate to which the plot applies is displayed as a rug plot at the foot of each plot of a 1-d smooth. Setting to \code{FALSE} will speed up plotting for large datasets.} 
  
%\item{se}{ when TRUE (default) upper and lower lines are added to the   1-d plots at 2 standard errors above and below the estimate of the smooth being plotted while for 2-d plots, surfaces at +1 and -1 standard errors are contoured and overlayed on the contour plot for the estimate. If a positive number is supplied then this number is multiplied by the standard errors when calculating standard error curves or surfaces. See also \code{shade}, below. }

\item{pages}{ (default 0) the number of pages over which to spread the output. For example, if \code{pages=1} then all terms will be plotted on one page with the layout performed automatically.  Set to 0 to have the routine leave all graphics settings as they are. }

\item{select}{Allows the  plot for a single model smooth term to be selected for printing. e.g. if you just want the plot for the second smooth term set \code{select = 2}. }

\item{grid}{number of covariate values used for each 1-d plot - for a nice smooth plot this needs to be several times the estimated degrees of freedom for the smooth. Default value 100.
}

%\item{n2}{Square root of number of points used to grid estimates of 2-d   functions for contouring.}

% \item{pers}{Set to \code{TRUE} if you want perspective plots for 2-d  terms.}

% \item{theta}{One of the perspective plot angles.}

%\item{phi}{The other perspective plot angle.}

%\item{jit}{Set to TRUE if you want rug plots for 1-d terms to be jittered.}

%\item{xlab}{If supplied then this will be used as the x label for all plots.}

%\item{ylab}{If supplied then this will be used as the y label for all plots.}

%\item{main}{Used as title (or z axis label) for plots if supplied.}

%\item{ylim}{If supplied then this pair of numbers are used as the y limits for each plot.}

%\item{xlim}{If supplied then this pair of numbers are used as the x limits for each plot.}

%\item{all.terms}{if set to \code{TRUE} then the partial effects of parametric model components are also plotted, via a call to \code{\link{termplot}}. Only terms of order 1 can be plotted in this way.}

\item{...}{ other graphics parameters to pass on to plotting commands. See details for smooth plot specific options.}

}
\details{Produces default plot showing the smooth and random components of a fitted SOP.

For smooth terms \code{plot.sop} actually calls plot method functions depending on the dimension of the smooth function.

For plots of smooths in one dimension, the x axis of each plot is labelled with the covariate name, while the y axis is labelled '\code{f(cov),edf}'  where \code{cov} is the covariate name, and \code{edf} the estimated degrees of freedom of the smooth.

Several smooth plots methods using \code{\link{image}} will accept a \code{colors} argument, which can be anything documented in \code{\link{topo.colors}} (in which case something like \code{colors=rainbow(50)} is appropriate), or the \code{\link{grey}} function (in which case something like \code{colors=grey(0:50/50)} is needed).
}

\value{The function main purpose is to generate plots. It also (silently) returns a list of the data used to produce the plots for the smooth terms. This function is inspired by the plot.gam function of the same name described in mgcv package (but is not a clone).
}
%\references{
%Rodriguez-Alvarez, M.X, Lee, D. J., Kneib, T., Durban, M., and Eilers, P. Fast smoothing parameter separation in multidimensional generalized P-splines: the SAP algorithm. \emph{Statistics and Computing}, \bold{25}(5), 941--957.}
%\author{ Maria Xose Rodriguez-Alvarez  \email{mxrodriguez@bcamath.org}, Manuel Oviedo de la Fuente%, Dea Jin Lee, Thomas Kneib, Maria Durban.}

\seealso{\code{\link{sop}}, \code{\link{predict.sop}}}

\examples{
library(SOP)
## Simulate the data
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
plot(m0)

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

m0 <- sop(formula = y ~ f(x1, x2, nseg = 10), data = dat, 
        control = list(trace = FALSE))

summary(m0)
plot(m0, col = topo.colors(100))

plot(m0, col = grey(0:100/100))

aux <- plot(m0)
names(aux)
}
\keyword{smooth} \keyword{regression} \keyword{hplot}


