summary.sop <- function(object, ...) {
	edf <- object$out$edf
	nam <- names(edf)
	edf <- c(edf, sum(edf), sum(edf) + length(object$b.fixed))
	names(edf) <- c(nam, "Total edf", "Total")
	nobs <- length(object$y) - length(object$na.action)  
	residual.df <- nobs - sum(edf) - length(object$b.fixed)
	w <- object$weights
	mean.y <- sum(w * object$y, na.rm = TRUE)/sum(w)
	w <- sqrt(w)
	r.sq <- if(inherits(object$family, "general.family") || !is.null(object$family$no.r.sq)) { 
		NULL
	} else {
		1-var(w*(as.numeric(object$y) - object$fitted.values), na.rm = TRUE)*(nobs - 1)/(var(w*(as.numeric(object$y) - mean.y), na.rm = TRUE)*residual.df)
    }
	dev.expl <- (object$null.deviance - object$deviance)/object$null.deviance    
	out <- list(call = object$call, b.random = object$b.random, 
				b.fixed = object$b.fixed, r.sq.adj = r.sq, deviance = object$deviance, 
				null.deviance = object$null.deviance,
				dev.expl = dev.expl, n = nobs, iter = object$out$it.ol,
				residual.df = residual.df, edf = edf,
				formula = object$formula, family = object$family,
				na.action = object$na.action)
	
    class(out) <- "summary.sop"
	out
}
#------------------------------------------------------------
#------------------------------------------------------------
#------------------------------------------------------------
# deviance.sop <- function(object,...)
# {
#   sum(object$dev.residuals)  
# }
#------------------------------------------------------------
#------------------------------------------------------------
# AIC.sop <- function (object, ..., k = 2, c = FALSE) 
# {
#   if (length(list(...))) {
#     object <- list(object, ...)
#     issop <- unlist(lapply(object, is, class2="sop"))
#     if (!any(issop)) 
#       stop("some of the objects are not sop")
#     
#     df <- as.numeric(lapply(object, function(x) x$df.fit))
#     N <- as.numeric(lapply(object, function(x) x$N))
#     Cor <- if ((k == 2) && (c == TRUE)) 
#       (2 * df * (df + 1))/(N - df - 1)
#     else rep(0, length(object))
#     AIC <- as.numeric(lapply(object, function(x) x$G.dev + 
#                                x$df.fit * k)) + Cor
#     val <- as.data.frame(cbind(df, AIC))
#     Call <- match.call()
#     Call$k <- NULL
#     Call$c <- NULL
#     row.names(val) <- as.character(Call[-1])
#     val <- val[order(AIC), ]
#     val
#   }
#   else {
#     val <- if (is.gamlss(object)) 
#       object$G.dev + object$df.fit * k
#     else stop(paste("this is not a gamlss object"))
#     if ((k == 2) && (c == TRUE)) 
#       val <- val + (2 * object$df.fit * (object$df.fit + 
#                                            1))/(object$N - object$df.fit - 1)
#     val
#   }
# }
