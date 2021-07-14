####################################################################
####################################################################
sop.fit <- function(y, X, Z, weights = NULL, G = NULL,
                  vcstart = NULL, etastart = NULL, mustart = NULL, 
                  offset = NULL,
                  family = gaussian(), control = sop.control()) {
#------------------------------------------------------
# local function
deviance <- function(C, G, w, sigma2, ssr, edf) {
    log_det_C <- determinant(C)$modulus
    log_det_G <- determinant(G)$modulus
    deviance <- log_det_C + log_det_G + sum(log(sigma2*1/w)) + ssr/sigma2 + edf
    deviance  
  }  
#-----------------------------------------------------  	
	control <- do.call("sop.control", control)  
	if (missing(X))  
		stop("Missing argument: 'X' must be provided")
	if (missing(y))  
		stop("Missing argument: 'y' must be provided")
	if (missing(Z)) 
		stop("Missing argument: 'Z' must be provided")  #

	if(!is.null(vcstart)) {
		if(length(vcstart) != (length(G) + 1)) {
			stop("The length of 'vcstart' should be equal to the length of 'G' + 1)")  #
		}
	}
	trace <- control$trace

	ynames <- if (is.matrix(y)) { 
		rownames(y)
	} else {
		names(y)
	}

	conv <- FALSE
	nobs <- NROW(y)
	nvars <- NCOL(X)
	EMPTY <- nvars == 0
	nc <- ncol(X)
	ncomp <- length(G)
	nq <- ncol(Z)

	if(!is.null(vcstart)) {
		la <- vcstart
	} else {
		la <- rep(1, len = ncomp + 1)
	}

	devold <- 1e10
	if (is.null(weights)) {
		weights <- rep.int(1, nobs)
	}

	prior.weights <- weights

	if (is.null(offset)) { 
		offset <- rep.int(0, nobs)
	}

	if (is.null(G)) { 
		stop("Missing argument: 'G' must be provided")
	}
 
	Xnames <- colnames(X)
	na.ind <- is.na(y)
	y.tmp <- y
	y[na.ind] <- 1 # TODO: they are not included in the fit, but needed to initialise mu and eta (if negative or zero, problems)
	weights <- weights*(!na.ind)
	X <- as.matrix(X)

	start <- NULL
	variance <- family$variance
	linkinv <- family$linkinv
	linkfun <- family$linkfun
	mu.eta <- family$mu.eta
	
	if (!is.function(variance) || !is.function(linkinv)) { 
		stop("illegal `family' argument")
	}
	valideta <- family$valideta
	if (is.null(valideta)) {
 		valideta <- function(eta) TRUE
 	}
 	validmu <- family$validmu

	if (is.null(validmu)) { 
		validmu <- function(mu) TRUE
	}
	if (is.null(mustart)) {
		eval(family$initialize)
	} else {
		mukeep <- mustart
		eval(family$initialize)
		mustart <- mukeep
	}
	if (NCOL(y) > 1) { 
		stop("y must be univariate")
	}
	eta <- if (!is.null(etastart)) { 
		etastart
	} else if (!is.null(start)) { 
		if (length(start) != nvars) { 
			stop(gettextf("Length of start should equal %d and correspond to initial coefs.", nvars))
		} else {
			coefold <- start
			offset + as.vector(if (nvars == 1) 
			c(X,Z) * start
			else c(X,Z) %*% start)
		}
	} else {
		family$linkfun(mustart)
	}
 	
 	mu <- linkinv(eta)
	if (!(validmu(mu) && valideta(eta))) { 
		stop("Can't find valid starting values: please specify some")
	}

	for (i in 1:control$maxit) {  
		deriv <- family$mu.eta(eta)
		z <- (eta - offset) + (y - mu)/deriv
		w <- as.vector(deriv^2/family$variance(mu))
		w <- w*weights
		z[!weights] <- 0
		mat <- construct.matrices(X, Z, z, w)
		
		if(trace) start1 <- proc.time()[3]

		for (it in 1:control$maxit) {
			Ginv <- 0   
			for (ii in 1:ncomp) {
				Ginv <- Ginv + 1/la[ii+1]*G[[ii]]    
			} 
			GG <- 1/Ginv
			V <- construct.block(mat$XtX., mat$XtZ., mat$ZtX., mat$ZtZ.)
			D <- diag(c(rep(0,nc), Ginv))
			H <- (1/la[1])*V + D
			Hinv <- try(solve(H))

			if(inherits(Hinv,"try-error")) {
				Hinv <- ginv(H)
			}
			b <- (1/la[1])*Hinv%*%mat$u
			
			b.fixed  <- b[1:nc]
			b.random <- b[-(1:nc)]

			aux <- GG - diag(Hinv[-(1:nc),-(1:nc)])
	 		ed.sum <- 0  
			ied <- taus <- NULL
			for (i.fun in 1:ncomp) {
				d1.d <- (1/la[i.fun+1])*G[[i.fun]]
				ed1  <- sum(aux*d1.d)
				ed1  <- ifelse(ed1 <= 1e-10, 1e-10,ed1)
				tau1 <- sum(b.random^2*G[[i.fun]])/ed1
				tau1 <- ifelse(tau1 <= 1e-10, 1e-10, tau1)
				taus <- c(taus,tau1)
				ied  <- c(ied,ed1)
			}
			ssr <- mat$yty. - t(c(b.fixed, b.random)) %*% (2*mat$u - V %*% b)
			dev <- deviance(H, diag(GG), w[w != 0], la[1], ssr, sum(b.random^2*Ginv))[1]
			if(family$family == "gaussian" | family$family == "Gamma" | family$family == "quasipoisson") {
				# Sigma^2
				sig2 <- as.numeric((ssr / (length(y[weights != 0]) - sum(ied) - nc)))
			} else {
				sig2 <- 1
			}
			lanew <- c(sig2, taus)
			dla <- abs(devold - dev)			
			if(trace) {
				cat(sprintf("%1$3d %2$10.6f", it, dev))
				cat(sprintf("%8.3f", ied), '\n')
			}
			if (dla < control$epsilon) break			
			la <- lanew
			devold <- dev
		}
		if (trace) {
			end1 <- proc.time()[3]
			cat("Timings:\nSOP", (end1-start1), "seconds\n")
		}
		eta.old <- eta
		eta <- X%*%b.fixed + Z%*%b.random + offset	
		mu <- linkinv(eta)
		tol <- sum((eta - eta.old)^2)/sum(eta^2) 
		if (tol < control$epsilon | (family$family == "gaussian" & family$link == "identity")) break
	} 
	end <- proc.time()[3]
	mu.eta <- family$mu.eta
	mu.eta.val <- mu.eta(eta)
	linear.predictor <- eta
	mu <- linkinv(eta)
	names(mu) <- ynames
	names(linear.predictor) <- ynames
	
	# ojo offsetsssssssssssssssssss
	# Deviance residuals
	residuals <- family$dev.resids(y, mu, weights)
	s <- attr(residuals, "sign")
	if (is.null(s)) { 
		s <- sign(y - mu)
	}

	residuals <- sqrt(pmax(residuals, 0))*s
	names(residuals) <- ynames
	residuals[na.ind] <- NA

	names(b.fixed)<- Xnames
	names(b.random)<- colnames(Z)
  
	names(la) <- c("ssr",names(G))
	names(ied) <- c(names(G))
    
	dev.residuals <- family$dev.resids(y, mu, weights)
	dev.residuals[na.ind] <- NA
	deviance <- sum(dev.residuals, na.rm = TRUE)
	null.deviance <- glm(y ~ offset(offset), 
						family = family, 
						weights = prior.weights)$deviance
	out <- list("tol.ol" = tol, "it.ol" = i, "tol.il" = dla, "it.in" = it, "vc" = la, "edf" = ied)
	fit <- list()
	fit$b.fixed  <- b.fixed
	fit$b.random <- b.random 
	fit$fitted.values <- mu
	fit$linear.predictor <- linear.predictor
	fit$residuals <- residuals # Deviance residuals
	fit$X <- X
	fit$Z <- Z
	fit$G <- G
	fit$y <- y.tmp
	fit$weights <- weights
	fit$family <- family
	fit$out <- out 
	#fit$dev.residuals <- dev.residuals
	fit$deviance <- deviance
	fit$null.deviance <- null.deviance
	fit$Vp <- Hinv
	
	class(fit) <- "sop"
	invisible(fit)
}
####################################################################
####################################################################
construct.matrices <- function(X, Z, z, w) {
  XtW. <- t(X*w)
  XtX. <- XtW.%*%X
  XtZ. <- XtW.%*%Z
  ZtX. <- t(XtZ.)
  ZtW. <- t(Z*w)
  ZtZ. <- ZtW.%*%Z
  Xty. <- XtW.%*%z
  Zty. <- ZtW.%*%z
  yty. <- sum((z^2)*w)
  ZtXtZ <- rbind(XtZ., ZtZ.)
  u <- c(Xty.,Zty.)
  res <- list(XtX. = XtX., XtZ. = XtZ., ZtX. = ZtX., ZtZ. = ZtZ., Xty. = Xty., Zty. = Zty., yty. = yty., ZtXtZ = ZtXtZ, u = u)
}

construct.block <- function(A1,A2,A3,A4) {
  block <- rbind(cbind(A1,A2), cbind(A3,A4))
  return(block)
}
####################################################################
####################################################################