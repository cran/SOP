###################################################################
###################################################################
construct.1D.prediction.matrix <- function(object, newdata) {
	x.coord <- newdata[,object$vars]    
	B <- spline.bbase(object$Xmat$terms$MM$knots, x.coord, object$degree)
	X <- B%*%object$Xmat$terms$MM$U.X
	Z <- B%*%object$Xmat$terms$MM$U.Z
	X <- X[,-1,drop = FALSE]
	# Center matrices
	X <- sweep(X, 2, object$Xmat$cm$X)
	Z <- sweep(Z, 2, object$Xmat$cm$Z)    
	res <- list(X = X, Z = Z)
	res  
}
###################################################################
###################################################################
construct.2D.prediction.matrix <- function(object, newdata) {
	x.coord <- newdata[,object$vars[1]]
	y.coord <- newdata[,object$vars[2]]    
	B1p <- spline.bbase(object$Xmat$terms$MM$MM1$knots, x.coord, object$degree[1])
	B2p <- spline.bbase(object$Xmat$terms$MM$MM2$knots, y.coord, object$degree[2])    
	X1p <- B1p%*%object$Xmat$terms$MM$MM1$U.X
	X2p <- B2p%*%object$Xmat$terms$MM$MM2$U.X  
	Z1p <- B1p%*%object$Xmat$terms$MM$MM1$U.Z
	Z2p <- B2p%*%object$Xmat$terms$MM$MM2$U.Z    
	Xp <- Rten2(X2p, X1p)
	Xp <- Xp[,-1,drop = FALSE]    
	Zp <- cbind(Rten2(X2p, Z1p), Rten2(Z2p, X1p), Rten2(Z2p, Z1p))
	# Center matrices
	Xp <- sweep(Xp, 2, object$Xmat$cm$X)
	Zp <- sweep(Zp, 2, object$Xmat$cm$Z)
	res <- list(X = Xp, Z = Zp)
	res	
}
####################################################################
####################################################################
construct.3D.prediction.matrix <- function(object, newdata) {
	x.coord <- newdata[,object$vars[1]]
	y.coord <- newdata[,object$vars[2]]    
	z.coord <- newdata[,object$vars[3]]        
	B1p <- spline.bbase(object$Xmat$terms$MM$MM1$knots, x.coord, object$degree[1])
	B2p <- spline.bbase(object$Xmat$terms$MM$MM2$knots, y.coord, object$degree[2])
	B3p <- spline.bbase(object$Xmat$terms$MM$MM3$knots, z.coord, object$degree[3])
	X1p <- B1p%*%object$Xmat$terms$MM$MM1$U.X
	X2p <- B2p%*%object$Xmat$terms$MM$MM2$U.X  
	X3p <- B3p%*%object$Xmat$terms$MM$MM3$U.X  
	Z1p <- B1p%*%object$Xmat$terms$MM$MM1$U.Z
	Z2p <- B2p%*%object$Xmat$terms$MM$MM2$U.Z    
	Z3p <- B2p%*%object$Xmat$terms$MM$MM3$U.Z    
	rx12<-Rten2(X1p,X2p)  			
	Xp <- Rten2(rx12,X3p)    
	Xp <- Xp[,-1,drop = FALSE]    
	Zp <- cbind(Rten2(Rten2(Z1p,X2p), X3p), 
				Rten2(Rten2(X1p,Z2p), X3p),
				Rten2(rx12, Z3p),
				Rten2(Rten2(Z1p,Z2p), X3p),
				Rten2(Rten2(Z1p,X2p), Z3p),
				Rten2(Rten2(X1p,Z2p), Z3p),
				Rten2(Rten2(Z1p,Z2p), Z3p))
	# Center matrices
	Xp <- sweep(Xp, 2, object$Xmat$cm$X)
	Zp <- sweep(Zp, 2, object$Xmat$cm$Z)
	res <- list(X = Xp, Z = Zp)
	res  
  }
####################################################################
####################################################################
construct.random.prediction.matrix <- function(object, newdata) {
	if(!is.null(object$random)) {
		Zp <- NULL
		mfp <- model.frame(object$random$terms, newdata, xlev = attr(object$random$terms, "xlev"), na.action = na.pass)
		Zp <- model.matrix(object$random$terms, data = mfp, contrasts.arg = attr(object$random$terms, "contrast"))
		Zp <- Zp[,-1,drop = FALSE]
		Zp[is.na(Zp)] <- 0
	} else {
		Zp <- NULL
	}
	Zp
}
####################################################################
####################################################################
construct.fixed.prediction.matrix <- function(object, newdata) {
	if(!is.null(object$lin$terms)) {
		mfp <-model.frame(object$lin$terms, newdata, xlev = attr(object$lin$terms, "xlev"))
		Xp <- model.matrix(object$lin$terms, data = mfp, contrasts.arg = attr(object$lin$terms, "contrast"))
		Xp <- Xp[,-1,drop = FALSE]
	} else {
		Xp <- NULL
	}
	Xp
}
####################################################################
####################################################################
predict.sop <- function(object, newdata, type = c("response", "link", "terms"), se.fit = FALSE, ...) {
	type <- match.arg(type)
	na.action <- na.pass
	terms <- NULL
	if (missing(newdata)) {
		na.ind <- apply(is.na(object$dat[, object$model.terms, drop = FALSE]), 1, any) | is.na(object$y)
		if(type == "link") {
			if(se.fit) {
				pred <- list()
				pred$fit <- object$linear.predictor
				mm <- cbind(object$X, object$Z)
				se <- sqrt(colSums(t(mm) * tcrossprod(object$Vp, mm)))
				pred$se.fit <- insert.na(se, na.ind)
				#se.fit <- rep(NA, length = length(object$linear.predictor))
				#se.fit[!na.ind] <- se
				#pred$se.fit <- se.fit
			} else {
				pred <- object$linear.predictor
			}
			return(pred)
		} else if (type == "response") {
			if(se.fit) {
				pred <- list()
				pred$fit <- object$fitted.values
				mm <- cbind(object$X, object$Z)
				se <- sqrt(colSums(t(mm) * tcrossprod(object$Vp, mm)))
				pred$se.fit <- insert.na(se, na.ind)
				#se.fit <- rep(NA, length = length(object$fitted.values))
				#se.fit[!na.ind] <- se
				#pred$se.fit <- se.fit
			} else {
				pred <- object$fitted.values
			}
			return(pred)
		} else {
			newdata <- object$dat
			newdata.na <- object$dat[!na.ind,] # Exclude NAs in the covariates
		}
	} else {
		if (inherits(newdata, what = 'data.frame')) {
	        newdata <- as.data.frame(newdata)
	    } else {
	        stop("The object specified in argument 'data' is not a data frame")
	    }

		if(sum(is.na(match(object$model.terms, names(newdata)))))
	    	stop("Not all needed variables are supplied in newdata")

	    # Exclude NAs in the covariates
	    na.ind <- apply(is.na(newdata[, object$model.terms, drop = FALSE]), 1, any)	    
	    newdata.na <- newdata[!na.ind,, drop = FALSE]
	}

	cnames <- object$names.terms
	if (type == "terms") {
		if (is.null(terms)) {
			terms <- cnames
		} else { 
			if (sum(!(terms %in% cnames))) {
				stop("Non-existent terms requested")
			}
		}
	}
	b.fixed <- object$b.fixed
	b.random <- object$b.random
	nlin <- object$nterms[1]
	nre  <- object$nterms[2]
	nfun <- object$nterms[3]

	Z <- X <- NULL  
	pterms <- type == "terms"
	y.terms <- NULL
	y.terms.se <- NULL
	names.terms <- object$names.terms

	dim.fixed <- c(if(nlin > 0) object$lin$dim, if(nfun > 0) unlist(lapply(object$f, function(x) x$Xmat$dim$fixed)))
	dim.random <- c(if(nre > 0) object$random$dim, if(nfun > 0) unlist(lapply(object$f, function(x) x$Xmat$dim$random)))

	e.f <- cumsum(dim.fixed)
	s.f <- e.f - dim.fixed + 1

	e.r <- cumsum(dim.random)
	s.r <- e.r - dim.random + 1

	if (nlin > 0) {
		X <- construct.fixed.prediction.matrix(object, newdata.na)
		lab <- object$lin$vars
		if (pterms) {
			for (i.lin in 1:length(lab)) {
				if ((cnames[i.lin] %in% terms))  {
					vcum <- s.f[i.lin]:e.f[i.lin]
					mm <- X[,vcum, drop = FALSE]
					yt <- as.vector(mm%*%b.fixed[vcum + 1])
					y.terms <- cbind(y.terms, yt)
					if(se.fit) {
						yt.se <- sqrt(colSums(t(mm) * tcrossprod(object$Vp[vcum + 1,vcum + 1], mm)))
						y.terms.se <- cbind(y.terms.se, yt.se)
					}
				}
			}
		}
	}
	if (nre > 0) {
		aa <- construct.random.prediction.matrix(object, newdata.na)
		lab <- object$random$vars
		if (pterms) {
			for (i.re in 1:length(lab)) {
				if ((cnames[i.re + nlin] %in% terms))  {
					vcum <- s.r[i.re]:e.r[i.re]
					mm <- aa[,vcum, drop = FALSE]
					yt <- as.vector(mm%*%b.random[vcum])
					y.terms <- cbind(y.terms, yt)
					if(se.fit) {
						yt.se <- sqrt(colSums(t(mm) * tcrossprod(object$Vp[vcum + length(b.fixed), vcum + length(b.fixed)], mm)))
						y.terms.se <- cbind(y.terms.se, yt.se)
					}  
				}
			}
		}
		Z <- cbind(Z, aa)
	}  
	if (nfun > 0){    
		for (i.fun in 1:nfun) {			
			lab <- object$f[[i.fun]]$label
			aa <- switch(as.character(object$f[[i.fun]]$dim),
			   "1" = {
				aa <- construct.1D.prediction.matrix(object$f[[i.fun]], newdata.na)					
			}, "2" = {      
				aa <- construct.2D.prediction.matrix(object$f[[i.fun]], newdata.na)				
			}, "3" = {
				aa <- construct.3D.prediction.matrix(object$f[[i.fun]], newdata.na)				
			}
			)
			if (pterms && (cnames[i.fun + nre + nlin] %in% terms)) {
				if(ncol(aa$X) == 0) {
					vcum.f <- NULL
				} else {											
					vcum.f <- s.f[i.fun + nlin]:e.f[i.fun + nlin]
				}
				vcum.r <- s.r[i.fun + nre]:e.r[i.fun + nre]
				mm <- cbind(aa$X, aa$Z)
				pos.coeff <- c(vcum.f + 1, vcum.r + length(b.fixed))
				yt <- aa$X%*%b.fixed[vcum.f + 1] + aa$Z%*%b.random[vcum.r]
				y.terms <- cbind(y.terms,yt)
				if(se.fit) {
					yt.se <- sqrt(colSums(t(mm) * tcrossprod(object$Vp[pos.coeff,pos.coeff], mm)))
					y.terms.se <- cbind(y.terms.se, yt.se)
				}				
			}
			X <- cbind(X, aa$X)
			Z <- cbind(Z, aa$Z)
		} 
	}

	if (is.null(newdata.na$offset)) {
		offset <- 0
	} else {
		offset <- newdata.na$offset
	}
	if (length(b.fixed) > 1) {
		eta <- b.fixed[1] + X%*%b.fixed[-1] + Z%*%b.random + offset
	} else { 
		eta <- b.fixed[1] + Z%*%b.random + offset
	}
	if (type == "link") {
		yp <- eta
		if(se.fit) {
			pred <- list()
			pred$fit <- insert.na(yp, na.ind)
			mm <- as.matrix(cbind(rep(1, lenght = nrow(X)), X, Z)) #  Intercept
			se <- sqrt(colSums(t(mm) * tcrossprod(object$Vp, mm)))
			pred$se.fit <- insert.na(se, na.ind)
		} else {
			pred <- insert.na(yp, na.ind)
		}
	}
	if (type == "response") {
		yp <- object$family$linkinv(eta)
		if(se.fit) {
			pred <- list()
			pred$fit <- insert.na(yp, na.ind)
			mm <- as.matrix(cbind(rep(1, lenght = nrow(X)), X, Z)) # Intercept
			se <- (sqrt(colSums(t(mm) * tcrossprod(object$Vp, mm))))*abs(object$family$mu.eta(eta))
			pred$se.fit <- insert.na(se, na.ind)
		} else {
			pred <- insert.na(yp, na.ind)
		}
	}
	if (pterms) {
		iterms <- which(cnames %in% terms)
		y.terms <- insert.matrix.na(y.terms, na.ind)
		colnames(y.terms) <- object$names.terms[iterms]
		rownames(y.terms) <- rownames(newdata)
		if(se.fit) {
			y.terms.se <- insert.matrix.na(y.terms.se, na.ind)
			colnames(y.terms.se) <- object$names.terms[iterms]
			rownames(y.terms.se) <- rownames(newdata)
			pred <- list()
			pred$fit <- y.terms
			pred$se.fit <- y.terms.se 
		} else {
			pred <- y.terms
		}
		attr(pred, "constant") <- b.fixed[1]
		pred
	} else {
		drop(pred)
	}
}
####################################################################
####################################################################