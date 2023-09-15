###############################################################
###############################################################
plot.sop <- function(x, rug = TRUE, pages = 0, select = NULL, grid , ...) {
	lin.plot <- FALSE # Do not plot linear effects
	if(!is.null(select) && length(select) != 1) {
		stop("'select' arguments should be a vector of length one")
	}
	if (missing(grid)) {
		grid <- 100
	}


	b.fixed <- x$b.fixed
	b.random <- x$b.random
	
	nlin <- x$nterms[1]
	nre <- x$nterms[2]
	nfun <- x$nterms[3]

	dim.fixed <- c(if(nlin > 0) x$lin$dim, if(nfun > 0) unlist(lapply(x$f, function(x) sum(x$Xmat$dim$fixed))))
	dim.random <- c(if(nre > 0) x$random$dim, if(nfun > 0) unlist(lapply(x$f, function(x) x$Xmat$dim$random)))

	e.f <- cumsum(dim.fixed)
	s.f <- e.f - dim.fixed + 1

	e.r <- cumsum(dim.random)
	s.r <- e.r - dim.random + 1

	jit <- FALSE
	n.plots <- nfun + nre
	if(lin.plot) {
		n.plots <- n.plots + nlin
	}

 	if (n.plots == 0) { 
		stop("No terms to plot - nothing for plot.sop() to do.")
	}
	if (pages > n.plots) { 
		pages <- n.plots
	}
	if (pages < 0) {
		pages <- 0
	}
	if (pages != 0 & is.null(select)) {
		ppp <- n.plots%/%pages
		if (n.plots%%pages != 0) {
			ppp <- ppp + 1
			while (ppp * (pages - 1) >= n.plots) pages <- pages - 1   
		}
		c <- r <- trunc(sqrt(ppp))
		if (c < 1)       r <- c <- 1
		if (c * r < ppp) c <- c + 1
		if (c * r < ppp) r <- r + 1
		#oldpar <- par(mfrow = c(r, c))
		oldpar <- par(no.readonly = TRUE)
		on.exit(par(oldpar))
		par(mfrow = c(r, c))
	} #else {
		#oldpar <- par()
	#}
	if ((pages == 0 && prod(par("mfcol")) < n.plots && dev.interactive()) || pages > 1 && dev.interactive()) {
		ask <- TRUE
	} else { 
		ask <- FALSE
	}
	if (!is.null(select)) {
		ask <- FALSE
	}
	#if (ask) {
	#	oask <- devAskNewPage(TRUE)
	#	on.exit(devAskNewPage(oask))
	#}
	output <- list()
	names.output <- NULL
	cnames <- x$model.terms
	if (!is.null(select)) {
		terms <- cnames[select]
	} else {
		terms <- cnames
	}	
	if (nlin > 0 & lin.plot & is.null(select)) {
		lab <- x$lin$vars
		for (i.lin in 1:length(lab)) {
			if ((cnames[i.lin] %in% terms) & attr(x$lin$terms,"dataClasses")[i.lin] != "factor")  {
				x0 <- x$dat[,cnames[i.lin]]
				xseq <- seq(min(x0 ),max(x0), l = grid)
				vcum <- s.f[i.lin]:e.f[i.lin]
				yt <- b.fixed[vcum + 1]*xseq
				plot(xseq , yt, xlab = lab[i.lin], type = "l", ylab = paste("Partial for ", lab[i.lin], sep=""), ...)
				if (rug) {
					if (jit) 
						 rug(jitter(as.numeric(x$dat[,cnames[i.lin]])), ...)
					else rug(as.numeric(x$dat[,cnames[i.lin]]), ...)
				}
				if (ask) {
                	oask <- devAskNewPage(TRUE)
                	on.exit(devAskNewPage(oask), after = TRUE, add = TRUE)
                	ask <- FALSE
            	}
			}			
		}
	}
	if (nre > 0 & is.null(select)) {
		aa <- construct.random.prediction.matrix(x, x$data) 
		lab <- x$random$vars
		for (i.re in 1:length(lab)) {
			if ((cnames[i.re + nlin] %in% terms))  {
				vcum <- s.r[i.re]:e.r[i.re]
				lab.main <- paste0("rae(", lab[i.re], "), edf = ", round(sum(x$out$edf[names(x$out$edf) == lab[i.re]]), 2))
				qqnorm(b.random[vcum], main = lab.main, xlab = "Gaussian quantiles", ylab = "Effects",...)
				qqline(b.random[vcum])   
			}
			if (ask) {
                oask <- devAskNewPage(TRUE)
                on.exit(devAskNewPage(oask), after = TRUE, add = TRUE)
                ask <- FALSE
            }
		}
	}
	if (nfun > 0) { 
		for (i.fun in 1:nfun) {
			nseg <- x$f[[i.fun]]$nseg
			bdeg <- x$f[[i.fun]]$degree
			pord <- x$f[[i.fun]]$pord
			lab <- x$f[[i.fun]]$label
			if (x$f[[i.fun]]$dim == 1) {   
				if ((is.null(select) || i.fun == select)) {
					names.output <- c(names.output, x$f[[i.fun]]$label)
					lab <- paste0(x$f[[i.fun]]$label,", edf = ", round(sum(x$out$edf[names(x$out$edf) == x$f[[1]]$Xmat$edflabel]) + x$f[[i.fun]]$Xmat$dim$fixed, 2))    
					xlab <- x$f[[i.fun]]$vars
					x0 <- x$dat[,xlab]
					xseq <- seq(min(x0), max(x0), l = grid)
					Bp <- spline.bbase(x$f[[i.fun]]$Xmat$terms$MM$knots, xseq, bdeg)
					Xp <- Bp%*%x$f[[i.fun]]$Xmat$terms$MM$U.X
					Zp <- Bp%*%x$f[[i.fun]]$Xmat$terms$MM$U.Z
					Xp <- Xp[,-1, drop = FALSE]
					
					# Center matrices
					Xp <- sweep(Xp, 2, x$f[[i.fun]]$Xmat$cm$X)
					Zp <- sweep(Zp, 2, x$f[[i.fun]]$Xmat$cm$Z)

					if(ncol(Xp) == 0) {
						vcum.f <- NULL
					} else {											
						vcum.f <- s.f[i.fun + nlin]:e.f[i.fun + nlin]
					}
					vcum.r <- s.r[i.fun + nre]:e.r[i.fun + nre]

					# Standard errors
					pos.coeff <- c(vcum.f + 1, vcum.r + length(b.fixed))
					mm <- cbind(Xp, Zp)
					se <- sqrt(colSums(t(mm) * tcrossprod(x$Vp[pos.coeff,pos.coeff], mm)))

					output[[i.fun]] <- list()
					output[[i.fun]]$fit <- Xp%*%b.fixed[vcum.f + 1] + Zp%*%b.random[vcum.r]
					output[[i.fun]]$se <- se
					res <- c(output[[i.fun]]$fit, output[[i.fun]]$fit + 1.96*se, output[[i.fun]]$fit - 1.96*se)
					range <- c(min(res) - 0.1*abs(min(res)), max(res) + 0.1*abs(max(res)))

					plot(xseq, output[[i.fun]]$fit, xlab = xlab, ylab = lab, type = "l", ylim = range, ...)
					lines(xseq, output[[i.fun]]$fit + 1.96*se, lty = 2)
					lines(xseq, output[[i.fun]]$fit - 1.96*se, lty = 2)
					if (rug) {
						if (jit) 
							rug(jitter(as.numeric(x$dat[,xlab])), ...)
						else rug(as.numeric(x$dat[,xlab]), ...)
					}
					if (ask) {
                		oask <- devAskNewPage(TRUE)
                		on.exit(devAskNewPage(oask), after = TRUE, add = TRUE)
                		ask <- FALSE
            		}
				}
			} else {
				if (x$f[[i.fun]]$dim == 2) { 
					if ((is.null(select) || i.fun == select)) {
						names.output <- c(names.output, lab)
						xlab <- x$f[[i.fun]]$vars
						x1 <- x$dat[,xlab[1]]
						x2 <- x$dat[,xlab[2]]
						x1seq <- seq(min(x1), max(x1), l = grid)
						x2seq <- seq(min(x2), max(x2), l = grid)
						B1p <- spline.bbase(x$f[[i.fun]]$Xmat$terms$MM$MM1$knots, x1seq, bdeg[1])
						B2p <- spline.bbase(x$f[[i.fun]]$Xmat$terms$MM$MM2$knots, x2seq, bdeg[2])    
						X1p <- B1p%*%x$f[[i.fun]]$Xmat$terms$MM$MM1$U.X
						X2p <- B2p%*%x$f[[i.fun]]$Xmat$terms$MM$MM2$U.X  
						Z1p <- B1p%*%x$f[[i.fun]]$Xmat$terms$MM$MM1$U.Z
						Z2p <- B2p%*%x$f[[i.fun]]$Xmat$terms$MM$MM2$U.Z    
						Xp <- X2p%x%X1p
						Xp <- Xp[,-1,drop = FALSE]    
						Zp <- cbind(X2p%x%Z1p, Z2p%x%X1p, Z2p%x%Z1p)
						
						# Center matrices
						Xp <- sweep(Xp, 2, x$f[[i.fun]]$Xmat$cm$X)
						Zp <- sweep(Zp, 2, x$f[[i.fun]]$Xmat$cm$Z)
						
						if(ncol(Xp) == 0) {
							vcum.f <- NULL
						} else {											
							vcum.f <- s.f[i.fun + nlin]:e.f[i.fun + nlin]
						}
						vcum.r <- s.r[i.fun + nre]:e.r[i.fun + nre]

						# Standard errors
						pos.coeff <- c(vcum.f + 1, vcum.r + length(b.fixed))
						mm <- cbind(Xp, Zp)
						se <- sqrt(colSums(t(mm) * tcrossprod(x$Vp[pos.coeff,pos.coeff], mm)))

						p <- Xp%*%b.fixed[vcum.f + 1] + Zp%*%b.random[vcum.r]
						output[[i.fun]] <- list()
						output[[i.fun]]$fit <-  matrix(p, grid, grid)
						output[[i.fun]]$se <-  matrix(se, grid, grid)
						image(x1seq, x2seq, output[[i.fun]]$fit, main = lab, xlab = xlab[1], ylab = xlab[2], ...)

						if (ask) {
                			oask <- devAskNewPage(TRUE)
                			on.exit(devAskNewPage(oask), after = TRUE, add = TRUE)
                			ask <- FALSE
            			}
					}
				} else {
					if ((is.null(select) || i.fun == select)) {
						print("No 3D plot implemented")
					}
				}
			}
		}
	}
	#if (pages > 0 & is.null(select)) 
    #    par(oldpar)

    names(output) <- names.output
	invisible(output)
}


