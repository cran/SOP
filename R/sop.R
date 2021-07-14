####################################################################
####################################################################
sop <- function(formula, data = list(),  family = gaussian(), weights = NULL, offset = NULL, control = sop.control(), fit = TRUE) {
	control <- do.call("sop.control", control)

	if (inherits(data, what = 'data.frame')) {
        data <- as.data.frame(data)
    } else {
        stop("The object specified in argument 'data' is not a data frame")
    }

    nobs <- nrow(data)
	if (is.null(offset)) {
		offset <- rep.int(0, nobs)
	}
	if (is.null(weights)) {
		weights <- rep.int(1, nobs)
	}

	tf <- terms.formula(formula, specials = c("f", "ad","rae"))
	terms <- attr(tf, "term.labels")
	nt <- length(terms)     
	if (attr(tf, "response") > 0) {
		response <- as.character(attr(tf, "variables")[2])
	}
	ifun <- sort(c(attr(tf,"specials")$f, attr(tf,"specials")$ad))
	ire  <- attr(tf,"specials")$rae
	if(response != 0) { 
		ifun <- ifun - 1
		ire <- ire - 1
	}
	nfun <- length(ifun)
	smooth <- terms[ifun]

	nre <- length(ire)
	random <- terms[ire]

	if(length(ifun) || length(ire)) {
		fixed <- terms[-c(ifun,ire)]
		ilin <- (1:nt)[-c(ifun,ire)]
		nlin <- length(fixed)
	} else {
		fixed <- terms
		ilin <- 1:nt
		nlin <- length(fixed)
	}

	nterms <- c(nlin,nre,nfun)
	names(nterms) <- c("linear", "random", "smooth")

	# Check if this is needed
	sopenv <- environment(formula)
	sopns <- loadNamespace("SOP")

	# Obtain variables
	vars.formula <- NULL
	for(i in 1:length(terms)) {
		if(nfun != 0 && i %in% ifun) {
			vars.formula <- c(vars.formula, eval(parse(text = terms[i]), enclos = sopenv, envir = sopns)$vars)
		} else if (nre != 0 && i %in% ire) {
			vars.formula <- c(vars.formula, as.character(eval(parse(text = terms[i]), enclos = sopenv, envir = sopns)))
		} else if (nlin != 0 && i %in% ilin) {
			vars.formula <- c(vars.formula, all.vars(formula(paste("~", paste(terms[ilin], collapse = "+")))))
		}
	}
	vars.formula <- unique(vars.formula)

	if(sum(is.na(match(c(response, vars.formula), names(data)))))
    	stop("Not all needed variables are supplied in data")

	l.f <- l.lin <- l.re <- NULL
	np <- Z <- X <- NULL
	gg <- init.var <- dim <- list()
	var.re <- character(nre)
	Xnames <- "(Intercept)"
	model.terms <- NULL
	names.terms <- terms[c(ilin,ire,ifun)]
	
	# TRUE/FALSE for observations with non permitted NA (NAs in the covariates and/or response)
	na.ind <- apply(is.na(data[, c(response, vars.formula), drop = FALSE]), 1, any)
	weights <- weights*(!na.ind)
	data.na <- data[!na.ind,]
	weights.na <- weights[!na.ind]
	offset.na <-  offset[!na.ind]
	y <- data.na[,response]
	na.action <- (1:nrow(data))[na.ind] # Observations deleted due to missingness	
	
	########################## Linear effects and factors ############################
	if (nlin != 0) {
		var.lin <- terms[ilin]
		formula.lin <- formula(paste("~", paste(var.lin, collapse = "+")))
		l.lin <- construct.fixed.part(formula.lin, data.na)
		l.lin$vars <- var.lin
		l.lin$form   <- formula.lin  
		
		X <- cbind(X, l.lin$X)
		Xnames <- c(Xnames, colnames(l.lin$X))
		model.terms <- c(model.terms, var.lin)
		attr(names.terms,"type") <- c(attr(names.terms,"type"), attr(l.lin$terms,"dataClasses"))
	} else {
		l.lin <- NULL
	}
	
	########################## Random effects #######################################
	if (nre != 0) {
		for (i in 1:nre) {
    		var.re[i] <- as.character(eval(parse(text = terms[ire[i]]), enclos = sopenv, envir = sopns))
  		}
		formula.re <- formula(paste("~", paste(var.re, collapse = "+")))
		l.re <- construct.random.part(formula.re, data.na)
		l.re$vars <- var.re
		l.re$form <- formula.re

		np <- c(np, l.re$dim)
		Z <- cbind(Z, l.re$Z)
		init.var <- c(init.var, list(l.re$init.var))
		dim <- c(dim, list(l.re$dim))
		gg <- list(l.re$g)
		model.terms <- c(model.terms, var.re)
		attr(names.terms,"type") <- c(attr(names.terms,"type"),rep("random", length = nre))
	} else {
		l.re <- NULL
	}
	
	############################# Smooth functions ################################  
	if (nfun > 0){
		names.fun <- NULL
		for (i in 1:nfun) {
			l.f[[i]] <- eval(parse(text = terms[ifun[i]]), enclos = sopenv, envir = sopns)
			names.fun <- c(names.fun, terms[ifun[i]])
			form <- as.formula(paste("~", terms[ifun[i]], sep = ""))
			aa <- switch(as.character(l.f[[i]]$dim),
			"3" = {  				
				l.f[[i]]$Xmat <- construct.3D.pspline(form, data.na)
			}, "2"= {            				
				l.f[[i]]$Xmat <- construct.2D.pspline(form, data.na)
			}, "1" = {				
				l.f[[i]]$Xmat <- construct.1D.pspline(form, data.na)
			})

			np <- c(np, l.f[[i]]$Xmat$dim$random)		
			X <- cbind(X,l.f[[i]]$Xmat$X)		
			Z <- cbind(Z,l.f[[i]]$Xmat$Z)
			gg <- c(gg, list(l.f[[i]]$Xmat$g))										
			init.var <- c(init.var, list(l.f[[i]]$Xmat$init.var))
			dim <- c(dim, list(l.f[[i]]$Xmat$dim))
			Xnames <- c(Xnames, colnames(l.f[[i]]$Xmat$X))	
			model.terms <- c(model.terms, l.f[[i]]$vars)
			attr(names.terms,"type") <- c(attr(names.terms,"type"), rep("smooth",length = length(l.f[[i]]$label)))
			names(l.f) <- names.fun			
		}
	}
	if(is.null(X) & is.null(Z)) {
		stop("Neither fixed nor random or smooth effects have been specified in formula.")
	}
	if(is.null(Z)) {
		stop("Neither random nor smooth effects have been specified in formula. Currently, a model with only fixed effects cannot be estimated using sop()")
	}
	if(!is.null(X)) {
		X <- as.matrix(cbind(rep(1, lenght = nrow(X)), X)) #  Intercept
	} else {
		X <- matrix(1, ncol = 1, nrow = nrow(Z)) #  Intercept
	}
	colnames(X) <- Xnames
	np.names <- names(np)
	np <- c(ncol(X), np)
	names(np) <- c("fixed effects", np.names)
	if (length(gg) == 1 & length(gg[[1]]) == 1) {
		G <- list(unlist(gg))
		names(G) <- names(gg[[1]])
	} else {
		G <- construct.capital.lambda(gg)
	}
	if(fit) {
		mustart <- etastart <- NULL
		res <- sop.fit(y, X, Z, weights = weights.na, G = G, vcstart = c(1, unlist(init.var)),
					etastart = etastart, mustart = mustart,
					offset = offset.na, 
					family = family, control = control)

		# NAs for those observation in the dataset not included in the fit (NAs in the covariate values and/or response)
		res$fitted.values <- insert.na(res$fitted.values, na.ind)
		res$linear.predictor <- insert.na(res$linear.predictor, na.ind)
		res$residuals <- insert.na(res$residuals, na.ind)
		res$y <- data[,response] # TODO: including NAs
		res$weights <- weights # TODO: including NAs (set to zero)
		tmp <- add_coeff_to_terms(object = list(lin = l.lin, random = l.re, f = l.f), nterms = nterms, b.fixed = res$b.fixed, b.random = res$b.random)

	} else {
		res <- list()
		res$X <- X
		res$Z <- Z
		res$G <- G
		res$y <- y # TODO: excludings NAs
		res$weights <- weights.na # TODO: excludings NAs
		res$family <- family
	}
	res$call<- match.call()
	res$data <- data	
	res$formula <- formula
	res$lin <- if(fit) {
		tmp$lin
	} else {		
		l.lin
	}
	res$random <- if(fit) {
		tmp$random
	} else {		
		l.re
	}
	res$f <- if(fit) {
		tmp$f
	} else {		
		l.f
	}
	res$na.action <- na.action 
	res$names.terms <- names.terms # The same as terms, but ordered: fixed, random and smooth
	res$model.terms <- model.terms # The name of the variables (not the terms), ordered: fixed, random and smooth  
	res$nterms <- nterms
	class(res) <- "sop"
	invisible(res)
}
####################################################################
####################################################################