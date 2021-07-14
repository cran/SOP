deviance2 <- function(C, w, sigma2, ssr, edf) {
  log_det_C <- determinant(C)$modulus
  deviance <- log_det_C + sum(log(sigma2*1/w)) + ssr/sigma2 + edf
  deviance  
}
deviance <- function(C, G, w, sigma2, ssr, edf) {
  log_det_C <- determinant(C)$modulus
  log_det_G <- determinant(G)$modulus
  deviance <- log_det_C + log_det_G + sum(log(sigma2*1/w)) + ssr/sigma2 + edf
  deviance  
}
###################################################### 
insert.na <- function(vec, na.ind) {
  aux <- rep(NA, l = length(na.ind))
  aux[!na.ind] <- vec
  aux
}
insert.matrix.na <- function(mat, na.ind) {
  aux <- matrix(NA, ncol = ncol(mat), nrow = length(na.ind))
  aux[!na.ind,] <- mat
  aux
}
add_coeff_to_terms <- function(object, nterms, b.fixed, b.random) {
  nlin <- nterms[1]
  nre <- nterms[2]
  nfun <- nterms[3]

  dim.fixed <- c(if(nlin > 0) object$lin$dim, if(nfun > 0) unlist(lapply(object$f, function(x) sum(x$Xmat$dim$fixed))))
  dim.random <- c(if(nre > 0) object$random$dim, if(nfun > 0) unlist(lapply(object$f, function(x) x$Xmat$dim$random)))

  e.f <- cumsum(dim.fixed)
  s.f <- e.f - dim.fixed + 1

  e.r <- cumsum(dim.random)
  s.r <- e.r - dim.random + 1

  if (nlin > 0) {
    object$lin$coeff <- b.fixed[2:(sum(object$lin$dim)+1)]
  }

  if (nre > 0) {
    object$random$coeff <- b.random[1:sum(object$random$dim)]
  }

  if (nfun > 0) {
    for (i.fun in 1:nfun) {
      if(ncol(object$f[[i.fun]]$Xmat$X) == 0) {
        vcum.f <- NULL
      } else {                      
        vcum.f <- s.f[i.fun + nlin]:e.f[i.fun + nlin]
      }
      vcum.r <- s.r[i.fun + nre]:e.r[i.fun + nre]
    }
    object$f[[i.fun]]$Xmat$coeff.fixed <- b.fixed[vcum.f + 1]
    object$f[[i.fun]]$Xmat$coeff.random <- b.random[vcum.r]
  }
  object
}
######################################################
interpret.f.formula <- function(formula) {
    env <- environment(formula) 
    if(inherits(formula, "character"))  {          
      formula <- as.formula(formula)
    }
    
    tf <- terms.formula(formula, specials = c("f"))
    terms <- attr(tf, "term.labels")
    nt <- length(terms)
    res <- eval(parse(text = terms[1]), envir = env)
    res
  }
######################################################
construct.1D.pspline <- function(formula, data) {
	env <- environment(formula) 
	if(inherits(formula, "character")) {          
		formula <- as.formula(formula)
	}	
  
	res <- interpret.f.formula(formula)
	x1 <- data[ ,res$vars]
	type <- res$type
	
  if(!is.null(res$type) && res$type == "adaptive") {
    MM1 <- MM.basis(x1, min(x1), max(x1), res$nseg[1], res$degree[1], res$pord[1], 5)
    smooth.comp <- paste("ad(", res$vars,")", sep = "")
  } else {
    MM1 <- MM.basis(x1, min(x1), max(x1), res$nseg[1], res$degree[1], res$pord[1], 4)
    smooth.comp <- paste("f(", res$vars,")", sep = "")
  }  
	X <- MM1$X[,-1,drop = FALSE]
  Z <- MM1$Z
  d <- MM1$d
  B <- MM1$B
	c1 <- ncol(B)    

	x.fixed <-  ""
	for(i in 0:(res$pord[1]-1)){
		if(i == 1) 
			x.fixed <- c(x.fixed, res$vars)
		else if( i > 1)
			x.fixed <- c(x.fixed, paste(res$vars ,"^", i, sep = ""))
	}
	names.fixed <- x.fixed	
  names.random <- paste(smooth.comp, c(res$vars), sep = "|")    

	dim.random <- (c1 -res$pord[1])
  dim <- list(fixed = rep(1, ncol(X)), random = sum(dim.random))
	
  names(dim$fixed) <- names.fixed[-1]
	names(dim$random) <- paste(smooth.comp, "Global")
	
  
  v.comp <- 1:ncol(Z)/ncol(Z)
  if (!is.null(res$type) && res$type == "adaptive") {
    C <-  bbase(v.comp, min(v.comp), max(v.comp), res$nseg.sp, res$degree.sp)$B
  } else {
    C <- matrix(MM1$d, ncol = 1)
  }
  G.comp <- t(C)
  g <- vector(mode = "list", length = nrow(G.comp)) 
  for (i in 1:nrow(G.comp)) {
    g[[i]] <- G.comp[i,]
  }

  names(g) <- rep(names.random, length(g))
	colnames(X) <- names.fixed[-1]
	colnames(Z) <- paste(smooth.comp, 1:ncol(Z), sep = ".")
	
	attr(dim$fixed, "random") <- rep(FALSE, length(dim$fixed))
	attr(dim$fixed, "smooth") <- rep(TRUE, length(dim$fixed))
  
	attr(dim$random, "random") <- rep(TRUE, length(dim$random)) 
	attr(dim$random, "smooth") <- rep(TRUE, length(dim$random))

	terms <- list()
	terms$MM <- MM1

	attr(terms, "term") <- smooth.comp    
	init.var <- rep(1, length(g))

  # Center the matrices
  cmX <- colMeans(X)
  cmZ <- colMeans(Z)

  X <- sweep(X, 2, cmX)
  Z <- sweep(Z, 2, cmZ)

	res <- list(X = X, Z = Z, g = g, init.var = init.var, dim = dim, terms = terms, edflabel = names.random, cm = list(X = cmX, Z = cmZ))
	res
}
######################################################
######################################################
construct.2D.pspline <- function(formula, data) {
  env <- environment(formula) 
  if(inherits(formula, "character")) {
    formula <- as.formula(formula)
  }
  
  res <- interpret.f.formula(formula)    
  x1 <- data[ ,res$vars[1]]
  x2 <- data[ ,res$vars[2]]    
  
  type <- res$type

  MM1 <- MM.basis(x1, min(x1), max(x1), res$nseg[1], res$degree[1], res$pord[1], 4)
  MM2 <- MM.basis(x2, min(x2), max(x2), res$nseg[2], res$degree[2], res$pord[2], 4)    

  X1 <- MM1$X; Z1 <- MM1$Z; d1 <- MM1$d; B1 <- MM1$B
  X2 <- MM2$X; Z2 <- MM2$Z; d2 <- MM2$d; B2 <- MM2$B    
  c1 <- ncol(B1); c2 <- ncol(B2)

  # Very dirty: names of the fixed part (null) part of the 2D P-spline:
  # Depend on the penatly order assume for each covariate.
  # For intance, if we assume pord = 2 for x1 (i.e. X1 = [1|x1] and pord = 3 for x2 (X2 = [1|x2|x2^2]), then we have that
  # X = [1 | x1 | x2 | x1:x2| x2^2 | x1:x2^2]  

  x.fixed <- y.fixed <- ""
  for(i in 0:(res$pord[1]-1)){
    if(i == 1) 
      x.fixed <- c(x.fixed, res$vars[1])
    else if( i > 1)
      x.fixed <- c(x.fixed, paste(res$vars[1], "^", i, sep = ""))
  }
  for(i in 0:(res$pord[2]-1)){
    if(i == 1) 
      y.fixed <- c(y.fixed, res$vars[2])
    else if( i > 1)
      y.fixed <- c(y.fixed, paste(res$vars[2], "^", i, sep = ""))
  }
  xy.fixed <- NULL
  for(i in 1:length(y.fixed)) {
    xy.fixed <- c(xy.fixed, paste(y.fixed[i], x.fixed, sep= ""))
  }
  xy.fixed <- xy.fixed[xy.fixed != ""]
  names.fixed <- xy.fixed
  
  smooth.comp <- paste("f(", res$vars[1],",", res$vars[2],")", sep = "")
  names.random <- paste(smooth.comp, c(res$vars[1], res$vars[2]), sep = "|")
  
  X <- Rten2(X2, X1)
  
  # Delete the intercept
  X <- X[,-1,drop = FALSE]
  Z <- cbind(Rten2(X2, Z1), Rten2(Z2, X1), Rten2(Z2, Z1))
  
  dim.random <- c((c1 -res$pord[1])*res$pord[2] , (c2 - res$pord[2])*res$pord[1] , (c1 - res$pord[1])*(c2 - res$pord[2]))  	
  dim <- list(fixed = rep(1, ncol(X)), random = sum(dim.random))
  names(dim$fixed) <- names.fixed
  names(dim$random) <- paste(smooth.comp, "Global")
  
  # Variance/Covariance components: two variances, one for each covariate
  g1u <- rep(1, res$pord[2])%x%d1
  g2u <- d2%x%rep(1, res$pord[1])
  g1b <- rep(1,c2 - res$pord[2])%x%d1
  g2b <- d2%x%rep(1,c1 - res$pord[1])  
  g <- list()	
  g[[1]] <- c(g1u, rep(0, dim.random[2]), g1b)
  g[[2]] <- c(rep(0, dim.random[1]), g2u, g2b)    
  names(g) <- names.random
  
  colnames(X) <- names.fixed
  colnames(Z) <- paste(smooth.comp, 1:ncol(Z), sep = ".")
  
  attr(dim$fixed, "random") <- rep(FALSE, length(dim$fixed))
  attr(dim$fixed, "smooth") <- rep(TRUE, length(dim$fixed))
  
  attr(dim$random, "random") <- rep(TRUE, length(dim$random)) 
  attr(dim$random, "smooth") <- rep(TRUE, length(dim$random))
  
  terms <- list()
  terms$MM <- list(MM1 = MM1, MM2 = MM2)    
  attr(terms, "term") <- smooth.comp
  
  # Initialize variance components
  init.var <- rep(1, length(g))

  # Center the matrices
  cmX <- colMeans(X)
  cmZ <- colMeans(Z)

  X <- sweep(X, 2, cmX)
  Z <- sweep(Z, 2, cmZ)
  
  res <- list(X = X, Z = Z, dim = dim, g = g, init.var = init.var, terms = terms, cm = list(X = cmX, Z = cmZ))	
  res
}
######################################################
######################################################
construct.3D.pspline <- function(formula, data) {
  env <- environment(formula) 
  if(inherits(formula, "character"))  {          
    formula <- as.formula(formula)
  }

  res <- interpret.f.formula(formula)
  x1 <- data[ ,res$vars[1]]
  x2 <- data[ ,res$vars[2]]    
  x3 <- data[ ,res$vars[3]]    

  type <- res$type
  intercept <- TRUE

  MM1 = MM.basis(x1, min(x1), max(x1), res$nseg[1], res$degree[1], res$pord[1], 4 , intercept)
  MM2 = MM.basis(x2, min(x2), max(x2), res$nseg[2], res$degree[2], res$pord[2], 4 , intercept)    
  MM3 = MM.basis(x3, min(x3), max(x3), res$nseg[3], res$degree[3], res$pord[3], 4 , intercept)

  X1 <- MM1$X; Z1 <- MM1$Z; d1 <- MM1$d; B1 <- MM1$B
  X2 <- MM2$X; Z2 <- MM2$Z; d2 <- MM2$d; B2 <- MM2$B    
  X3 <- MM3$X; Z3 <- MM3$Z; d3 <- MM3$d; B3 <- MM3$B

  c1 <- ncol(B1)
  c2 <- ncol(B2)
  c3 = ncol(B3)      
  
  x.fixed <- y.fixed <- z.fixed <- ""

  for(i in 0:(res$pord[1]-1)){
    if(i == 1) 
      x.fixed <- c(x.fixed, res$vars[1])
    else if( i > 1)
      x.fixed <- c(x.fixed, paste(res$vars[1], "^", i, sep = ""))
  }
  for(i in 0:(res$pord[2]-1)){
    if(i == 1) 
      y.fixed <- c(y.fixed, res$vars[2])
    else if( i > 1)
      y.fixed <- c(y.fixed, paste(res$vars[2], "^", i, sep = ""))
  }
  for(i in 0:(res$pord[3]-1)){
    if(i == 1) 
      z.fixed <- c(z.fixed, res$vars[3])
    else if( i > 1)
      z.fixed <- c(z.fixed, paste(res$vars[3], "^", i, sep = ""))
  }	
  xy.fixed <- NULL
  for(i in 1:length(x.fixed)) {
		for(j in 1:length(y.fixed)) {
			xy.fixed <- c(xy.fixed, paste(z.fixed, y.fixed[j], x.fixed[i], sep= ""))
		}
	}
  xy.fixed <- xy.fixed[xy.fixed != ""]	
  names.fixed <- xy.fixed
      
  smooth.comp <- paste("f(", res$vars[1],",", res$vars[2],",", res$vars[3],")", sep = "")
  names.random <- paste(smooth.comp, c(res$vars[1], res$vars[2], res$vars[3]), sep = "|")
  
  rx12 <- Rten2(X1,X2)			
  X <- Rten2(rx12, X3)
  # Delete the intercept
  X <- X[,-1,drop = FALSE]
  colnames(X) <- names.fixed

  Z <- cbind(Rten2(Rten2(Z1,X2), X3), 
             Rten2(Rten2(X1,Z2), X3),
             Rten2(rx12, Z3),
             Rten2(Rten2(Z1,Z2), X3),
             Rten2(Rten2(Z1,X2), Z3),
             Rten2(Rten2(X1,Z2), Z3),
             Rten2(Rten2(Z1,Z2), Z3))
  
  dim.random <- c((c1 -res$pord[1])*res$pord[2] ,
                  (c2 - res$pord[2])*res$pord[1] ,
                  (c1 - res$pord[1])*(c2 - res$pord[2]))		
  dim <- list(fixed = rep(1, ncol(X)), random = sum(dim.random))
  names(dim$fixed)  <- names.fixed
  names(dim$random) <- paste(smooth.comp, "Global")
  np <-  c(prod(res$pord),
           (c1-res$pord[1])*res$pord[2]*res$pord[3],
           (c2-res$pord[2])*res$pord[1]*res$pord[3],
           (c3-res$pord[3])*res$pord[1]*res$pord[2],
           (c1-res$pord[1])*(c2-res$pord[2])*res$pord[3],
           (c1-res$pord[1])*(c3-res$pord[3])*res$pord[2],
           (c2-res$pord[2])*(c3-res$pord[3])*res$pord[1],
           (c1-res$pord[1])*(c2-res$pord[2])*(c3-res$pord[3])) 
  
  # Variance/Covariance components: two variances, one for each covariate
  r1<-rep(1,res$pord[1])
  r2<-rep(1,res$pord[2])
  r3<-rep(1,res$pord[3])
  d1u <- d1%x%r1%x%r2
  d2u <- r1%x%d2%x%r3
  d3u <- r1%x%r2%x%d3
  
  r01<-rep(1,c1-res$pord[1])
  r02<-rep(1,c2-res$pord[2])
  r03<-rep(1,c3-res$pord[3])
  d11b <- d1%x%r02%x%r3
  d12b <- d1%x%r2%x%r03
  
  d21b <- r01%x%d2%x%r3
  d22b <- r1%x%d2%x%r03
  
  d31b <- r01%x%r2%x%d3
  d32b <- r1%x%r02%x%d3
  
  d1t <- d1%x%r02%x%r03
  d2t <- r01%x%d2%x%r03
  d3t <- r01%x%r02%x%d3
  
  # Number of parameters in each part: 8, 44, 44, 44, 242, 242, 242, 1331
  
  D <- diag(c(rep(0,np[1]), rep(1,sum(np[-1]))))
  
  G1inv.n <- c(d1u, rep(0, sum(np[3:4])), d11b, d12b, rep(0, np[7]), d1t)
  G2inv.n <- c(rep(0, np[2]), d2u, rep(0, np[4]), d21b, rep(0, np[6]), d22b, d2t)
  G3inv.n <- c(rep(0, sum(np[2:3])), d3u, rep(0, np[5]), d31b, d32b, d3t)
  
  g <- list(G1inv.n, G2inv.n, G3inv.n)
  names(g) <- names.random
	colnames(X) <- names.fixed
  colnames(Z) <- paste(smooth.comp, 1:ncol(Z), sep = ".")
  
  attr(dim$fixed, "random") <- rep(FALSE, length(dim$fixed))
  attr(dim$fixed, "smooth") <- rep(TRUE, length(dim$fixed))
  
  attr(dim$random, "random") <- rep(TRUE, length(dim$random)) 
  attr(dim$random, "smooth") <- rep(TRUE, length(dim$random))
  
  terms <- list()
  terms$MM <- list(MM1 = MM1, MM2 = MM2, MM3 = MM3)
  attr(terms, "term") <- smooth.comp
  
  # Initialize variance components
  init.var <- rep(1, length(g))

  # Center the matrices
  cmX <- colMeans(X)
  cmZ <- colMeans(Z)

  X <- sweep(X, 2, cmX)
  Z <- sweep(Z, 2, cmZ)
  
  res <- list(X = X, Z = Z, dim = dim, g = g, init.var = init.var, terms = terms, cm = list(X = cmX, Z = cmZ))	
  res
}
######################################################
######################################################
construct.random.part <- function(formula, data) {
  env <- environment(formula) 
  if(inherits(formula, "character"))          
    formula <- as.formula(formula)
  
  mf <- model.frame(formula, data=data, drop.unused.levels = TRUE, na.action = NULL)
  mt <- terms(mf)    
  f.terms <- attr(mt, "term.labels")[attr(mt,"dataClasses") == "factor"]
  Z <- model.matrix(mt, data = mf, contrasts.arg = lapply(mf[,f.terms, drop = FALSE], contrasts, contrasts=FALSE))
  Z[is.na(Z)] <- 0
  
  attr(mt, "contrast") <- attr(Z,"contrast")
  attr(mt, "xlev") <- .getXlevels(mt, mf)
  
  
  dim <- table(attr(Z,"assign"))[-1]
  
  e <- cumsum(dim)
  s <- e - dim + 1
  
  g <- list()
  for(i in 1:length(dim)) {
    g[[i]] <- rep(0,sum(dim))
    g[[i]][s[i]:e[i]] <- 1
  }
  names(g) <- names(dim) <- attr(mt,"term.labels")
  attr(dim, "random") <- rep(TRUE, length(dim)) 
  attr(dim, "smooth") <- rep(FALSE, length(dim))
  
  # Initialize variance components
  init.var <- rep(0.01, length(g))
  
  res <- list(Z = Z[,-1, drop = FALSE], dim = dim, g = g, init.var = init.var, terms = mt)
  res
}
########################################################
########################################################	
construct.fixed.part <- function(formula, data) {
  env <- environment(formula) 
  if(inherits(formula, "character"))          
    formula <- as.formula(formula)

  mf <- model.frame(formula, data, drop.unused.levels = TRUE)
  mt <- terms(mf)   
  X <- model.matrix(mt, mf)

  dim <- table(attr(X,"assign"))[-1]
  names(dim) <- attr(mt, "term.labels")

  attr(mt, "contrast") <- attr(X,"contrast")
  attr(mt, "xlev") <- .getXlevels(mt, mf)

  var.aux <- attr(mt, "term.labels")[attr(mt, "order") == 1]
  dataClasses <- attr(mt, "dataClasses")
  int <- attr(mt, "term.labels")[attr(mt, "order") != 1]
  for (i in int){
    var.aux <- c(var.aux, i)
    if(any(dataClasses[strsplit(i, ":")[[1]]] == "factor")) {
      dataClasses <- c(dataClasses, "factor")     
    } else {
      dataClasses <- c(dataClasses, "numeric")
      names(dataClasses) <- var.aux
    }
  }
  attr(mt, "dataClasses") <- dataClasses

  attr(dim, "random") <-  attr(dim, "smooth") <- rep(FALSE, length(dim)) 	
  res <- list(X = X[,-1, drop = FALSE], dim = dim, terms = mt)
  res	
}
########################################################
########################################################
########################################################	
rae <- function(x) {
  args <- match.call()
  res <- args$x
  res    	
}
######################################################## 
########################################################  
f <- function (..., nseg = 10, pord = 2, degree = 3) {
	vars <- as.list(substitute(list(...)))[-1]
	d <- length(vars)

	if (length(nseg)<d) nseg=rep(nseg,d)
	if (length(pord)<d) pord=rep(pord,d)
	if (length(degree)<d) degree=rep(degree,d)

	term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
	if (term[1] == ".") {
		stop("f(.) not yet supported.")
	}	
	if (d > 1) { 
		for (i in 2:d) {
			term[i] <- deparse(vars[[i]], backtick = TRUE, width.cutoff = 500)
			if (term[i] == ".") { 
				stop("f(.) not yet supported.")
			}
		}
	}
	for (i in 1:d){
		term[i] <- attr(terms(reformulate(term[i])), "term.labels")
	} 
	nseg.new <- round(nseg)
	if (all.equal(nseg.new,nseg) != TRUE) {
		warning("argument nseg of f() should be integer and has been rounded")
	}
	nseg <- nseg.new
	pord.new <- round(pord)
	if (all.equal(pord.new,pord) != TRUE) {
		warning("argument pord of f() should be integer and has been rounded")
	}
	pord <- pord.new
	
	if (length(unique(term)) != d) { 
		stop("Repeated variables as arguments of a smooth are not permitted")
	}	
	full.call <- paste("f(", term[1], sep = "")
	if (d > 1) { 
		for (i in 2:d) {
			full.call <- paste(full.call, ",", term[i], sep = "")
		}
	}	
	label <- paste(full.call, ")", sep = "")
	ret <- list(vars = term, nseg = nseg, pord = pord, degree = degree, dim = d, label = label)
	ret
}
########################################################
########################################################  
ad <- function (..., nseg = 10, pord = 2 , degree = 3, nseg.sp = 5, degree.sp = 3) {
  vars <- as.list(substitute(list(...)))[-1]
  d <- length(vars)
  if (d > 1) {
    stop("Adaptive B-splines for more than one dimension are not yet supported.")
  }  
  if (length(nseg) < d) nseg <- rep(nseg, d)
  if (length(pord) < d) pord <- rep(pord, d)
  if (length(degree) < d) degree <- rep(degree, d)

  if (length(nseg.sp) < d) nseg.sp <- rep(nseg.sp, d)
  if (length(degree.sp) < d) degree.sp <- rep(degree.sp, d)  

  term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  if (term[1] == ".") {
    stop("f(.) not yet supported.")
  } 

  for (i in 1:d){
    term[i] <- attr(terms(reformulate(term[i])), "term.labels")
  }

  nseg.new <- round(nseg)
  if (all.equal(nseg.new, nseg) != TRUE) {
    warning("argument nseg of ad() should be integer and has been rounded")
  }
  nseg <- nseg.new

  nseg.sp.new <- round(nseg.sp)
  if (all.equal(nseg.sp.new, nseg.sp) != TRUE) {
    warning("argument nseg.sp of ad() should be integer and has been rounded")
  }
  nseg.sp <- nseg.sp.new

  pord.new <- round(pord)
  if (all.equal(pord.new, pord) != TRUE) {
    warning("argument pord of ad() should be integer and has been rounded")
  }
  pord <- pord.new

  label <- paste0("ad(", term[1], ")")
  ret <- list(vars = term, nseg = nseg, pord = pord, degree = degree, nseg.sp = nseg.sp, degree.sp = degree.sp, dim = d, label = label, type = "adaptive")
  ret
}
#######################################################################
#######################################################################
tpower <- function(x, t, p) {
  # Function for truncated p-th power function
  return((x - t) ^ p * (x > t))
}

#######################################################################
#######################################################################
spline.bbase<-function (knots, X., bdeg, eps = 1e-05) {
  if(is.null(attributes(knots)$dx))
    dx <- mean(diff(knots))
  else dx <- attributes(knots)$dx
  P <- outer(X., knots, tpower, bdeg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = bdeg + 1)/(gamma(bdeg + 1) * 
                                          dx^bdeg)
  B <- (-1)^(bdeg + 1) * P %*% t(D)
  B[B < eps] = 0
  B
}

##################################################################
##################################################################
bbase <- function(x, xl = min(x), xr = max(x), ndx = 10,   bdeg = 3, eps = 1e-5) {
  # Function for B-spline basis
  dx <- (xr - xl) / ndx
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)  
  attr(knots,"dx")<-dx
  B <- spline.bbase(knots, x, bdeg, eps)  
  res <- list(B = B, knots = knots)
  res 
}
#######################################################################
#######################################################################
MM.basis <- function (x, xl, xr, ndx, bdeg, pord, decom = 1, intercept = TRUE) {
#print("MM.basis")  
  Bb = bbase(x,xl,xr,ndx,bdeg)
  knots <- Bb$knots
  B = Bb$B
  m = ncol(B)
  n = nrow(B)
  D = diff(diag(m), differences=pord)
  P.svd = svd(crossprod(D))
  U.Z = (P.svd$u)[,1:(m-pord)]  # eigenvectors
  d = (P.svd$d)[1:(m-pord)]     # eigenvalues
  Z = B%*%U.Z
  U.X = NULL
  if(decom == 1) {
    U.X = ((P.svd$u)[,-(1:(m-pord))])
    X = B%*%U.X
  } else if (decom == 2){
    X = NULL
    for(i in 0:(pord-1)){
      X = cbind(X,x^i)
    }
  } else if(decom == 3) {
    U.X = NULL
    for(i in 0:(pord-1)){
      U.X = cbind(U.X, knots[-c((1:(bdeg - 1)),(length(knots)- (bdeg - 1) + 1):length(knots))]^i)
    }
    X = B%*%U.X
  } else if(decom == 4) { # Wood's 2013
    X = B%*%((P.svd$u)[,-(1:(m-pord))])
    #id.v <- rep(1, nrow(X))
    #D.temp = X - ((id.v%*%t(id.v))%*%X)/nrow(X)
    D.temp <- sweep(X, 2, colMeans(X))
    Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
    X <- X%*%Xf
    U.X = ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf
  } else if(decom == 5) { # Paul's parameterization
    U.Z <- (t(D)%*%solve(D%*%t(D)))
    Z <- B%*%U.Z
    X = B%*%((P.svd$u)[,-(1:(m-pord))])
    #id.v <- rep(1, nrow(X))
    #D.temp = X - ((id.v%*%t(id.v))%*%X)/nrow(X)
    D.temp <- sweep(X, 2, colMeans(X))
    Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
    X <- X%*%Xf
    U.X = ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf					
  } else if(decom == 6) { # martin's approach
    U.Z <- t(D)
    Z <- B%*%U.Z  
    X = B%*%((P.svd$u)[,-(1:(m-pord))])
    #id.v <- rep(1, nrow(X))
    #D.temp = X - ((id.v%*%t(id.v))%*%X)/nrow(X)
    D.temp <- sweep(X, 2, colMeans(X))
    Xf <- svd(crossprod(D.temp))$u[,ncol(D.temp):1]
    X <- X%*%Xf
    U.X = ((P.svd$u)[,-(1:(m-pord)), drop = FALSE])%*%Xf			
  }
  if (!intercept) X <- X[,-1,drop = FALSE]
  list(X = X, Z = Z, d = d, B = B, m = m, D = D, knots = knots, U.X = U.X, U.Z = U.Z)
}
#######################################################################
#######################################################################
construct.block <- function(A1,A2,A3,A4) {
  block <- rbind(cbind(A1,A2), cbind(A3,A4))
  return(block)
}
construct.capital.lambda <- function(g) {
  length.eq <- all(sapply(g, function(x) {
    diff(range(unlist(lapply(x, length)))) < .Machine$double.eps ^ 0.5
  }))
  if(length.eq) {
    l <- length(g)
    if(l == 1) {
      if(length(g[[1]]) == 1) {
        res <- g
      } else {
        res <- do.call("c", lapply(g, function(x) x))
      }
    } else {
      dim <- sapply(g, function(x) {
        if(is.list(x))
          unlist(lapply(x, length))[1]
        else
          length(x)
      })		
      end <- cumsum(dim)
      init <- end - dim + 1
      
      res <- do.call("c", lapply(1:length(g), function(x, g, init, end, dim) {
        temp <- g[[x]]
        if(is.list(temp)) {
          lapply(temp, function(y, x, dim) {
            aux <- rep(0, l = sum(dim))
            aux[init[x]:end[x]] <- y
            aux
          }, x = x, dim = dim)
        } else {
          aux <- rep(0, l = sum(dim))
          aux[init[x]:end[x]] <- temp
          list(aux)
        }
      }, g = g, init = init, end = end, dim = dim))	
    }
  } else {
    stop("Error in construct.capital.lambda")
  }	
  res
}
###################################################################
###################################################################
# Model matrices and GLAM
###################################################################
###################################################################
Rten <- function(X) {
  one <- matrix(1, 1, ncol(X))
  kronecker(X,one)*kronecker(one,X)
}
###################################################################
###################################################################
Rten2 <- function(X1,X2) {
  one.1 <- matrix(1,1,ncol(X1))
  one.2 <- matrix(1,1,ncol(X2))
  kronecker(X1,one.2)*kronecker(one.1,X2)
}
###################################################################
###################################################################
H <- function(X,A) {
  d <- dim(A)
  M <- matrix(A, nrow = d[1])
  XM <- X%*%M
  array(XM, c(nrow(XM),d[-1]))
}
###################################################################
###################################################################
Rotate <- function(A) {
  d <- 1:length(dim(A))
  d1 <- c(d[-1],d[1])
  aperm(A, d1)
}
###################################################################
###################################################################
RH <- function(X,A) {
  Rotate(H(X,A))
}
###################################################################
###################################################################
A1.form <- function(l, w = NULL){
  d <- length(l)
  n <- rev(sapply(l,nrow))
  c <- rev(sapply(l,ncol))
  if (is.null(w)) {
    W <- array(1, n)
  } else {
    W <- array(w, n)
  }
  tmp <- RH(t(Rten(l[[d]])), W)
  for (i in (d-1):1) {
    tmp <- RH(t(Rten(l[[i]])),tmp)
  }
  dim(tmp)<- rep(c, rep(2,d))
  Fast1 <- aperm(tmp, as.vector(matrix(1:(d*2), byrow = TRUE, ncol = 2)))
  Fast <- if(prod(c)) matrix(Fast1, nrow = prod(c)) else Fast1
  return(Fast)
}
###################################################################
###################################################################
A2.form <- function(l1, l2, w = NULL) {
  d1 <- length(l1)
  d2 <- length(l2)
  if(!(d1 == d2)) {
    stop("l1 and l2 should have the same dimension")
  }
  n <- rev(sapply(l1, nrow))
  d <- rev(sapply(l1, ncol))
  c <- rev(sapply(l2, ncol))
  
  if (is.null(w)) {
    W <- array(1, n)
  } else {
    W <- array(w, n)
  }
  tmp <- RH(t(Rten2(l2[[d1]], l1[[d1]])), W)
  for (i in (d1-1):1) {
    tmp <- RH(t(Rten2(l2[[i]], l1[[i]])),tmp)
  }
  dim(tmp)<- as.vector(rbind(d,c))
  Fast1 <- aperm(tmp, as.vector(matrix(1:(d1*2), byrow = TRUE, ncol = 2)))
  Fast <- if(prod(d)) matrix(Fast1, nrow = prod(d)) else Fast1
  return(Fast)
}
###################################################################
###################################################################
XtX <- function(X, w = NULL) {
  A1.form(X, w)
}
###################################################################
###################################################################
XtZ <- function(X, Z, w = NULL) {
  d <- length(Z)
  res <- NULL
  for (i in 1:d) {
    res <- cbind(res, A2.form(X, Z[[i]], w))
  }
  res
}
###################################################################
###################################################################
ZtZ <- function(Z, w = NULL) {
  d <- length(Z)
  upper <- list()
  for(i in 1:d) {
    upper[[i]] <- list()
    upper[[i]][[i]] <- A1.form(Z[[i]], w)
  }
  # Obtain the elements of the matrix
  for(i in 1:(d-1)) {
    for(j in (i+1):d) {
      upper[[i]][[j]] <- A2.form(Z[[i]], Z[[j]], w)
    }
  }
  # Create the matrix
  res <- NULL
  for (i in 1:d) {
    if( i == 1) {
      res <- do.call("cbind", upper[[1]])
    } else {
      tmp <- do.call("cbind", upper[[i]])
      for(j in (i-1):1) {
        if(length(upper[[j]][[i]]))
          tmp <- cbind(t(upper[[j]][[i]]), tmp)
      }
      if(nrow(tmp))
        res <- rbind(res, tmp)
    }
  }
  res
}
###################################################################
###################################################################
Xty <- function(X, y, w = NULL) {
  d <- length(X)
  n <- rev(sapply(X, nrow))
  if(is.null(w)) {
    Y <- array(y, n)
  } else {
    Y <- array(w*y, n)
  }
  tmp <- RH(t(X[[d]]),Y)
  for(i in (d-1):1) {
    tmp <- RH(t(X[[i]]), tmp)
  }
  as.vector(tmp)
}
###################################################################
###################################################################
Zty <- function(Z, y, w = NULL) {
  d <- length(Z)
  n <- rev(sapply(Z[[1]], nrow))
  if(is.null(w)) {
    Y <- array(y, n)
  } else {
    Y <- array(w*y, n)
  }
  res <- NULL
  for(i in 1:d) {
    k <- length(Z[[i]])
    tmp <- RH(t(Z[[i]][[k]]),Y)
    for(j in (k-1):1) {
      tmp <- RH(t(Z[[i]][[j]]), tmp)
    }
    res <- c(res, as.vector(tmp))
  }
  res
}
###################################################################
###################################################################
Xtheta <- function(X, theta) {
  d <- length(X)
  n <- rev(sapply(X, ncol))
  Theta <- array(theta, n)
  tmp <- RH(X[[d]], Theta)
  for(i in (d-1):1) {
    tmp <- RH(X[[i]], tmp)
  }
  as.vector(tmp)
}
###################################################################
###################################################################
Ztheta <- function(Z, theta, np) {
  d <- length(Z)
  for(i in 1:d) {
    if (i == 1) {
      res <- Xtheta(Z[[i]], theta[1:(np[1])])
    } else {
      init <- sum(np[1:(i-1)])
      fin  <- np[i]
      if(fin) res <- res + Xtheta(Z[[i]], theta[(init+1):(init+fin)])
    }
  }
  res
}
###################################################################
###################################################################
construct.matrices <- function(X, Z, z, w, GLAM) {
  if(GLAM) {
    XtX. <- XtX(X,w) 
    XtZ. <- XtZ(X,Z,w)
    ZtX. <- t(XtZ.)
    ZtZ. <- ZtZ(Z,w)
    Zty. = Zty(Z,z,w)
    Xty. = Xty(X,z,w)
    yty. <- sum((z^2)*w)
    ZtXtZ = rbind(XtZ., ZtZ.)
    u <- c(Xty.,Zty.)
  } else {
    XtW. = t(X*w)
    XtX. = XtW.%*%X
    XtZ. = XtW.%*%Z
    ZtX. = t(XtZ.)
    ZtW. =  t(Z*w)
    ZtZ. = ZtW.%*%Z
    Xty. = XtW.%*%z
    Zty. = ZtW.%*%z
    yty. <- sum((z^2)*w)
    ZtXtZ = rbind(XtZ., ZtZ.)
    u <- c(Xty.,Zty.)
  }
  res <- list(XtX. = XtX., XtZ. = XtZ., ZtX. = ZtX., ZtZ. = ZtZ., Xty. = Xty., Zty. = Zty., yty. = yty., ZtXtZ = ZtXtZ, u = u)
}
###################################################################
