####################################################################
####################################################################  
sop.control<-function( maxit = 200, epsilon = 1e-6, trace = FALSE) {
   if (!is.numeric(epsilon) || epsilon <= 0) 
      stop("value of 'epsilon' must be > 0")

   if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")     
   list(maxit = maxit, epsilon = epsilon, trace = trace)
}
####################################################################
####################################################################