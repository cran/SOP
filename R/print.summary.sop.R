####################################################################
####################################################################
print.summary.sop <- function(x, ...) {
	print(x$family)
	cat("\nFormula:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
	
	#cat("\nRandom terms: \n")
	#print(round(x$b.random, 8), quote = FALSE, justify = "right")
	cat("\n\nFixed terms: \n")
	print(round(x$b.fixed, 8), quote = FALSE, justify = "right")
	# class(x) <- c("summary.sop")
	
	cat("\n\nEstimated degrees of freedom:\n")
	print(round(x$edf,4), quote = FALSE, justify = "right")
  
	cat("\n")
	if (!is.null(x$r.sq.adj)) { 
		cat("R-sq.(adj) = ", formatC(x$r.sq, digits = 3, width = 5),"  ")
	}
	if (length(x$dev.expl) > 0) { 
		cat("Deviance explained = ", formatC(x$dev.expl * 100, digits = 3, width = 4, ), "%  ", "n = ", x$n, "\n", sep = "")
	}
	if (!is.null(x$na.action)) {
		cat("(",length(x$na.action)," observations deleted due to missingness)\n", sep = "")
	}
	cat("\n")
	cat("Number of iterations: ", x$iter, sep = "")
	cat("\n")
	invisible(x)
}
####################################################################
####################################################################