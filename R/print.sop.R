####################################################################
####################################################################
print.sop <- function(x, ...) {
  print(x$family)
  cat("Formula:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
	cat("\nEstimated degrees of freedom:\n")
  edf <- x$out$edf
  pord <- 0
  if (length(x$f)>0) {
    for (i in 1:length(x$f)) {
      pord<-pord+sum(x$f[[i]]$pord)
    }
  }
  nam <- names(edf)
	edf <- c(edf, sum(edf), sum(edf) + prod(pord))
  names(edf) <- c(nam, "Total edf", "Total")
	print(round(edf,4), quote = FALSE, justify = "right")
}
####################################################################
####################################################################
