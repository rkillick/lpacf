summary.lpacf <-
function (object, ...) 
{
    nlags <- dim(object$lpacf)
    ntime <- nlags[1]
    nlags <- nlags[2]
    cat("Number of times: ", ntime, "\n")
    cat("Number of lags: ", nlags, "\n")
    cat("Range of times from: ", min(object$the.x), " to ", max(object$the.x), "\n")
    if (object$the.vacc==TRUE)
		cat("Complete time series was included (alltimes=TRUE)\n")
    else
		cat("Part series was analyzed (alltimes=FALSE)\n")
    cat("Smoothing binwidth used was: ", object$binwidth, "\n") 
    if (object$AutoBinWidth==TRUE)
	cat("\tBinwidth was chosen automatically\n")
    else
	cat("\tBinwidth was specified by the user\n")
}
