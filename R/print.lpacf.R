print.lpacf <-
function (x, ...) 
{
    cat("Class 'lpacf' : Localized Partial Autocorrelation Object:\n")
    cat("       ~~~~  : List with", length(x), "components with names\n")
    cat("             ", names(x), "\n\n")
    cat("\nsummary(.):\n----------\n")
    lpacf::summary.lpacf(x)
}
