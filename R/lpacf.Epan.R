lpacf.Epan= function (x, binwidth, lag.max = NULL, filter.number = 10, family = "DaubLeAsymm", 
          smooth.dev = var, AutoReflect = TRUE, tol = 0.1, maxits = 5, 
          ABBverbose = 0, lapplyfn = lapply, allpoints = FALSE, verbose = FALSE, 
          ...) 
{
  TT <- length(x)
  filter <- filter.select(filter.number, family)
  Jp <- ceiling(logb(TT, 2))
  add <- 2^Jp - TT
  Nh <- length(filter$H != 0)
  xa <- c(rep(0, times = add), x)
  lxa <- length(xa)
  dsname = deparse(substitute(x))
  AutoBinWidth <- FALSE
  if (missing(binwidth) || binwidth == 0) {
    binwidth <- lpacf::AutoBestBW.Epan(x = xa, filter.number = filter.number, 
                                family = family, smooth.dev = smooth.dev, AutoReflect = AutoReflect, 
                                tol = tol, maxits = maxits, plot.it = FALSE, verbose = ABBverbose)
    AutoBinWidth <- TRUE
  }
  if (binwidth >= TT) {
    binwidth = lpacf::AutoBestBW.Epan(x = x[(TT - 2^{
      Jp - 1
    } + 1):TT], filter.number = filter.number, family = family, 
    smooth.dev = smooth.dev, AutoReflect = AutoReflect, 
    ...)
  }
  n <- length(x)
  if (is.null(lag.max)) 
    lag.max <- floor(10 * (log10(n)))
  if(binwidth<=lag.max){ # RK add
    lag.max=binwidth-1
    warning(paste("lag.max is larger than binwidth, setting lag.max to",binwidth-1))
  }
  start <- 1:(n - binwidth + 1)
  end <- binwidth:n
  the.x <- (start + end)/2
  if (allpoints == TRUE) {
    vacc.lo <- the.x[1]
    vacc.hi <- the.x[length(the.x)]
    end1=NULL
    start1=NULL
    if((binwidth-lag.max)>1){
      end1=(lag.max + 1):(binwidth-1) # RK: modified from global (binwidth-1)
      start1 <- rep(1, length(end1)) # RK modified from global rep(1,length(end1))
    }
    if (verbose == TRUE) {
      message("Length start1 is: ", length(start1))
      message("Length end1 is: ", length(end1))
    }
    the.x1 <- (start1 + end1)/2
    start2=NULL
    end2=NULL
    if((binwidth-lag.max)>1){
      start2 <- (n - binwidth + 2):(n - lag.max) # RK: modified from global (n-binwidth+2):(n-lag.max-1)
      end2 <- rep(n, length(start2)) # RK: modified from global rep(n, length(start2))
    }
    if (verbose == TRUE) {
      message("Length start2 is: ", length(start2))
      message("Length end2 is: ", length(end2))
    }
    the.x2 <- (start2 + end2)/2
    start <- c(start1, start, start2)
    end <- c(end1, end, end2)
    if (verbose == TRUE) {
      message("Length start is: ", length(start))
      message("Length end is: ", length(end))
    }
    the.x <- c(the.x1, the.x, the.x2)
  }
  se.df <- data.frame(t(cbind(start, end)))
  dvals <- function(sv, x) return(x[seq(from = sv[1], to = sv[2])])
  dv.rep <- lapplyfn(se.df, dvals, x = x)
  mypacf <- function(x, lag.max, plot) return(pacf(x = Epanechnikov(x)*x, lag.max = lag.max, 
                                                   plot = plot)$acf[, , 1])
  dv.pacf <- lapplyfn(dv.rep, mypacf, plot = FALSE, lag.max = lag.max)
  dv.pacf <- matrix(as.vector(unlist(dv.pacf)), ncol = lag.max, 
                    byrow = TRUE)
  if (allpoints == TRUE) {
    out <- list(the.x = the.x, lpacf = dv.pacf, the.x1 = the.x1, 
                the.x2 = the.x2, vacc = c(vacc.lo, vacc.hi), the.vacc = TRUE, 
                binwidth = binwidth, AutoBinWidth = AutoBinWidth)
    class(out) <- "lpacf"
    return(out)
  }
  else {
    out = list(the.x = the.x, lpacf = dv.pacf, the.vacc = FALSE, 
               binwidth = binwidth, AutoBinWidth = AutoBinWidth)
    class(out) = "lpacf"
    return(out)
  }
}
