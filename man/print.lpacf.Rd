\name{print.lpacf}
\alias{print.lpacf}
\title{Prints a \code{lpacf} object.
}
\description{Prints a \code{lpacf} object, basically telling you what's
	there.
}
\usage{
\method{print}{lpacf}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The \code{lpacf} object to print.
}
  \item{\dots}{Other arguments
}
}
\details{Prints a \code{lpacf} object.
}
\value{None
}
\references{
Killick, R., Knight, M. I., Nason, G.P. and Eckley, I. A. (2020) The local partial autocorrelation function and some applications. \emph{Electron. J. Statist.} \bold{14} (2), 3268-3314. DOI: \href{https://www.doi.org/}{10.1214/20-EJS1748}.
}
\author{
Guy Nason
}

\seealso{\code{\link{lpacf}}, \code{\link{summary.lpacf}}
}
\examples{
#
# Simulate an example
#
x.test <- tvar2sim()
#
# Compute the lpacf
#
x.lpacf <- lpacf(x.test)
#
# Print it out - note, can normally just type name of object
#
print(x.lpacf)
#Class 'lpacf' : Localized Partial Autocorrelation Object:
#       ~~~~  : List with 5 components with names
#              the.x lpacf the.vacc binwidth AutoBinWidth 
#
#
#summary(.):
#----------
#Number of times:  220 
#Number of lags:  27 
#Range of times from:  147  to  366 
#Part series was analyzed (alltimes=FALSE)
#Smoothing binwidth used was:  293 
#        Binwidth was chosen automatically
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
