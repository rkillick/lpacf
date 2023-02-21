\name{tvar2sim}
\alias{tvar2sim}
\title{
	Simulate a realization from a particular TVAR(2) model.
}
\description{Simulates a realization from a TVAR(2) model where
	both parameters move from -1.1 to 0.5 in equal steps
	over 512 time points. The realization is of length 512.
	The innovations are normally distributed with mean zero and
	standard deviation of \code{sd}. 
}
\usage{
tvar2sim(sd = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sd}{This is the standard deviation of the Gaussian innovation.
}
}
\details{This function is easily converted into one that does the
	same thing but for a different sample size.
}
\value{A realization of the aforementioned TVAR(2) process.
}

\author{
Guy Nason.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{lpacf}}
}
\examples{
#
# Generate realization from the TVAR(2) process
#
x <- tvar2sim()
#
# Maybe plot it
#
\dontrun{ts.plot(x)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}
