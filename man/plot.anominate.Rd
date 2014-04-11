\name{plot.anominate}
\alias{plot.anominate}
\title{alpha-NOMINATE Coordinate Plot Function}

\description{
  \code{plot.anominate} reads a \code{anominate} object and plots the estimated ideal points.
}

\usage{
  \method{plot}{anominate}(x, ...)
}

\arguments{
  \item{x}{ output of class \code{anominate}, generated as output from function \code{anominate}.}
  \item{...}{ other arguments to \code{plot}}
}

\value{
   A coordinate plot, with 95\% credible interval bars.
}

\references{
  Carroll, Royce, Jeffrey B. Lewis, James Lo, Keith T. Poole and Howard Rosenthal. 2013. ``The Structure of
  Utility in Spatial Models of Voting.'' American Journal of Political Science 57(4): 1008--1028.
}


\author{
Christopher Hare, Royce Carroll, Jeffrey B. Lewis, James Lo, Keith T. Poole and Howard Rosenthal
}

\seealso{
  '\link{summary.anominate}','\link{anominate}','\link{densplot.anominate}','\link{traceplot.anominate}'.
}


\examples{

data(sen111)

### This command conducts estimates, which we instead load using data()
#sen111_anom <- anominate(sen111, dims=1, polarity=2, nsamp=200, thin=1, 
#	burnin=100, random.starts=FALSE, verbose=TRUE)

data(sen111_anom)

summary(sen111_anom)

## Graphical summaries
plot(sen111_anom)
densplot.anominate(sen111_anom)
traceplot.anominate(sen111_anom)

}

\keyword{ideal point estimation, NOMINATE, Bayesian latent variable models}
