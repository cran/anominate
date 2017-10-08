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
Carroll, Royce, Jeffrey B. Lewis, James Lo, Keith T. Poole and Howard Rosenthal. 2013. ``The Structure of Utility in Spatial Models of Voting.'' \emph{American Journal of Political Science} 57(4): 1008--1028.
}


\author{
Christopher Hare, Royce Carroll, Jeffrey B. Lewis, James Lo, Keith T. Poole and Howard Rosenthal
}

\seealso{
  '\link{summary.anominate}','\link{anominate}','\link{densplot.anominate}','\link{traceplot.anominate}'.
}


\examples{
data(sen111_anom)
plot(sen111_anom)
}

\keyword{ideal point estimation, NOMINATE, Bayesian latent variable models}
