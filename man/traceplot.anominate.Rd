\name{traceplot.anominate}
\alias{traceplot.anominate}
\title{alpha-NOMINATE Traceplot Function}

\description{
  \code{densplot.anominate} reads a \code{anominate} object and plots a traceplot of the alpha parameter using the \code{traceplot} function in \code{coda}.
}

\usage{
  traceplot.anominate(x, ...)
}

\arguments{
  \item{x}{ output of class \code{anominate}, generated as output from function \code{anominate}.}
  \item{...}{ other arguments to \code{traceplot}}
}

\value{
   A traceplot of alpha.
}

\references{
Carroll, Royce, Jeffrey B. Lewis, James Lo, Keith T. Poole and Howard Rosenthal. 2013. ``The Structure of Utility in Spatial Models of Voting.'' \emph{American Journal of Political Science} 57(4): 1008--1028.
}

\author{
Christopher Hare, Royce Carroll, Jeffrey B. Lewis, James Lo, Keith T. Poole and Howard Rosenthal
}

\seealso{
  '\link{summary.anominate}','\link{anominate}','\link{plot.anominate}','\link{traceplot.anominate}'.
}


\examples{
data(sen111)
data(sen111_anom)

traceplot.anominate(sen111_anom)
}

\keyword{ideal point estimation, NOMINATE, Bayesian latent variable models}
