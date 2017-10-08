\name{summary.anominate}
\alias{summary.anominate}
\title{alpha-NOMINATE Summary Function}

\description{
  \code{summary.anominate} reads a \code{anominate} object and prints a summary.
}

\usage{
  \method{summary}{anominate}(object, ...)
}

\arguments{
  \item{object}{ output of class \code{anominate}, generated as output from function \code{anominate}.}
  \item{...}{ other arguments to \code{print}}
}

\value{
   A summary of the estimate, with number of legislators, votes, dimensions, and summary of alpha estimate.
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
data(sen111)
data(sen111_anom)

summary(sen111_anom)
}

\keyword{ideal point estimation, NOMINATE, Bayesian latent variable models}
