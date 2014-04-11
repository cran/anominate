\name{quad_anom}
\alias{quad_anom}
\title{alpha-NOMINATE Estimate, Simulated Roll Call Matrix using Quadratic Utility}

\description{
alpha-NOMINATE estimates from simulated roll call matrix using quadratic utility. Estimates in one dimension.
}

\usage{
data(quad_anom)
}

\value{
  An object of class \code{anominate}, a list with the following elements:
\item{alpha}{An object of class \code{mcmc} with the sampled values of the alpha parameter}
\item{beta}{An object of class \code{mcmc} with the sampled values of the beta parameter}
\item{legislators}{A object of class \code{mcmc} with the sampled values of the legislator ideal points, with each dimension stored in a separate list (e.g., the first dimension coordinates are stored in \code{legislators[[1]]}, the second dimension coordinates in \code{legislators[[2]]}, etc.)}
\item{yea.locations}{A object of class \code{mcmc} with the sampled values of the Yea locations (midpoint - spread in W-NOMINATE) for each vote, with each dimension stored in a separate list (e.g., the first dimension coordinates are stored in \code{yea.locations[[1]]}, the second dimension coordinates in \code{yea.locations[[2]]}, etc.)}
\item{nay.locations}{A object of class \code{mcmc} with the sampled values of the Nay locations (midpoint + spread in W-NOMINATE) for each vote, with each dimension stored in a separate list (e.g., the first dimension coordinates are stored in \code{nay.locations[[1]]}, the second dimension coordinates in \code{nay.locations[[2]]}, etc.)}
\item{wnom.result}{An object of class \code{nomObject} with the W-NOMINATE results}
}


\seealso{
  '\link{anominate.sim}','link{anominate}','\link{summary.anominate}','\link{plot.anominate}','\link{densplot.anominate}','\link{traceplot.anominate}'.
}

\references{
Carroll, Royce, Jeffrey B. Lewis, James Lo, Keith T. Poole and Howard Rosenthal. 2013. ``The Structure of Utility in Spatial Models of Voting.'' \emph{American Journal of Political Science} 57(4): 1008--1028.
}


\author{
Christopher Hare, Royce Carroll, Jeffrey B. Lewis, James Lo, Keith T. Poole and Howard Rosenthal
}

\examples{

quadratic.data <- anominate.sim(utility="quadratic") 
### This command conducts estimates, which we instead load using data()
#quad_anom <- anominate(quadratic.data, dims=1, polarity=2, nsamp=200, thin=1, 
#	burnin=100, random.starts=FALSE, verbose=TRUE)
data(quad_anom)
summary(quad_anom)

}

\keyword{ideal point estimation, NOMINATE, Bayesian latent variable models}
