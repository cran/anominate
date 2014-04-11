\name{anominate}
\alias{anominate}
\title{alpha-NOMINATE: Ideal Point Estimator}

\description{
Fits the Bayesian (MCMC-based) version of Poole and Rosenthal's NOMINATE model described in Carroll, Lewis, Lo, Poole and Rosenthal, ``The Structure of Utility in Spatial Models of Voting,'' \emph{American Journal of Political Science} 57(4): 1008--1028.
The method estimates the structure of political actors' utility functions from roll call data with the inclusion of a separate parameter denoted as alpha.
alpha values of 1 indicate normal (Gaussian) utility, while alpha values of 0 indicate quadratic utility.
}

\usage{
anominate(rcObject, dims=1, nsamp=1000, thin=1, burnin=500, minvotes=20, 
	lop=0.025, polarity=1, random.starts=TRUE, verbose=FALSE, constrain=FALSE)
}

\arguments{
  \item{rcObject}{An roll call matrix of class \code{rollcall}, from Simon Jackman's \code{pscl} package}
  \item{dims}{Number of dimensions to estimate}
  \item{nsamp}{Total number of iterations for the sampler. \code{nsamp} divided by \code{thin} must be larger than \code{burnin}.}
  \item{thin}{Thinning interval}
  \item{burnin}{Number of iterations to be discarded}
  \item{minvotes}{Minimum number of votes required for a legislator to be included in the analysis}
  \item{lop}{A proportion between 0 and 1, the cut-off used for excluding lopsided votes, expressed as the proportion of non-missing votes on the minority side. The default, \code{lop=0.025}, eliminates votes where the minority is smaller than 2.5 percent}
  \item{polarity}{A vector specifying the row number of the legislator(s) constrained to have a positive (i.e., right-wing or conservative) score on each dimension}
  \item{random.starts}{If \code{TRUE}, initial values for the legislator and bill parameters are randomly drawn from a uniform distribution. If \code{FALSE}, the W-NOMINATE estimates are used as the initial values}
  \item{verbose}{If \code{TRUE}, the progress of the sampler at each 100th iteration is printed to the screen}
  \item{constrain}{If \code{FALSE}, this constrains alpha=1 and does not estimate it.}
}

\value{
A list with the following elements:
\item{alpha}{An object of class \code{mcmc} with the sampled values of the alpha parameter}
\item{beta}{An object of class \code{mcmc} with the sampled values of the beta parameter}
\item{legislators}{A object of class \code{mcmc} with the sampled values of the legislator ideal points, with each dimension stored in a separate list (e.g., the first dimension coordinates are stored in \code{legislators[[1]]}, the second dimension coordinates in \code{legislators[[2]]}, etc.)}
\item{yea.locations}{A object of class \code{mcmc} with the sampled values of the Yea locations (midpoint - spread in W-NOMINATE) for each vote, with each dimension stored in a separate list (e.g., the first dimension coordinates are stored in \code{yea.locations[[1]]}, the second dimension coordinates in \code{yea.locations[[2]]}, etc.)}
\item{nay.locations}{A object of class \code{mcmc} with the sampled values of the Nay locations (midpoint + spread in W-NOMINATE) for each vote, with each dimension stored in a separate list (e.g., the first dimension coordinates are stored in \code{nay.locations[[1]]}, the second dimension coordinates in \code{nay.locations[[2]]}, etc.)}
\item{wnom.result}{An object of class \code{nomObject} with the W-NOMINATE results}
}

\references{
Carroll, Royce, Jeffrey B. Lewis, James Lo, Keith T. Poole and Howard Rosenthal. 2013. ``The Structure of Utility in Spatial Models of Voting.'' \emph{American Journal of Political Science} 57(4): 1008--1028.
}


\author{
Christopher Hare, Royce Carroll, Jeffrey B. Lewis, James Lo, Keith T. Poole and Howard Rosenthal
}

\seealso{
  '\link{summary.anominate}','\link{plot.anominate}','\link{densplot.anominate}','\link{traceplot.anominate}'.
}


\examples{

data(sen111)

### This command conducts estimates, which we instead load using data()
#sen111_anom <- anominate(sen111, dims=1, polarity=2, nsamp=200, thin=1,
#	 burnin=100, random.starts=FALSE, verbose=TRUE)

data(sen111_anom)

summary(sen111_anom)

## Graphical summaries
plot(sen111_anom)
densplot.anominate(sen111_anom)
traceplot.anominate(sen111_anom)

}

\keyword{ideal point estimation, NOMINATE, Bayesian latent variable models}
