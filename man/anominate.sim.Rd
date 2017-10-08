\name{anominate.sim}
\alias{anominate.sim}
\title{alpha-NOMINATE Estimation on Simulated Roll Call Matrices Generated with Normal (Gaussian) or Quadratic Utility}

\description{
Generates a matrix of roll call votes based on the assumption that legislators possess either normal (Gaussian) or quadratic utility functions. 
The roll call votes are then analyzed using the ideal point model described in Carroll, Lewis, Lo, Poole and Rosenthal, ``The Structure of Utility in Spatial Models of Voting,'' \emph{American Journal of Political Science} 57(4): 1008--1028.
The estimated value of the alpha parameter can be compared to the true value (1 for normal (Gaussian) utility, 0 for quadratic utility).
}

\usage{
anominate.sim(nvotes=500, nlegis=101, seed=123345, utility="normal")
}

\arguments{
  \item{nvotes}{Number of roll call votes to be simulated}
  \item{nlegis}{Number of legislators to be simulated}
  \item{seed}{Seed for the random number generator}
  \item{utility}{The utility function used to generate the roll call data (``normal'' or ``quadratic'')}
}

\value{
  An object of class \code{rollcall}, a list with the following components:
  
  \item{votes}{a \code{\link{matrix}} containing voting decisions, with
    rows corresponding to legislators (test subjects) and columns to
    roll call votes (test items). Legislators (test subjects)
    and items (or votes) have been labeled in the
    \code{\link{dimnames}} attribute of this matrix, using
    the \code{legis.names} and/or \code{vote.names} arguments, respectively.}

  \item{codes}{a \code{\link{list}} with named components \code{yea},
    \code{nay}, \code{notInLegis} and \code{missing}, each component a
    numeric vector (possibly of length 1 and possibly \code{NA}),
    indicating how entries in the \code{votes} component of the
    \code{rollcall} object should be considered.  This list simply gathers
    up the values in the \code{yea}, \code{nay}, \code{notInLegis} and
    \code{missing} arguments passed to the function.}
  
  \item{n}{numeric, number of legislators, equal to \code{dim(votes)[1]}}

  \item{m}{numeric, number of votes, equal to \code{dim(votes)[2]}}

  \item{legis.data}{user-supplied data on legislators/test-subjects.}
  
  \item{vote.data}{user-supplied data on rollcall votes/test-items.}
  
  \item{desc}{any user-supplied description.  If no description was provided,
    defaults \code{desc} defaults to \code{NULL}.}
  
  \item{source}{any user-supplied source information (e.g., a url or a
    short-form reference).  If no description is provided, \code{source}
    defaults to \code{NULL}.}
  
}

\seealso{
  \code{\link{rollcall}} for the full documentation of a roll call object from Simon Jackman's \code{pscl} package.

  Output from this function is intended for use with \code{anominate}.
}

\references{
Carroll, Royce, Jeffrey B. Lewis, James Lo, Keith T. Poole and Howard Rosenthal. 2013. ``The Structure of Utility in Spatial Models of Voting.'' \emph{American Journal of Political Science} 57(4): 1008--1028.
}


\author{
Christopher Hare, Royce Carroll, Jeffrey B. Lewis, James Lo, Keith T. Poole and Howard Rosenthal
}



\examples{

\dontshow{
quadratic.data <- anominate.sim(utility="quadratic") 
data(quad_anom)
summary(quad_anom)

normal.data <- anominate.sim(utility="normal") 
data(norm_anom)
summary(norm_anom)
}

\dontrun{
quadratic.data <- anominate.sim(utility="quadratic") 
quad_anom <- anominate(quadratic.data, dims=1, polarity=2, nsamp=200, thin=1,
	burnin=100, random.starts=FALSE, verbose=TRUE)
summary(quad_anom)

normal.data <- anominate.sim(utility="normal") 
norm_anom <- anominate(normal.data, dims=1, polarity=2, nsamp=200, thin=1, 
	burnin=100, random.starts=FALSE, verbose=TRUE)
summary(norm_anom)
}

}


\keyword{ideal point estimation, NOMINATE, Bayesian latent variable models}
