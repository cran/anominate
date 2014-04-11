\name{sen111}
\alias{sen111}
\title{111th U.S. Senate Roll Call Vote Matrix}

\description{
  This dataframe contains a matrix of votes cast by U.S. Senators
  in the 111th Congress.  The data are formatted consistent with the
  \code{rollcall} object format in Simon Jackman's \code{pscl} package.}

\usage{
  data(sen111)
}


\value{
  The dataframe contains roll call data for all Senators in the 111th Senate.
  The data is formatted as a \code{rollcall} object with the following elements.

  \item{votes}{ data frame, containing all data from the old \code{nom31.dat} file about
        legislators. For a typical W-NOMINATE object run with an ORD file read using
        \code{readKH}, it will contain the following:
    \itemize{
    \item{\code{state}}{ State name of legislator.}
    \item{\code{icpsrState}}{ ICPSR state code of legislator.}
    \item{\code{cd}}{ Congressional District number.}
    \item{\code{icpsrLegis}}{ ICPSR code of legislator.}
    \item{\code{party}}{ Party of legislator.}
    \item{\code{partyCode}}{ ICPSR party code of legislator.}
    }
    }

  \item{codes}{ list of four vectors. \code{yea} shows the codes in \code{votes} that are yea
        votes, 'nay' shows nay codes, \code{notInLegis} shows absences, and \code{missing}
        shows the missing codes. }
  \item{n}{ numeric, number of legislators }
  \item{m}{ numeric, number of roll calls }
  \item{legis.data}{ data frame, containing the following information on legislators:
    \itemize{
    \item{\code{state}}{ State name of legislator.}
    \item{\code{icpsrState}}{ ICPSR state code of legislator.}
    \item{\code{cd}}{ Congressional District number.}
    \item{\code{icpsrLegis}}{ ICPSR code of legislator.}
    \item{\code{party}}{ Party of legislator.}
    \item{\code{partyCode}}{ ICPSR party code of legislator.}
   }
   }
  \item{vote.data}{ null, would otherwise be a data frame containing data on the votes. }
  \item{desc}{ null, would otherwise be a string describing the data set. }
  \item{source}{ string, describing where data set was read from. }

}

\source{
  Keith Poole. 2014. \emph{111th Senate Roll Call Vote Data}.
  \url{http://www.voteview.com/}.
}

\references{
  Carroll, Royce, Jeffrey B. Lewis, James Lo, Keith T. Poole and Howard Rosenthal. 2013. ``The Structure of
  Utility in Spatial Models of Voting.'' American Journal of Political Science 57(4): 1008--1028.
}


\author{
Christopher Hare, Royce Carroll, Jeffrey B. Lewis, James Lo, Keith T. Poole and Howard Rosenthal
}

\seealso{
  '\link{summary.anominate}','\link{anominate}','\link{densplot.anominate}','\link{traceplot.anominate}','\link{plot.anominate}.
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
