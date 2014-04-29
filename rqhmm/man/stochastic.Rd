\name{stochastic.Rd}
\alias{stochastic.backtrace.qhmm}

\title{Stochastic backtrace}
\description{Function generates a sample path through the HMM conditional on parameters and observations.}

\usage{

stochastic.backtrace.qhmm(hmm, emissions, covars = NULL, missing = NULL, fwdmatrix = NULL)

}

\arguments{
  \item{hmm}{QHMM instance object}
  \item{emissions}{numeric vector or matrix with observed sequence.}
  \item{covars}{numeric vector or matrix with covariate values for each position in the sequence.}
  \item{missing}{integer vector or matrix with missing data indicator for each position (0 for present, 1 for missing).}
  \item{fwdmatrix}{forward probability matrix, obtained from \code{forward.qhmm}. If not supplied, it will be computed internally.}
}

\details{
Sample path is computed using the 'forward' algorithm results. Starting with the last position in the sequence, a state is sampled proportionally to the 'forward' state probabilities at that position. Remaining positions (iterating from end to start) are sampled proportionally to the product of the 'forward' state probability at the current position with the transition probabiltiy from that state to the previously sampled state (next position in the sequence).
}

\value{
Integer vector with sampled state numbers at each position.
}


\author{André Luís Martins}

\seealso{forward.qhmm}

\keyword{qhmm}
\keyword{sample}
