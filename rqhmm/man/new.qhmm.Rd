\name{new.qhmm}
\alias{new.qhmm}
\title{Create QHMM Instance}
\description{Creates a new HMM instance with default parameters}

\usage{new.qhmm(data.shape, valid.transitions, transition.functions, emission.functions, transition.groups = NULL, emission.groups = NULL, support.missing = FALSE, enable.debug = FALSE)
}

\arguments{
  \item{data.shape}{list with two vectors defining the dimensions of the observation and covariate data tracks, respectively. Covariates are optional, use NULL when no covariates are needed.}
  \item{valid.transitions}{matrix defining the valid transitions and their per state labels. invalid transitions are marked with 0, and valid transition labels are consecutive integers starting from 1.}
  \item{transition.functions}{character vector of transition function names, one per state}
  \item{emission.functions}{list of emission function name vectors; each state has one name per "slot" (observation track)}
  \item{transition.groups}{list of vectors of state numbers, each vector defining a group of transition distributions that share the same parameters (must have the same transition distributions and an identical number of outgoing transitions.}
  \item{emission.groups}{emission group object defining which state/slot sets share the same emission parameteres (must have the same emission distributions)}
  \item{support.missing}{logical value; if TRUE, resulting HMM will accept a missing data matrix (one row per observation track) in the various HMM function calls.}
  \item{enable.debug}{logical value; if TRUE, emission and transition distribution instances are augmented with code that checks for various errors, providing additional debug information at the cost of slower execution.}
}

\details{

\subsection{Data Shape}{

QHMM instances can make use of two distinct classes of data: observations and covariates. Observations are required, while covariates are entirely optional. In either case, each class is formed by a set of "slots" (or tracks), where each one may be uni- or multi-dimensional. When passing data to one of the qhmm instance functions, tracks are concatenated into a single matrix (per class).
}

\subsection{Valid Transitions}{
The \code{valid.transitions} matrix is a key piece in the description of the HMM structure. It is a square matrix in the number of states, where each row defines and labels the valide transitions out of the corresponding state. Positions corresponding to invalid transitions are labeled with 0, and those corresponding to valid transitions are labeled with integers from 1 to N (where N is the count of valid transitions for that state). Labels are defined on a per state basis and play two roles. First, they define an ordering for outgoing transitions which the transition distribution implementations may leverage to match parameters to transitions. Second, when a set of states is part of a transition group, these labels serve to match equivalent transitions across the group for purposes of parameter inference in EM.
}

}

\value{List of class \code{qhmm} with three fields:
\item{n.states}{Number of states in the HMM}
\item{n.emission.slots}{Number of "slots" (independent emission tracks) required by the HMM}
\item{valid.transitions}{Valid transition table for the HMM}

Note that, although these values are part of the definition of the HMM, they should be treated as read-only as any changes to them will be ignored or cause unexpected behaviour.
}

\note{Return value contains a pointer to a C allocated structure stored
  as an attribute (handle_pte). Do not change this value directly!
}

\author{André Luís Martins}

\seealso{new.emission.groups, add.emission.groups}

\examples{
## TO DO
}

\keyword{qhmm}
