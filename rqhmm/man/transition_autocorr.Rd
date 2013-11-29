\name{transition_autocorr}
\docType{other}
\alias{autocorr}
\title{Transition: Auto-correlation Distribution}
\description{
Transition probabilities are defined in terms of an auto-correlation parameter (self transition probability) and optionally weights for outgoing transitions.}
\usage{autocorr}
\details{
Requires that the first valid transition be the self-transition.

\subsection{Arguments}{
\tabular{lll}{
Name \tab Type \tab Count \cr
------ \tab ------ \tab ------ \cr
\eqn{\alpha} \tab probability \tab 1 \cr
\eqn{w_{ij}} \tab probability \tab N-1 (optional)\cr
}
where \eqn{\alpha} is the self-transition probability. \eqn{w_{ij}} is the weight of the transition from state \eqn{i} (the current state) to state \eqn{j \neq i} (state labeled with \eqn{j} in the \code{valid_transition} argument of \code{new.qhmm}); \eqn{N} is the number of valid transitions from state \eqn{i}.

Weights, if present, are fixed, and should sum to one. They must either all be present or none. When weights are present, outgoing transitions are defined as \eqn{P(Z = j | Z = i) = w_{ij} * (1 - \alpha)}. Otherwise, \eqn{P(Z = j | Z = i) = (1 - \alpha)}.
}
\subsection{Options}{
\tabular{lll}{
Name \tab Type \tab Description \cr
------ \tab ------ \tab ------ \cr
pseudo_count \tab numeric (>= 0)  \tab transition pseudo count used during EM \cr
}
}

\subsection{Covariates}{
\emph{None.}
}

}
\seealso{autocorr_covar, autocorr_wcovar}
\keyword{transitions}
\keyword{qhmm}
