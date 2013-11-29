\name{transition_autocorr_wcovar}
\docType{other}
\alias{autocorr_wcovar}
\title{Transition: Auto-correlation Distribution (Covariate weights version)}
\description{
Transition probabilities are defined in terms of an auto-correlation parameter (self transition probability) and per position weights for outgoing transitions.}
\usage{autocorr_wcovar}
\details{
Requires that the first valid transition be the self-transition.

\subsection{Arguments}{
\tabular{lll}{
Name \tab Type \tab Count \cr
------ \tab ------ \tab ------ \cr
\eqn{\alpha} \tab probability \tab 1 \cr
}
where \eqn{\alpha} is the self-transition probability.
}
\subsection{Options}{
\tabular{lll}{
Name \tab Type \tab Description \cr
------ \tab ------ \tab ------ \cr
pseudo_count \tab numeric (>= 0)  \tab transition pseudo count used during EM \cr
}
}

\subsection{Covariates}{
\tabular{llll}{
Name \tab Type \tab Dimension \tab Description \cr
------ \tab ------ \tab ------ \tab ------ \cr
\eqn{w_ij} \tab log probability \tab N \tab weight of transition from \eqn{i} to \eqn{j}\cr
}
where \eqn{w_ij} is the weight of the transition from state \eqn{i} (the current state) to state \eqn{j \neq i} (state labeled with \eqn{j} in the \code{valid_transition} argument of \code{new.qhmm}); \eqn{N} is the number of valid transitions from state \eqn{i}. 
}

}
\seealso{autocorr, autocorr_covar}
\keyword{transitions}
\keyword{qhmm}
