\name{transition_wacpmix}
\docType{other}
\alias{wacpmix}
\title{Transition: Auto-correlation and Prior Mixture Distribution With Weights}
\description{
Transition probabilities are defined in terms of the mixture of an auto-correlation distribution and a weighted prior distribution (taken from covariate data).}
\usage{acpmix}
\details{
\subsection{Arguments}{
\tabular{lll}{
Name \tab Type \tab Count \cr
------ \tab ------ \tab ------ \cr
\eqn{\alpha} \tab probability \tab 1 \cr
\eqn{\gamma} \tab probability \tab 1 \cr
}
where \eqn{\alpha} is the self-transition probability and \eqn{\gamma} is the mixture weight for the auto-correlation distribution. \eqn{gamma} is always fixed.

Self transitions are defined as \eqn{P(Z = i | Z = i) = \gamma \alpha + (1 - \gamma) p_{self}}.

Outgoing transitions are defined as \eqn{P(Z = j | Z = i) = (\gamma (1 - \alpha) + (1 - \gamma) p_{other}) w_{ij}}.

Prior probabilities, \eqn{p_{self}} and \eqn{p_{other}} taken from (prior) covariate at target position. Weights \eqn{w_{ij}} also taken from covariates

}
\subsection{Options}{
\tabular{lll}{
Name \tab Type \tab Description \cr
------ \tab ------ \tab ------ \cr
maxIters \tab integer (> 0)  \tab maximum number of numeric optimization iterations per EM step\cr
tolerance \tab numeric (> 0) \tab boundary barriers for \eqn{\alpha}\cr
}

\eqn{\alpha} is constrained to the range \eqn{[tolerance, 1 - tolerance]}.
}

\subsection{Covariates}{
\tabular{llll}{
Name \tab Type \tab Dimension \tab Description \cr
------ \tab ------ \tab ------ \tab ------ \cr
priors \tab log probability \tab 2 \tab prior for self and outgoing transitions\cr
weights \tab log probability \tab N \tab relative weights \eqn{w_{ij}} for outgoing transitions\cr
}
where weights \eqn{w_{ij}} are the relative probability of going from state \eqn{i} (the current state) to state \eqn{j != i} (state labeled with \eqn{j} in the \code{valid_transition} argument of \code{new.qhmm}); \eqn{N} is the number of valid transitions from state \eqn{i}.
}

}
\seealso{acpmix, autocorr, autocorr_covar, autocorr_wcovar}
\keyword{transitions}
\keyword{qhmm}
