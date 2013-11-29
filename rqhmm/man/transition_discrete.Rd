\name{transition_discrete}
\docType{other}
\alias{discrete}
\title{Transition: Discrete Distribution}
\description{
Discrete (sometimes called multinomial) distribution with a separate parameter for each transition.}
\usage{discrete}
\details{
\subsection{Arguments}{
\tabular{lll}{
Name \tab Type \tab Count \cr
------ \tab ------ \tab ------ \cr
\eqn{p_{ij}} \tab probability \tab \eqn{N} \cr
}
where \eqn{p_{ij}} is the probability of going from state \eqn{i} (the current state) to state \eqn{j} (state labeled with \eqn{j} in the \code{valid_transition} argument of \code{new.qhmm}; \eqn{N} is the number of valid transitions from state \eqn{i}.
}
\subsection{Options}{
\tabular{lll}{
Name \tab Type \tab Description \cr
------ \tab ------ \tab ------ \cr
pseudo_count \tab numeric (>= 0) \tab transition pseudo count used during EM \cr
}
}

\subsection{Covariates}{
\emph{None.}
}

}
\keyword{transitions}
