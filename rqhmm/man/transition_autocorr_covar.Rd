\name{transition_autocorr_covar}
\docType{other}
\alias{autocorr_covar}
\title{Transition: Auto-correlation Distribution (Covariate version)}
\description{
Transition probabilities are defined in terms of an auto-correlation parameter (self transition probability) and optionally weights for outgoing transitions.}
\usage{autocorr_covar}
\details{
Requires that the first valid transition be the self-transition.

\subsection{Arguments}{
\emph{None.}
}

\subsection{Options}{
\emph{None.}
}

\subsection{Covariates}{
\tabular{llll}{
Name \tab Type \tab Dimension \tab Description \cr
------ \tab ------ \tab ------ \tab ------ \cr
alpha \tab probability  \tab 1 \tab self-transition probability\cr
}
}

}
\seealso{autocorr, autocorr_wcovar}
\keyword{transitions}
