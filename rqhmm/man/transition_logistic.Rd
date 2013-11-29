\name{transition_logistic}
\docType{other}
\alias{logistic}
\title{Transition: Logistic Function}
\description{
Transition probabilities are defined in terms of a set of logistic functions of a single covariate}
\usage{logistic}
\details{
Requires that the first valid transition be the self-transition.

\subsection{Arguments}{
\tabular{lll}{
Name \tab Type \tab Count \cr
------ \tab ------ \tab ------ \cr
\eqn{\beta} \tab numeric \tab 2*(N-1) \cr
}
where \eqn{\beta} is the vector of logistic function parameters, two for each valid outgoing transition from state \eqn{i}.
}
\subsection{Options}{
\tabular{lll}{
Name \tab Type \tab Description \cr
------ \tab ------ \tab ------ \cr
maxIters \tab integer (> 0)  \tab maximum number of numeric optimization iterations per EM step\cr
tolerance \tab numeric (> 0) \tab tolerance used in numeric optimization per EM step\cr
}
}

\subsection{Covariates}{
\tabular{llll}{
Name \tab Type \tab Dimension \tab Description \cr
------ \tab ------ \tab ------ \tab ------ \cr
x \tab numeric \tab N \tab covariate used as input for the logistic functions\cr
}
}

}
\keyword{transitions}
