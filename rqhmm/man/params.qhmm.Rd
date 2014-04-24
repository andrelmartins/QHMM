\name{params.qhmm.Rd}
\alias{params.qhmm}
\alias{collect.params.qhmm}
\alias{restore.params.qhmm}
\alias{get.emission.params.qhmm}
\alias{set.emission.params.qhmm}
\alias{get.transition.params.qhmm}
\alias{set.transition.params.qhmm}

\title{QHMM Parameter Functions}
\description{Functions provide a way to manipulate QHMM instance parameters.}

\usage{

collect.params.qhmm(hmm)
restore.params.qhmm(hmm, saved) 

get.emission.params.qhmm(hmm, state, slot = 1)
set.emission.params.qhmm(hmm, state, params, slot = 1, fixed = attr(params, "fixed"))

get.transition.params.qhmm(hmm, state)
set.transition.params.qhmm(hmm, state, params, fixed = attr(params, "fixed"))

}

\arguments{
  \item{hmm}{QHMM instance object}
  \item{saved}{result from \code{collect.params.qhmm}}
  \item{state}{integer representing the HMM state number}
  \item{slot}{integer representing the HMM emission slot (track) number}
  \item{params}{numeric vector of parameter values (length depends on corresponding distribution)}
  \item{fixed}{logical vector indicating which of the parameter values are fixed (unchanged) when running EM}
}

\value{
  \code{set.*.params.qhmm} functions do not return anything, neither does the \code{restore.params.qhmm} function.
  
  \code{get.*.params.qhmm} functions return a numeric vector with parameter values and, as an attribute "fixed", a logical vector of the same length indicating which paramters are currently set to "fixed".
  
  \code{collect.params.qhmm} function returns a nested list of parameter values. First level splits between emissions and transitions, next level splits by states. For transitions, the third level are the parameter vectors (same as those obtainable via \code{get.transition.params.qhmm}). For emissions, the third level splits by slots and the fourth level contains the parameter vectors (same as those obtainable via \code{get.emission.params.qhmm}).
}

\author{André Luís Martins}
\keyword{params}
