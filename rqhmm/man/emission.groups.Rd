\name{emission.groups}
\alias{new.emission.groups}
\alias{add.emission.groups}

\title{Create QHMM Emission Groups}
\description{Functions provide a way to define which emission state/slot pairs should share parameters.}

\usage{

new.emission.groups(nstates, nslots)
add.emission.groups(emission_sharing_matrix, group=NULL, states=NULL, slots=NULL)

}

\arguments{
  \item{nstates}{integer, number of states in the HMM}
  \item{nslots}{integer, number of slots (emission tracks) in the HMM}
  \item{emission_sharing_matrix}{output of \code{new.emission.groups} and \code{add.emission.groups}, i.e., matrix that defines the state/slot pairs should share emission distribution parameters.}
  \item{group}{Group is defined in one of two alternative ways: (1) an integer vector (slot number followed by two or more state numbers); (2) a list with two integer vectors (slots) and (states)}
  \item{states}{integer vector of state numbers (ranging from 1 to 'nstates')}
  \item{slots}{integer vector of slots numbers (ranging from 1 to 'nslots')}
}

\details{
  Either supply the group definition using the \code{group} argument, or both the \code{states} and the \code{slots} arguments.
  If following the second approach, \code{states} and the \code{slots} arguments define an ordered set of pairs (states[i], slots[i]) that define the elements of the group.

  \subsection{Constraints}
  
  For a group to be valid, it must also have the following (non-enforced) properties:
  \item{Emission functions MUST be set the same within groups.}
  \item{Slot dimensions MUST be the same.}
}

\value{
  Matrix of States (rows) X Slices (cols) with the following properties:
  \item{0s specify no grouping}
  \item{Groups are specified as sequential integers (1:NGroups).}
  
  See "Details" for other properties of the HMM definition that must be true for valid groups.
}

\author{André Luís Martins}

\seealso{new.qhmm}

\keyword{qhmm}
