\name{estimateIncidence}
\alias{estimateIncidence}
\title{estimates incedence via back calculation}
\usage{
  estimateIncidence(y, pid, gamma = 0, tol = 10^-5,
    verbose = FALSE)
}
\arguments{
  \item{y}{diagnosis counts}

  \item{pid}{a function giving the probability of diagnosis
  given time units from infection}

  \item{gamma}{smoothing penulty}

  \item{tol}{The tolerance for the EM algorithm}

  \item{verbose}{level of verbosity}
}
\description{
  estimates incedence via back calculation
}

