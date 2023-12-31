% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eigen_IRAM.R
\name{Arnoldi_converge}
\alias{Arnoldi_converge}
\title{A implementation of checking convergence of Arnoldi Decomposition}
\usage{
Arnoldi_converge(A, b, V0, H0, r0, n, m, k, First = FALSE, tol = 1e-06)
}
\arguments{
\item{A}{input matrix, to decompose it into AV = VH + r*t(e)}

\item{b}{vector for Krylov-subspace}

\item{V0}{start V}

\item{H0}{start H}

\item{r0}{start r}

\item{n}{dimension of A}

\item{m}{dimension of H}

\item{k}{need k largest eigenvalues}

\item{First}{whether it is the first Arnoldi, before any restart}

\item{tol}{tolerance of convergence}
}
\value{
A list containing convergence, H, V, r or eigenpairs
}
\description{
A implementation of checking convergence of Arnoldi Decomposition
}
