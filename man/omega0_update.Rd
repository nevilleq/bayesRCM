% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/omega0_update.R
\name{omega0_update}
\alias{omega0_update}
\title{Function to update group-level sparse covariance matrix Omega_0 / precision matrix Omega_0}
\usage{
omega0_update(omega_0, D, nu, lambda_3)
}
\arguments{
\item{omega_0}{Matrix, current precision matrix Omega}

\item{D}{Matrix, sum_k Omega_k / tau_k}

\item{nu}{Numeric, sum_k tau_k}

\item{lambda_3}{Numeric, penalty on L1 Norm of omega_0}
}
\description{
Function to update group-level sparse covariance matrix Omega_0 / precision matrix Omega_0
}
