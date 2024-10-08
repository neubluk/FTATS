% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reconcile.R
\name{reconcile_forecasts.cv}
\alias{reconcile_forecasts.cv}
\title{reconcile_forecasts.cv}
\usage{
reconcile_forecasts.cv(
  y,
  yhat,
  k,
  var_type = c("var", "struc", "series"),
  cov_type = c("cov_shrink", "inv_glasso", "corr_shrink", "spectral"),
  tune_grid,
  yhat_test = NULL,
  ytest = NULL,
  full_cov = NULL,
  S = NULL,
  test_recon_cum = FALSE,
  err_measure = function(x, xtilde) mean((x - xtilde)^2),
  ...
)
}
\arguments{
\item{y}{data matrix representing the hierarchy of time series}

\item{yhat}{data matrix of base forecasts with same dimension as y}

\item{k}{vector of aggregation levels}

\item{var_type}{type of variance to be used on each aggregation level}

\item{cov_type}{type of covariance to be used between various aggregation levels}

\item{tune_grid}{matrix with hyperparameters as variables, and combinations in rows}

\item{yhat_test}{optional base test forecasts with same number of columns as y and yhat}

\item{ytest}{data matrix of test data in order to be able to use test_recon_cum=TRUE}

\item{full_cov}{if not given, it is estimated using data y and yhat}

\item{S}{summing matrix. if not given, it is computed according to k}

\item{test_recon_cum}{should the reonciliation method be re-estimated cumulatively on the test data (T/F)?}

\item{err_measure}{error function f(x,xtilde) to be minimized and choose optimal hyperparameters. default MSE.}

\item{...}{additional parameters for reconcile_forecasts()}
}
\value{
list of optimal recon result, and cv results
}
\description{
Cross-Validation for selected hyperparameters for recon. methods where necessary
}
\details{
cov_shrink and corr_shrink types need additional lambda parameter between 0 and 1.
inv_glasso requires rho>0 parameter.
spectral scaling requires lambda and neig parameters.
}
\seealso{
\code{\link[=reconcile_forecasts]{reconcile_forecasts()}}
}
