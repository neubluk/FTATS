% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/agg_arima.R
\name{fitted.agg_arima}
\alias{fitted.agg_arima}
\title{Obtain fitted values for agg. ARIMA models}
\usage{
\method{fitted}{agg_arima}(object, h = 1)
}
\arguments{
\item{object}{result of fit_agg_arima}

\item{h}{horizon used to internally compute fitted values}
}
\value{
list of data.frames with variables index and x
}
\description{
Obtain fitted values for agg. ARIMA models
}
