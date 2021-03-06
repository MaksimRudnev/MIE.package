% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/server.functions.R
\name{MGCFAparameters}
\alias{MGCFAparameters}
\title{Computes MGCFA models for all available groups in the data}
\usage{
MGCFAparameters(
  model = NULL,
  data,
  group = "cntry",
  parameters = "loadings",
  extra.options = NULL,
  shiny = FALSE
)
}
\arguments{
\item{model}{Model in lavaan syntax}

\item{data}{The data}

\item{group}{Character. Grouping variable.}

\item{parameters}{Character. If "loadings" then configural model is fitted and loadings are returned. If "intercepts" then metric invariance model is fitted and intercepts are returned.}

\item{extra.options}{Currently not used}

\item{shiny}{Logical. If it is evaluated in a shiny context. Default is TRUE.}
}
\description{
Computes MGCFA models for all available groups in the data
}
