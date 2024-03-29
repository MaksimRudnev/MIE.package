% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{stratifiedMI}
\alias{stratifiedMI}
\title{Run clustered measurement invariance tests}
\usage{
stratifiedMI(
  model,
  group,
  data,
  strata,
  parameters = c("loadings", "intercepts"),
  ref = "configural",
  ...
)
}
\arguments{
\item{model}{character, lavaan syntax model}

\item{group}{character, grouping variable}

\item{data}{data frame}

\item{strata}{A list of character vectors of the group names to create strata.}

\item{parameters}{character vector, "all", "loadings", "thresholds", or "intercepts". strata are applied to this subset of parameters.}
}
\description{
Run clustered measurement invariance tests
}
\details{
This function builds a single model with constraints applied to subsets of groups, and compares it to the reference model (less constrained) as well as to global invariance tests.
}
\examples{
stratifiedMI("F =~ v1 + v2 + v3 + v4", 
         group = "country", 
         data = Dat1, 
         strata = list(North = c("Norway", "Denmark", "Finland"), 
                         South = c("Spain", "Portugal", "Italy")
                         )
          )

}
