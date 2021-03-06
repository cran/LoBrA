% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train.peaksplines.R
\name{lobraModelSelection}
\alias{lobraModelSelection}
\title{Evaluation of different spline variants.}
\usage{
lobraModelSelection(
  selectedLDO,
  potentialBreaks = c(),
  nknots = c(0, 1, 2),
  splinetype = "linear",
  qualityMeasure = c("AIC", "BIC", "logLik")
)
}
\arguments{
\item{selectedLDO}{\code{\link{LDO}} containing all selected metabolites to be used for the model selection.}

\item{potentialBreaks}{Vector of all possible knots to be used for the spline modeling.}

\item{nknots}{Vector of number of spline knots to be used.  Therefore, 0 ~ no spline, 1 ~ spline with one knot, 2 ~ spline with two knots, etc.}

\item{splinetype}{spline type default is 'linear'. (Currently only linear is supported.)}

\item{qualityMeasure}{Vector of quality measures to be used. Possible options are 'AIC', 'BIC', and 'logLik'.}
}
\value{
\code{\link{LDOmodelselection}} Object.
For each quality measure the model list contains a list of models for each spline tested. Additionally, the output contains a matrix of qualities for each Spline Component pair. And finally there is a list of breaks for each spline tested.
}
\description{
The model selection method evaluates which spline models achieve the best quality among all tested metabolites.
}
\examples{
\dontrun{
} 
  data(LoBraExample)
  potentialBreaks <- c(8,12)
  selectedLDO <- selectComponents(ldo, components)
  ldoSelect<- lobraModelSelection(selectedLDO, potentialBreaks, nknots=c( 1, 2))
  length(ldoSelect@ldo@peaknames)
  
  
}
