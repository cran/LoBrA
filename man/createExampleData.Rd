% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initExampleData.R
\name{createExampleData}
\alias{createExampleData}
\title{Create example data set for 'LoBrA'}
\usage{
createExampleData(
  components = c(70, 10, 10, 10),
  samples = 10,
  classes = 2,
  bg = 3,
  timepoints = rep(5, 3),
  myfile = NA
)
}
\arguments{
\item{components}{vector numbers of background and informative components to be created.}

\item{samples}{number of experiments per class}

\item{classes}{number of classes}

\item{bg}{number of background measurements}

\item{timepoints}{number of sample measurements for each spline}

\item{myfile}{filename of the pdf file created. Note: '.pdf' is added automatically.}
}
\value{
final matrix of example data.
}
\description{
Real signals and background noise originating from experimental settings or random events
}
\examples{
\dontrun{
} 
  components = c(21,3,3,3)
  samples = 10
  classes = 2;
  bg = 3; 
  timepoints = rep(5,3)
  p=TRUE;
  longDataExample <- createExampleData(components, samples, classes, bg, 
                                       timepoints)
  dim(longDataExample)
  

}
