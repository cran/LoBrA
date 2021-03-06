% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initExampleData.R
\name{createInformativeComponents}
\alias{createInformativeComponents}
\title{Simulate informative peaks}
\usage{
createInformativeComponents(
  components,
  samples,
  labels,
  timepoints = c(5, 5, 5),
  bg = 3,
  mean = 5,
  sd = 3,
  segment = 1,
  slopeSD = 2,
  randomnoise = 0.5,
  plotting = FALSE
)
}
\arguments{
\item{components}{number of background components to be created}

\item{samples}{number of experiments}

\item{labels}{label of each experiment}

\item{timepoints}{number of sample measurements}

\item{bg}{number of background measurements}

\item{mean}{mean value of noise for the intercept of this components}

\item{sd}{standard deviation value of noise for the intercept of this component}

\item{segment}{indicating the segment, that will have an informative event (changing slope for one class)}

\item{slopeSD}{standard deviation value for the generated slope of for this component (mean is zero, therefore, the slope can be either negative or positive)}

\item{randomnoise}{random variation changing at each time point}

\item{plotting}{logical value, (default is FALSE), if TRUE the function will plot the created time series.}
}
\value{
matrix of informative components
}
\description{
This function simulates signals correlated to different informative events.
}
