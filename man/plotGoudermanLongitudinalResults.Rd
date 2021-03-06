% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gaudermanModelling.R
\name{plotGoudermanLongitudinalResults}
\alias{plotGoudermanLongitudinalResults}
\title{Plotting the 'Gouderman' LME Model and Results.}
\usage{
plotGoudermanLongitudinalResults(
  evaluationresult,
  main = "Mixed Effect Spline Model Evaluation",
  ylab = "Value",
  xlab = "Time",
  peaknames = NULL
)
}
\arguments{
\item{evaluationresult}{'GaudermanModelEvaluation' data object, created by the modelGoudermanLongitudinal function.}

\item{main}{title of the plot}

\item{ylab}{y axis label}

\item{xlab}{x axis label}

\item{peaknames}{selection of peaks to be plotted}
}
\value{
No return value
}
\description{
Plotting the 'Gouderman' LME Model and Results.
}
\examples{
  
  wd <- tempdir()
  data(LoBraExample)
  selectedLDO <- selectComponents(ldo, components)
  gaudermanLDOexample <- createGoudermanData(selectedLDO, breaks=c(8, 12), center=12, timeperiod=2)
  evalResult<- modelGoudermanLongitudinal(gaudermanLDOexample)
  # Plot all peaks
  filename<- file.path(wd, "finalModelEvaluation.pdf") ;
  oldpar <- par("mfrow")
  grDevices::pdf(filename, width=16, height=8);
    graphics::par(mfrow=c(1,1));
    plotGoudermanLongitudinalResults(evalResult);
  par(mfrow = oldpar)
  grDevices::dev.off();
  
  #Plot a selection of Peaks
  peaknames<- evalResult@gaudermanLDO@peaknames;
  filename<- file.path(wd, "finalModelEvaluation-components.pdf") ;
  oldpar <- par("mfrow")
  grDevices::pdf(filename, width=20, height=8);
    graphics::par(mfrow=c(2,5));
    plotGoudermanLongitudinalResults(evalResult, main="", peaknames=peaknames);
  par(mfrow = oldpar)
  grDevices::dev.off();
 
}
