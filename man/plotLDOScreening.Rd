% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/screening.R
\name{plotLDOScreening}
\alias{plotLDOScreening}
\title{Plotting the screening results.}
\usage{
plotLDOScreening(
  ldoscreen,
  plotAll = FALSE,
  correctionmethod = "levene",
  decs = 3,
  ask = FALSE,
  peaknames = rownames(ldoscreen@selectedPeaks)
)
}
\arguments{
\item{ldoscreen}{LDO screening result}

\item{plotAll}{Select all components to be plotted. Default plots only the selected peaks using the correction method.}

\item{correctionmethod}{Version of correction method to be used to select the peaks. Valid values are 'bf', 'levene', and 'bartlett'.}

\item{decs}{decimal numbers of p-values to be plotted.}

\item{ask}{logical. Modifies the graphical parameter \code{ask} in \code{par} (If TRUE (and the R session is interactive) the user is asked for input, before a new figure is drawn. As this applies to the device, it also affects output by packages grid and lattice. It can be set even on non-screen devices but may have no effect there.)}

\item{peaknames}{Defining a list of peaks to be plotted. By default all peaks will be plotted.}
}
\value{
No return value
}
\description{
For each peak two box plots are created. The first plot shows a boxplot of the Sample Intercept Comparison of the sample and the background, and the corresponding p-values. The second plot shows a boxplot of the Residual Comparison of the sample and the background, and the corresponding p-values.
}
\examples{
\dontrun{
} 

  wd <- tempdir()
  data(LoBraExample)
  ldos<-screening(ldo, method= c('levene'), alpha =0.05, criteria=c(1,1))
  filename<- file.path(wd, "screeningresults.pdf") 
  grDevices::pdf(filename, width=16, height=8)
  plotLDOScreening(ldos)
  grDevices::dev.off();

}
