% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LoBrA.R
\name{selectComponents}
\alias{selectComponents}
\title{Create a new 'LDO' Object that only contains the selected components.}
\usage{
selectComponents(ldo, components, name = paste(ldo@name, " selected"))
}
\arguments{
\item{ldo}{Longitudinal Data Object}

\item{components}{Component names to select for the new ldo object. Only elements from this list that overlap with the peak names in the given ldo, are utilized.}

\item{name}{Name of newly created 'LDO' object.}
}
\value{
new ldo object only containing the selected components.
}
\description{
Create a new 'LDO' Object that only contains the selected components.
}
