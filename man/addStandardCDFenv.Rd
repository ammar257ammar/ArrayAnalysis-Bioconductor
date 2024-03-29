% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_processingQC.R
\name{addStandardCDFenv}
\alias{addStandardCDFenv}
\title{Check that a cdf environment is loaded}
\usage{
addStandardCDFenv(Data, overwrite = FALSE)
}
\arguments{
\item{Data}{(Status: required) The raw data object (datatype: AffyBatch)}

\item{overwrite}{(Status: optional, Default: FALSE) Should the cdfName be overwritten if there is
already a value assigned to the object passed to the function (datatype: logical)}
}
\value{
The object with a cdf annotation assigned if found (datatype: AffyBatch)
}
\description{
This function (from functions_processing) makes sure
that a cdf environment is loaded for the current chiptype.
In some cases a cdf environment will already be available
after reading the data with the ReadAffy function,
then the function will detect this and return the object
as is (unless overwrite is set to TRUE).
In case no cdf environment has been assigned,
it will try to search for a suitable one, and add this if found.
If no suitable cdf can be found, a warning will be generated
and the object returned as is.
}
\examples{
#rawData <- addStandardCDFenv(rawData)
}
