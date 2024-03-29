% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_processingQC.R
\name{deduceSpecies}
\alias{deduceSpecies}
\title{Find the species of the chiptype}
\usage{
deduceSpecies(descr = NULL)
}
\arguments{
\item{descr}{(Status: required, Default: NULL) A string indicating the chiptype,
which can be obtained by getting the @annotation slot from an AffyBatch object (datatype: character)}
}
\value{
The species associated with the current chiptype, or ""
if detection was unsuccessful (datatype: character)
}
\description{
This function (from functions_processing) tries to determine the species
related to the current chiptype, if the species has not been provided
by the user (and as such is set to ""). If the descr parameter
is not provided or is empty, an empty string is returned as species.
In other cases the function tries to load an annotation library
depending on the chiptype to find the species. If not successful,
it will be set by hand for some predefined chiptypes.
If still not successful, the empty string is returned.
}
\examples{
#species <- deduceSpecies(rawData@annotation)
}
