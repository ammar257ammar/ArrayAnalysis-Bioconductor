% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_Stat.R
\name{computeAdvancedStatistics}
\alias{computeAdvancedStatistics}
\title{Compute more advanced statistics}
\usage{
computeAdvancedStatistics(
  normDataTable,
  descriptionFile,
  covariates_string,
  interaction_string,
  paired_string,
  plotVolcanoPlot,
  plotVennPlot,
  matfileName = NULL,
  keepAnnotation = FALSE,
  defaultContr = TRUE
)
}
\arguments{
\item{normDataTable}{normalized data object (Status=required)}

\item{descriptionFile}{(Status: required) The data.frame containing the description
file information (column 1: file names; column 2: names
to be used in the plots; column 3: experimental
groups the samples belong to)(datatype: data.frame)}

\item{covariates_string}{(Status=required)}

\item{interaction_string}{(Status=required)}

\item{paired_string}{(Status=required)}

\item{plotVolcanoPlot}{(Status=required)}

\item{plotVennPlot}{(Status=required)}

\item{matfileName}{(datatype=character, Default=NULL)}

\item{keepAnnotation}{(datatype=logical, Default=FALSE)}

\item{defaultContr}{(datatype=logical, Default=TRUE)}
}
\value{
the contrast matrix file
}
\description{
Compute more advanced statistics
}
\examples{
#example here
}
