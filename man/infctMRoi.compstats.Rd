% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mergeROI.r
\name{infctMRoi.compstats}
\alias{infctMRoi.compstats}
\title{Compute stats for each group of EIC}
\usage{
infctMRoi.compstats(arois, indicatorVec)
}
\arguments{
\item{arois}{Data frame with ROI definition for each sample}

\item{grp}{EIC group column indicator}

\item{sid}{Sample id column indicator}

\item{mz}{average mz column indicator}

\item{mzmin}{mzmin column indicator}

\item{mzmax}{mzmax column indicator}

\item{rtmin}{rtmin column indicator}

\item{rtmax}{rtmax column indicator}
}
\value{
data.frame with summary of each group of EIC
}
\description{
Compute stats for each group of EIC
}
\keyword{internal}
