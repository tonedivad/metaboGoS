% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/molformUtils.R
\name{getMF.sirius}
\alias{getMF.sirius}
\title{Interface to Sirius}
\usage{
getMF.sirius(precmz = NULL, ms1spec = NULL, ms2spec = NULL,
  ionProd = "[M+H]1+", ppmPrecMF = 11, maxMF = NA, limint = 10^-5,
  exec = "/media/D01/Metabo/MetSoft/sirius3.4.1/bin/sirius3",
  tempdir = "./", listmf = c(), retAnn = FALSE, clean = TRUE)
}
\arguments{
\item{precmz}{MS1 spectra}

\item{ms1spec}{MS1 spectra}

\item{ms2spec}{MS/MS spectra}

\item{ionProd}{Ionisation product}

\item{ppmPrecMF}{ppmPrecMF}

\item{maxMF}{Maximum number of MF to be returned}

\item{limint}{Minimum intensity perc.}

\item{exec}{Sirius executable}

\item{tempdir}{Directory containing the output}

\item{listmf}{Annotation of a single list of mol formula}

\item{retAnn}{Add annotation in the original MS/MS}

\item{clean}{clean-up the mess}
}
\value{
List of stuff
}
\description{
Interface to Sirus
}