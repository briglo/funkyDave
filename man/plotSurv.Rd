% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.r
\name{plotSurv}
\alias{plotSurv}
\title{plotSurv}
\usage{
plotSurv(genes, df = clindat_tcga, expdat = cdat_tcga, trimFirst = F)
}
\arguments{
\item{genes}{vector of gene symbols can be mouse or human, toupper'd to "convert" to human}

\item{df}{an output of cgdsr::getClinicalData with some modified fields}

\item{expdat}{an output of cgdsr::getProfileData with same ID as df}

\item{trimFirst}{logical, perform a hard trim of data first, default=F}
}
\value{
a plot and invisible list of stuff
}
\description{
takes bioportal inputs and tests outcome relative to an elf5 split
}
\examples{
load("r_objects/181206_survivalObjects.rdata")
load("r_objects/180928_HALLMARKandGORSMYgeneListsforGSVA.rdata")
plotSurv(ngs[[sample(1:length(ngs),1)]],clindat,cdat)

}
