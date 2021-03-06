% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/booty_utility.R
\name{surv_calc}
\alias{surv_calc}
\title{Estimate survivals, detections, and tally adult returns}
\usage{
surv_calc(ch, i, nocc, wt, wt_i, phi_p_only, fpc, ...)
}
\arguments{
\item{ch}{Input file made by \code{format_dat()} function.}

\item{i}{Iteration number used by the bootstrap function \code{bootystrapper()}.}

\item{nocc}{Total detection events including the trawl.}

\item{wt}{Indicates whether to weight the sampling probability.}

\item{wt_i}{Indicates whether to calculate the original estimates using weighted probability.}

\item{phi_p_only}{Option to only calculate survivals and detection and not do the adult counts.}

\item{fpc}{Option to adapt finite population correction for survival and detection calculations.}
}
\value{
Survivals, detection and returing adult counts
}
\description{
Estimate survivals, detections, and tally adult returns
}
\examples{
surv_calc(detect_data, i= 1, nocc= 8, wt= 'n', wt_i= 'n', phi_p_only= 'y', fpc= 'y')

}
