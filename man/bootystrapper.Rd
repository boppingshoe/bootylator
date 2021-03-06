% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/booty_utility.R
\name{bootystrapper}
\alias{bootystrapper}
\title{Bootstrap using surv_calc and organize output}
\usage{
bootystrapper(d, fn, iter, wgt = "n", wgt_init = "n",
  phi_p_only = "n", fpc = "y", logit_link = "n",
  save_name = "none", ...)
}
\arguments{
\item{d}{Input file made by \code{format_dat()}.}

\item{fn}{Function to run the bootstrap on.}

\item{iter}{Amount of bootstrap iterations.}

\item{wgt}{Indicates whether to weight the sampling probability. Default is no ("n").}

\item{wgt_init}{Indicates whether to calculate the original estimates using weighted probability. Default is no ("n").}

\item{phi_p_only}{Indicate to turn off the phi_p_only option in \code{curv_calc()}. Default is no ("n").}

\item{fpc}{Indicate to turn off the finite population correction option in \code{curv_calc()}. Default is yes ("y").}

\item{logit_link}{Indicate to use "RMark" or "marked" and estimate using logit link. The default here is none ("n").}

\item{save_name}{Name to save bootstrap output in CSSOUTPUT in SQL server. No results will be saved if nothing is specified.}
}
\value{
Estimates in a data frame with original estimate as the first row and bootstrap results in the remaining rows.
}
\description{
Bootstrap using surv_calc and organize output
}
\examples{
# To conduct standard CSS bootstrap procedures
out1<- bootystrapper(detect_data, surv_calc, iter= 1000, fpc= 'n', save_name='SR HCH 2008 CATH')
# To conduct weighted bootstrap and produce only survival and detection estimates
out2<- bootystrapper(detect_data, surv_calc, iter= 100, wgt= 'y', wgt_init= 'y', phi_p_only= 'y')
head(out2)

}
