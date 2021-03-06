% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regression.R
\name{solveWtCrit}
\alias{solveWtCrit}
\title{Correlation between weighted criterion composite and predictors.}
\usage{
solveWtCrit(ryy, rxy, wt)
}
\arguments{
\item{ryy}{A matrix of criterion intercorrelations.}

\item{rxy}{A vector of predictor criterion intercorrelations, or a
matrix with one row per criterion.}

\item{wt}{A vector of criterion weights, or a matrix with one set of 
criterion weights per row.}
}
\value{
A matrix of correlation coefficent with one row per weight vector
        and one column per predictor.
}
\description{
Correlation between weighted criterion composite and predictors.
}
\examples{
library(iopsych)
data(dls2007)
dat <- dls2007[1:6, 2:7]
ryy <- dat[5:6, 5:6]
rxy <- dat[5:6, 1:4]

wt1 <- c(.25, .75)
wt2 <- c(.75, .25)
wt_mat <- rbind(wt1, wt2)

solveWtCrit(ryy=ryy, rxy=rxy, wt=wt_mat)
}
\author{
Allen Goebl Jeff Jones
}
