% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covar_mtrx.R
\name{covar_mtrx}
\alias{covar_mtrx}
\title{Create Co-variance Matrix}
\usage{
covar_mtrx(data)
}
\arguments{
\item{data}{A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
representing the x and y coordinates of the data points and excludes any rank,
ID, or column not for analysis.}
}
\value{
A non-normalized co-variance matrix of the klovan data.
}
\description{
This function creates a non-normalized co-variance matrix
from the given klovan dataset. For further details on klovan datasets,
refer to the README.
}
\examples{
data("Klovan_Row80")
covar_mtrx(Klovan_Row80) # view co-variance matrix
}
