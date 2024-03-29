% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inv_dis_wt.R
\name{inv_dis_wt}
\alias{inv_dis_wt}
\title{Inverse Distance Weighting}
\usage{
inv_dis_wt(data, num_fac = 3)
}
\arguments{
\item{data}{A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
representing the x and y coordinates of the data points and excludes any rank,
ID, or column not for analysis, see README for details}

\item{num_fac}{a numeric value for how many factors to analyze. Recommended to use 3 and default to 3.}
}
\value{
a data frame with interpolated data for the whole graph. Data frame has collumns: "C_X"   "C_Y"   "value" "FA": C_X, C_Y are coordinates and "value" is the value for the "FA" the relevant factor.
}
\description{
This function applies the Inverse Distance Weighting interpolation algorithm
}
\examples{
data("Klovan_Row80")
inv_dis_wt(Klovan_Row80, 4)
inv_dis_wt(Klovan_Row80, 3)

}
