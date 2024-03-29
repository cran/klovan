% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rgeo_vario_model.R
\name{Rgeo_vario_model}
\alias{Rgeo_vario_model}
\title{Ominidirectional Variogram Model using Rgeostats}
\usage{
Rgeo_vario_model(db, num_fac, property, lag, nlag = 10, model)
}
\arguments{
\item{db}{A db-class object. This is the dataset used to calculate the experimental variogram.
The variogram is calculated for the set of "z*" variables present in the db.}

\item{num_fac}{A numeric value indicating how many factors to analyze. This helps to limit the
scope of the analysis to a specific number of factors. Default is 3.}

\item{property}{A string indicating which factor (or property) to build the variogram from.
For example, it can be "RC1" or "RC2".}

\item{lag}{A numeric value or an array containing the distance lags for each calculation
direction.
If the lag is not defined, set as NA. A default lag is calculated so that the maximum distance is equal
to half of the field diagonal.}

\item{nlag}{A numeric value or an array containing the number of lags for each calculation
direction. If nlag is not defined, set it as NA. If the number of lags is not defined, it
defaults to 10.}

\item{model}{A numeric value indicating what type of model to use in the variogram.
This parameter corresponds to the model types provided by the RGeostats package.
Run the line 'melem.name()' in RGeostats to see the number corresponding to each model.}
}
\value{
An object of class 'model'. This is a plottable Rgeostats omnidirectional variogram
model. It can be used for further geostatistical analysis or for visualizing the spatial
structure of the data.
}
\description{
This function calculates an omnidirectional variogram model for a given dataset using
'RGeostats' package. The variogram model can be useful for understanding the spatial dependence
structure of the data. The function also allows the user to choose the number of factors to
analyze, the property to build the variogram from, and the type of model to use for the
variogram.
}
\examples{
if(requireNamespace("RGeostats")){
    library(RGeostats)
    data("Klovan_Row80", package = "klovan")
    db <- Rgeo_database(Klovan_Row80, 3, "RC3")
    model <- Rgeo_vario_model(db, 3, "RC3", lag = 500, model = 13)
}
}
