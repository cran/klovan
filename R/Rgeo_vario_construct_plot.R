#' Experimental Ominidirectional Variogram Plot using Regostats
#'
#' @description
#' This function constructs an Experimental Ominidirectional Variogram using 'Rgeostats'.
#'
#' @param db
#' The db-class containing the data information used to calculate the experimental variogram.
#' The variograms are calculated for the set of "z*" variables present in the db.
#' @param num_fac
#' A numeric value indicating how many factors to analyze. Default is 3.
#' @param property
#' A string indicating which factor to build variogram from e.g. "RC1" or "RC2"
#' @param lag
#' Array containing the distance lags for each calculation direction.
#' If the lag is not defined, set as NA.
#' A default lag is calculated so that the maximum distance is equal to half of the field diagonal
#' @param nlag
#' Array containing the number of lags for each calculation direction
#' If nlag not defined, set as NA.
#' If the number of lags is not defined, it defaults to 10.
#'
#' @return a plottable Rgeostats Experimental Ominidirectional Variogram model
#' @export
#'
#' @examples
#' if(requireNamespace("RGeostats")){
#'     library(RGeostats)
#'     data("Klovan_Row80", package = "klovan")
#'     db <- Rgeo_database(Klovan_Row80, 3, "RC3")
#'     Rgeo_vario_construct_plot(db, 3, "RC3", lag = 500)
#' }
Rgeo_vario_construct_plot <- function(db, num_fac, property, lag, nlag = 10){
  RC1_vario.omni <- RGeostats::vario.calc(db, lag = lag, nlag = nlag)

  RGeostats::plot(
    RC1_vario.omni,
    type = "p",
    pch = 19,
    npairpt = T,
    npairdw = F,
    size = 2,
    title = paste(property, "Experimental Ominidirectional Variogram"),
    xlab = "Lag distance",
    ylab = expression(paste("Variance (", gamma, "(h))",
                            sep = "")))
}
