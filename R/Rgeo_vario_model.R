#' Ominidirectional Variogram Model using Rgeostats
#'
#' @description
#' This function calculates an omnidirectional variogram model for a given dataset using
#' 'RGeostats' package. The variogram model can be useful for understanding the spatial dependence
#' structure of the data. The function also allows the user to choose the number of factors to
#' analyze, the property to build the variogram from, and the type of model to use for the
#' variogram.
#'
#' @param db A db-class object. This is the dataset used to calculate the experimental variogram.
#' The variogram is calculated for the set of "z*" variables present in the db.
#' @param num_fac A numeric value indicating how many factors to analyze. This helps to limit the
#' scope of the analysis to a specific number of factors. Default is 3.
#' @param property A string indicating which factor (or property) to build the variogram from.
#' For example, it can be "RC1" or "RC2".
#' @param lag A numeric value or an array containing the distance lags for each calculation
#' direction.
#' If the lag is not defined, set as NA. A default lag is calculated so that the maximum distance is equal
#' to half of the field diagonal.
#' @param nlag A numeric value or an array containing the number of lags for each calculation
#' direction. If nlag is not defined, set it as NA. If the number of lags is not defined, it
#' defaults to 10.
#' @param model A numeric value indicating what type of model to use in the variogram.
#' This parameter corresponds to the model types provided by the RGeostats package.
#' Run the line 'melem.name()' in RGeostats to see the number corresponding to each model.
#'
#' @return An object of class 'model'. This is a plottable Rgeostats omnidirectional variogram
#' model. It can be used for further geostatistical analysis or for visualizing the spatial
#' structure of the data.
#'
#' @export
#'
#' @examples
#' if(requireNamespace("RGeostats")){
#'     library(RGeostats)
#'     data("Klovan_Row80", package = "klovan")
#'     db <- Rgeo_database(Klovan_Row80, 3, "RC3")
#'     model <- Rgeo_vario_model(db, 3, "RC3", lag = 500, model = 13)
#' }
Rgeo_vario_model <- function(db, num_fac, property, lag, nlag = 10, model){
  RC1_vario.omni <- RGeostats::vario.calc(db, lag = lag, nlag = nlag)

  struct <- c(model)

  RC1_data.model.omni <-
    RGeostats::model.auto(
      RC1_vario.omni,
      struct = struct,
      title = paste(property, "Model Omnidirectional"),
      pos.legend = 1,
      xlab = "Log distance",
      ylab = expression(paste("Variance (", gamma, "(h))", sep = ""))
    )


  return(RC1_data.model.omni)
}
