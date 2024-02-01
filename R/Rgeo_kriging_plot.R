#' Plotting Kriged Database
#'
#' @description
#' This function creates a plot of the kriged database.
#' The function is designed specifically for visualizing geostatistical data.
#'
#' @param krig_db
#' A kriged database object, resulting from the `Rgeo_kriging` function.
#' @param db
#' The original database object that was used to generate the kriged database.
#' @param property
#' A character string representing the property (or column name) from the database that you want to visualize e.g. "RC1" or "RC2".
#'
#' @return
#' A plot comparing the specified `property` in the original and kriged databases.
#'
#' @details
#' The function takes a kriged database and the original database, then generates a comparative plot for a specific property.
#' This helps in understanding the effect of kriging on the selected property.
#'
#' @export
#'
#' @examples
#' if(requireNamespace("RGeostats", quietly = TRUE)){
#'     library(RGeostats)
#'     data("Klovan_Row80", package = "klovan")
#'     db <- Rgeo_database(Klovan_Row80, 3, "RC3")
#'     model <- Rgeo_vario_model(db, 3, "RC3", lag = 500, model = 13)
#'     krig_db <- Rgeo_kriging(db, model)
#'     Rgeo_kriging_plot(krig_db, db, "RC3")
#' }
Rgeo_kriging_plot <- function(krig_db, db, property) {

  #find min/max
  coord_set <- db@items %>% dplyr::select(tidyselect::starts_with("C_"))
  if (min(coord_set[1]) < 0) {
    min_x <- min(coord_set[1]) * 1.15
  }else{
    min_x <- min(coord_set[1]) * .85
  }
  max_x <- max(coord_set[1]) * 1.15

  if (min(coord_set[2]) < 0) {
    min_y <- min(coord_set[2]) * 1.15
  }else{
    min_y <- min(coord_set[2]) * .85
  }
  max_y <- max(coord_set[2]) * 1.15


   krig_db <- RGeostats::db.locate(krig_db, 4, "z")

   RGeostats::plot(
    krig_db,
    pos.legend = 3,
    cex = .5,
    xlim = c(min_x,max_x),
    ylim = c(min_y,max_y),
    title =
      paste(property, "Kriging with omni-directional Model")
  )

  # Adding the contours
   RGeostats::plot(
    krig_db,
    name.contour = colnames(krig_db@items[4])[1],
    nlevels = 7,
    col = ("black") ,
    add = TRUE
  )

  # Adding the data locations
   RGeostats::plot(
    db,
    add = TRUE,
    pch = 21,
    name.color = property,
    col = "black")
}
