#' Kriging Interpolation for Databases
#'
#' @description
#' Performs kriging interpolation on a provided database using 'RGeostats' methods.
#' The data is interpolated over a grid covering the entire area of interest.
#'
#' @param db
#' A db-class object. Should be created using the `Rgeo_database()` function.
#' @param model
#' An S4 plottable Rgeostats omnidirectional variogram model.
#' Should be created using the `Rgeo_vario_model()` function.
#' @param dx
#' Optional. The grid cell size in the x-direction. If not provided, it is calculated
#' as the average of the ranges in x and y directions divided by 50.
#' @param dy
#' Optional. The grid cell size in the y-direction. If not provided, it is calculated
#' as the average of the ranges in x and y directions divided by 50.
#'
#' @details
#' The `Rgeo_kriging()` function performs kriging interpolation based on the provided
#' database (db) and variogram model (model). The grid cell sizes `dx` and `dy`
#' can be optionally specified or will be automatically determined based on the data.
#' Results can be visualized with the `Rgeo_kriging_plot()` function and summary statistics
#' can be printed by simply calling the returned kriged object.
#'
#' @return
#' A S4 plottable Rgeostats kriged database. Can be plotted using the
#' `Rgeo_kriging_plot` function. Summary statistics for the kriging process
#' can be printed by simply calling the returned dbgrid3 object.
#' @export
#'
#' @examples
#' if(requireNamespace("RGeostats", quietly = TRUE)){
#'     library(RGeostats)
#'     data("Klovan_Row80", package = "klovan")
#'     db <- Rgeo_database(Klovan_Row80, 3, "RC3")
#'     model <- Rgeo_vario_model(db, 3, "RC3", lag = 500, model = 13)
#'     krig <- Rgeo_kriging(db, model)
#'     krig # prints summary statistics for kriging
#' }
#'
Rgeo_kriging <- function(db, model, dx = NA, dy = NA){


  if (is.na(dx) || is.na(dy)) {
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
  }
  # Define the grid size
  if (is.na(dx)){
    dx <- round((ave((max_y - min_y), (max_x - min_x)))/50)
  }else{
    dx <- dx
  }
  if (is.na(dy)){
    dy <- round((ave((max_y - min_y), (max_x - min_x)))/50)
  }else{
    dy <- dy
  }


  neigh.unique <- RGeostats::neigh.create(type = 0, ndim = 2)

  ## Setting the grid cell size
  dbgrid2 <-
    RGeostats::db.grid.init(
      db,
      dcell = c(dx, dy),
    )
  RGeostats::migrate(
    dbin = db,
    dbout = dbgrid2,
    names = c("R*", "D*", "P*"),
    radix = ""
  )

  ## Kriging the variable of interest
  #Code below performs kriging of the selected property on a grid
  dbgrid3 <-
    RGeostats::kriging(
      db,
      dbgrid2,
      model,
      neigh.unique,
      radix = "Omni")

   return(dbgrid3)
}
