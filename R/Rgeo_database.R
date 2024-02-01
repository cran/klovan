#' Database Construction
#'
#' @description
#' Constructs a database from a provided dataset using specified factors.
#' For more details on the dataset format, see the package README.
#'
#' @param data
#' A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis, see README for details
#' @param num_fac
#' A numeric value indicating the number of factors to analyze. Default is 3.
#' @param property
#' A string indicating which factor to build variogram from e.g. "RC1" or "RC2"
#'
#' @details
#' The `Rgeo_database` function constructs a db-class object from the provided
#' dataset using the number of factors specified by `num_fac` and made for use with `property`.
#'
#' @return
#' A db-class object containing the factors selected with `num_fac` and made for use with `property`.
#' @export
#'
#' @examples
#' if(requireNamespace("RGeostats", quietly = TRUE)){
#'     library(RGeostats)
#'     data("Klovan_Row80", package = "klovan")
#'     Rgeo_database(Klovan_Row80, 3, "RC3")
#' }
#'
Rgeo_database <- function(data, num_fac = 3, property) {

  temp <- data %>% dplyr::select(!tidyselect::starts_with("C_")) %>%
    select_if(~{
      is.numeric(.) &&
        length(unique(.)) > 1 &&
        !all(diff(sort(.)) == 1)
    })

  #preform PCA
  pca <- stats::prcomp(temp, center = TRUE, scale. = TRUE)

  #preform FA
  rawLoadings <- pca$rotation[, 1:num_fac] %*% diag(pca$sdev, num_fac, num_fac)
  rotatedLoadings <- varimax(rawLoadings)$loadings
  invLoadings <- t(pracma::pinv(rotatedLoadings))
  scores <- scale(temp) %*% invLoadings %>% dplyr::as_tibble()

  cols_to_select <- paste0("V", 1:num_fac)
  scores_to_select <- scores[, cols_to_select]
  data <- cbind(data, scores_to_select)

  Final_Scores <- dplyr::select(data, "C_X", "C_Y", rep(paste0("V", 1:num_fac), num_fac))
  new_col_names <- rep(paste0("RC", 1:num_fac), 1)
  colnames(Final_Scores)[3:length(colnames(Final_Scores))] <- new_col_names

  xlon <- "C_X"
  ylat <- "C_Y"

  db <-
    Final_Scores %>%

    RGeostats::db.create() %>%
    RGeostats::db.locate(c(xlon, ylat), "x") %>%

    RGeostats::db.locate(names = property, loctype = "z")
  return(db)
}

