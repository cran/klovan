#' Inverse Distance Weighting
#'
#' This function applies the Inverse Distance Weighting interpolation algorithm
#'
#' @param data A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis, see README for details
#'
#' @param num_fac a numeric value for how many factors to analyze. Recommended to use 3 and default to 3.
#'
#' @return a data frame with interpolated data for the whole graph. Data frame has collumns: "C_X"   "C_Y"   "value" "FA": C_X, C_Y are coordinates and "value" is the value for the "FA" the relevant factor.
#' @export
#'
#' @examples
#' data("Klovan_Row80")
#' inv_dis_wt(Klovan_Row80, 4)
#' inv_dis_wt(Klovan_Row80, 3)
#'
inv_dis_wt <- function(data, num_fac = 3){
  C_X <- C_Y <- var1.pred <- pc.names <- EigenValues <- NULL

  temp <- data %>% select(!tidyselect::starts_with("C_")) %>%
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
  scores <- scale(temp) %*% invLoadings %>% as_tibble()

  cols_to_select <- paste0("V", 1:num_fac)
  scores_to_select <- scores[, cols_to_select]
  data <- cbind(data, scores_to_select)

  # Ensure the data_scores data.frame is a SpatialPointsDataFrame
  coord_names <- colnames(data  %>% select(starts_with("C_")))
  data_scores <- data
  coord_formula <- as.formula(paste("~", paste(coord_names[1:length(coord_names)], collapse = "+")))
  coordinates(data_scores) <- coord_formula
  message(paste("Using", paste(coord_names[1:length(coord_names)], collapse = " & "), "to make grid"))

  # Calculate min and max for each coordinate, with a 10% buffer
  if (min(data[[coord_names[1]]]) < 0) {
    min_x <- min(data[[coord_names[1]]]) * 1.1
  }else{
    min_x <- min(data[[coord_names[1]]]) * .9
  }
  max_x <- max(data[[coord_names[1]]]) * 1.1

  if (min(data[[coord_names[2]]]) < 0) {
    min_y <- min(data[[coord_names[2]]]) * 1.1
  }else{
    min_y <- min(data[[coord_names[2]]]) * .9
  }
  max_y <- max(data[[coord_names[2]]]) * 1.1

  # Define the grid size
  gridsize <- round((ave((max_y - min_y), (max_x - min_x)))/50)

  # Generate the grid
  grid <- expand.grid(
    seq(from = min_x, to = max_x, by = gridsize),
    seq(from = min_y, to = max_y, by = gridsize)
  )

  # Name the grid columns
  names(grid) <- coord_names[1:2]

  coordinates(grid) <- coord_formula
  gridded(grid) <- TRUE

  for (i in 1:num_fac) {
    points <- gstat::idw(get(paste0("V", i)) ~ 1, locations = data_scores, newdata = grid)
    assign(paste0("FA", i, "_points"), points)
  }

  Interp_Data_IDW <- data.frame(C_X = numeric(), C_Y = numeric(), value = numeric(), FA = character())  # Initialize an empty dataframe
  for (i in 1:num_fac) {
    factor_points <- get(paste0("FA", i, "_points"))  # Retrieve the interpolated points for the current factor

    factor_data <- as_tibble(factor_points) %>%
      dplyr::select(C_X, C_Y, value = var1.pred) %>%
      dplyr::mutate(FA = paste0("FA", i))  # Add a column indicating the factor name (e.g., RC1, RC2, etc.)

    Interp_Data_IDW <- dplyr::union(Interp_Data_IDW, factor_data)  # Union the current factor data with the existing data
  }
  return(Interp_Data_IDW)
}
