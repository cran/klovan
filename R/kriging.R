#' Perform Kriging Interpolation
#'
#' This function performs kriging interpolation on spatial data using ridge regression to calculate
#' the kriging weights. It uses either regular inverse or generalized inverse with ridge regression
#' based on the availability of regular inverse for the given covariance matrix.
#'
#' @param data
#' A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis, see README for details
#'
#' @param factor The target factor (FAC) to be interpolated using kriging.
#' @param grid_cell_size
#' The desired cell size for the grid. Default is NA, which will calculate the cell size based on the average distance between data points.
#' @param nugget The nugget effect parameter for the variogram model.
#' @param sill The sill parameter for the variogram model.
#' @param range_val The range parameter for the variogram model.
#' @param a Additional parameter (depends on the variogram model) use NA if not needed.
#' @param model_name The name of the model to use for variogram fitting and kriging. Options include
#' "Sph1", "Exp1", "Gau1", "Mat1", "Pow1", "Quad1", "Card1", "Gam1", "Cau1", "Sta1", "Ord1", "Tri1", and "Cos1".
#' use function`print_model_names()` for more information
#'
#' @return A data frame containing the interpolated values for the target factor (FAC).
#'
#' @examples
#' data(Klovan_Row80)
#' # Perform kriging interpolation for FAC1
#' kriging_results <- kriging(Klovan_Row80, factor = 1, grid_cell_size = NA,
#' nugget=.0001, sill=2.5, range_val=1000, a=NA, model_name="Sph1")
#'
#' @importFrom dplyr select
#' @importFrom dplyr %>% select_if
#' @importFrom stats prcomp
#' @importFrom pracma pinv
#' @importFrom MASS ginv
#'
#' @export
kriging <- function(data, factor, grid_cell_size = NA, nugget, sill, range_val, a, model_name) {

  temp <- data %>% dplyr::select(-tidyselect::starts_with("C_")) %>%
    select_if(~{
      is.numeric(.) &&
        length(unique(.)) > 1 &&
        !all(diff(sort(.)) == 1)
    })

  #preform PCA
  pca <- stats::prcomp(temp, center = TRUE, scale. = TRUE)

  #preform FA
  rawLoadings <- pca$rotation[, 1:length(temp)] %*% diag(pca$sdev, length(temp), length(temp))
  rotatedLoadings <- varimax(rawLoadings)$loadings
  invLoadings <- t(pracma::pinv(rotatedLoadings))
  scores <- scale(temp) %*% invLoadings %>% as_tibble()
  colnames(scores) <- paste0("FAC", 1:length(temp))


  data_scores <- cbind(data, as_tibble(scores[, paste0("FAC", 1:length(temp))]))
  coord_names <- colnames(data)[grepl("^C_", colnames(data))]

  # Define the grid size
  if (is.na(grid_cell_size)) {
    dist_x <- diff(range(data[[coord_names[1]]]))
    dist_y <- diff(range(data[[coord_names[2]]]))
    gridsize <- max(dist_x, dist_y) / 50
  } else {
    gridsize <- grid_cell_size
  }

  # Generate the grid
  grid <- expand.grid(
    C_X = seq(from = min(data[[coord_names[1]]]), to = max(data[[coord_names[1]]]), by = gridsize),
    C_Y = seq(from = min(data[[coord_names[2]]]), to = max(data[[coord_names[2]]]), by = gridsize)
  )

  # Corrected lines for n_rows and n_cols
  n_rows <- length(unique(grid$C_Y))
  n_cols <- length(unique(grid$C_X))

  target_var <- paste0("FAC", factor)

  model1 <- get_model(model_name)

  get_model_values <- function(h, nugget, sill, range_val, a, model_name) {
    if (model_name %in% c("Mat1", "Quad1", "Cau1", "Sta1", "Gam1")) {
      return(model1(h, nugget, sill, range_val, a))
    } else if (model_name == "Pow1") {
      return(model1(h, nugget, sill, a))
    } else if (model_name == "Card1") {
      return(model1(h, sill))
    } else if (model_name == "Cos1") {
      return(model1(h, sill, a))
    } else{
      return(model1(h, nugget, sill, range_val))
    }
  }

  # Variogram optimization for each target factor
  kriging_results_df <- data.frame(value = numeric(), C_X = numeric(), C_Y = numeric(), FA = character())

  # Calculate kriging weights using ridge regression
  coordinates <- as.matrix(data_scores[, c("C_X", "C_Y")])
  n <- nrow(coordinates)
  C <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      h <- sqrt(sum((coordinates[i,] - coordinates[j,])^2))
      C[i, j] <- get_model_values(h, nugget, sill, range_val, a, model_name)
    }
  }

  y <- as.matrix(data_scores[target_var])
  y <- rbind(y, 0) # Append 0 to y
  C <- rbind(cbind(C, rep(1, n)), c(rep(1, n), 0))

  lambda_values <- 10^seq(-10, 0, length.out=10)

  train_idx <- sample(1:n, 0.8 * n)
  valid_idx <- setdiff(1:n, train_idx)

  compute_validation_error <- function(lambda_ridge) {
    C_train <- C[train_idx, train_idx] + lambda_ridge * diag(length(train_idx))
    y_train <- y[train_idx]
    weights <- solve(C_train, y_train)

    C_valid <- C[valid_idx, train_idx]
    y_valid <- y[valid_idx]
    predictions <- C_valid %*% weights
    return(mean((y_valid - predictions)^2)) # MSE as validation error
  }

  validation_errors <- sapply(lambda_values, compute_validation_error)
  optimal_lambda_ridge <- lambda_values[which.min(validation_errors)]

  C_ridge <- C + optimal_lambda_ridge * diag(n + 1)

  tryCatch({
    # Try to solve using regular inverse
    lambda <- solve(C_ridge, y)
  }, error = function(error_condition) {
    # If error, use generalized inverse with ridge regression
    lambda <- MASS::ginv(C_ridge) %*% y
  })

  weights <- lambda[1:n]

  # Perform interpolation on the grid
  prediction_data <- data.frame(value = numeric(n_rows * n_cols),
                                C_X = numeric(n_rows * n_cols),
                                C_Y = numeric(n_rows * n_cols),
                                FA = target_var)
  for (i in 1:nrow(grid)) {
    new_location <- as.numeric(grid[i,])
    covariances <- sapply(1:n, function(j) {
      h <- sqrt(sum((new_location - coordinates[j,])^2))
      get_model_values(h, nugget, sill, range_val, a, model_name)
    })
    prediction <- sum(weights * covariances) + lambda[n + 1]

    prediction_data$value[i] <- prediction
    prediction_data$C_X[i] <- new_location[1]
    prediction_data$C_Y[i] <- new_location[2]
  }
  return(prediction_data)
}
