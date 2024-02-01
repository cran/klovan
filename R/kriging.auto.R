#' Automatic Kriging Interpolation with Factor Analysis Preprocessing
#'
#' @description
#' This function performs automatic kriging interpolation with factor analysis
#' preprocessing on input data. The optimization may not work as intended use higher
#' num_init_test and num_fin_test values or run the function multiples times to
#' ensure an accurate result.
#'
#' @param data
#' A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis.
#'
#' @param num_fac
#' A numeric value indicating the number of factors to analyze. Default is 3.
#'
#' @param grid_cell_size
#' The desired cell size for the grid. Default is NA, which will calculate the cell size based on the average distance between data points.
#'
#' @param num_init_test
#' The number of random starts for initial model optimization. Default is 8
#'
#' @param num_fin_test
#' The number of random starts for final model optimization. Default is 200
#'
#' @param nugget_bounds
#' A numeric vector specifying the lower and upper bounds for the nugget parameter during optimization. Default is c(0, .2).
#'
#' @param sill_bounds
#' A numeric vector specifying the lower and upper bounds for the sill parameter during optimization. Default is c(0, 20000).
#'
#' @param range_bounds
#' A numeric vector specifying the lower and upper bounds for the range parameter during optimization. Default is c(0, 25000).
#'
#' @return
#' A data frame with interpolated data for the whole grid. Data frame has columns: "C_X",
#' "C_Y", "value", "FA". "C_X" and "C_Y" are the coordinates, "value" is the interpolated value,
#' and "FA" indicates the relevant factor the value corresponds to.
#' @export
#'
#' @examples
#' \donttest{
#' data("Klovan_Row80")
#' kriging.auto(Klovan_Row80)
#' }
#'
#' @importFrom dplyr select_if
#' @importFrom stats prcomp
#' @importFrom stats as.formula
#' @importFrom dplyr select
#' @importFrom tibble as_tibble
#' @importFrom pracma pinv
#' @importFrom sp coordinates
#' @importFrom gstat gstat variogram fit.variogram krige
#' @importFrom sp coordinates<-
#' @importFrom sp gridded<-
#' @importFrom stats ave dist optim rnorm runif
#' @importFrom utils data
#'
kriging.auto <- function(data, num_fac = 3, grid_cell_size = NA, num_init_test = 8, num_fin_test = 200, nugget_bounds=c(0, .2), sill_bounds=c(0, 20000), range_bounds=c(0, 25000)) {

  temp <- data %>% dplyr::select(-tidyselect::starts_with("C_")) %>%
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
  colnames(scores) <- paste0("FAC", 1:num_fac)


  data_scores <- cbind(data, as_tibble(scores[, paste0("FAC", 1:num_fac)]))
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

  # Define a list of models to test
  models_to_test <- list("Sph1" = spherical_model, "Exp1" = exponential_model, "Gau1" = gaussian_model, "Mat1" = matern_model, "Pow1" = power_model, "Quad1" = quadratic_exponential_model, "Card1" = cardinal_sine_model, "Gam1" = gamma_model, "Cau1" = cauchy_model, "Sta1" = stable_model, "Ord1" = order_1_gc_model, "Tri1" = triangle_model, "Cos1" = cosinus_model)

  bin_change <- 0

  # Perform Variogram calculation once for all models
  variogram_calculate <- function(spdf, target_var) {
    dists <- as.matrix(dist(spdf[, c("C_X", "C_Y")]))

    #use Sturges' formula
    n_bins <- (1 + log2(length(dists))) * (.8 - bin_change)

    values <- spdf[[target_var]]
    pairs <- which(upper.tri(dists), arr.ind = TRUE)
    lags <- seq(from = 0, to = max(dists), length.out = n_bins)
    variogram <- data.frame(dist = double(), gamma = double())

    for (i in 1:(length(lags) - 1)) {
      in.bin <- which(dists[pairs] >= lags[i] & dists[pairs] < lags[i + 1])

      if (length(in.bin) > 0) {
        bin.pairs <- pairs[in.bin, , drop = FALSE]
        gamma <- sum((values[bin.pairs[, 1]] - values[bin.pairs[, 2]])^2) / (2 * length(in.bin))
        variogram <- rbind(variogram, data.frame(dist = lags[i], gamma = gamma))
      }
    }

    return(variogram)
  }

  # Function for optimizing variogram parameters
  optimize_variogram_random_starts <- function(target_var, spdf, model, model_name, n_starts, nugget_bounds, sill_bounds, range_bounds) {
    best_result <- Inf
    best_params <- c(NA, NA, NA)
    best_extra_params <- NULL
    variogram <- variogram_calculate(spdf, target_var)
    for (i in 1:n_starts) {
      nugget_init <- runif(1, nugget_bounds[1], nugget_bounds[2])
      sill_init <- runif(1, nugget_bounds[2], sill_bounds[2])
      if (model_name == "Quad1") sill_init <- runif(1, max(variogram$gamma)*.8, max(variogram$gamma)*1.5)
      range_init <- runif(1, sill_init, range_bounds[2])
      if (model_name == "Mat1") range_init <- runif(1, 100, range_bounds[2])

      # Set initial values for a based on the model
      if (model_name == "Quad1") {
        a_values <- c(rnorm(4, mean = 0.1, sd = 0.02), rnorm(4, mean = 1, sd = 0.2), rnorm(4, mean = 10, sd = 2), rnorm(4, mean = 100, sd = 20))
      } else if (model_name == "Mat1") {
        a_values <- c(rnorm(2, mean = 0.5, sd = 0.1), rnorm(2, mean = 1.5, sd = 0.15), rnorm(2, mean = 2.5, sd = 0.25), runif(2, 5, 10), runif(2, 15, 20), runif(2, 20, 50), runif(2, 70, 100))
      } else if (model_name == "Gam1" || model_name == "Pow1") {
        a_values <- runif(16, 1, 2)
      } else if (model_name == "Cau1" || model_name == "Sta1") {
        a_values <- runif(16, 0, 2)
      } else {
        a_values <- 1
      }

      # Define bounds dynamically based on the model_name
      if (model_name == "Quad1") {
        lower_bounds <- c(nugget_bounds[1], max(variogram$gamma)*.8, sill_init, 0.1)
        upper_bounds <- c(nugget_bounds[2], max(variogram$gamma)*1.5, range_bounds[2], 150)
      } else if (model_name == "Mat1") {
        lower_bounds <- c(nugget_bounds[1], nugget_bounds[2], 100, 0.5)
        upper_bounds <- c(nugget_bounds[2], sill_bounds[2], range_bounds[2], 150)
      } else if (model_name == "Gam1") {
        lower_bounds <- c(nugget_bounds[1], sill_bounds[1], range_bounds[1], 1)
        upper_bounds <- c(nugget_bounds[2], sill_bounds[2], range_bounds[2], 2)
      } else if (model_name == "Pow1") {
        lower_bounds <- c(nugget_bounds[1], sill_bounds[1], 0)
        upper_bounds <- c(nugget_bounds[2], sill_bounds[2], 2)
      } else if (model_name == "Card1") {
        lower_bounds <- c(sill_bounds[1])
        upper_bounds <- c(sill_bounds[2])
      } else if (model_name == "Cos1") {
        lower_bounds <- c(sill_bounds[1], 1)
        upper_bounds <- c(sill_bounds[2], 2)
      } else if (model_name == "Cau1" || model_name == "Sta1") {
        lower_bounds <- c(nugget_bounds[1], nugget_bounds[2], range_bounds[1], 0)
        upper_bounds <- c(nugget_bounds[2], sill_bounds[2], range_bounds[2], 2)
      } else {
        lower_bounds <- c(nugget_bounds[1], nugget_bounds[2], range_bounds[1])
        upper_bounds <- c(nugget_bounds[2], sill_bounds[2], range_bounds[2])
      }

      obj_func <- function(params) {
        if (model_name %in% c("Mat1", "Quad1", "Cau1", "Sta1", "Gam1")) {
          model_value <- do.call(model, c(list(h=variogram$dist, nugget=params[1], sill=params[2], range=params[3]), a=params[4]))
        } else if (model_name == "Pow1") {
          model_value <- do.call(model, c(list(h=variogram$dist, nugget=params[1], sill=params[2], a=params[3])))
        } else if (model_name == "Card1") {
          model_value <- do.call(model, c(list(h=variogram$dist, sill=params[1])))
        } else if (model_name == "Cos1") {
          model_value <- do.call(model, c(list(h=variogram$dist, sill=params[1], a=params[2])))
        } else {
          model_value <- do.call(model, list(h=variogram$dist, nugget=params[1], sill=params[2], range=params[3]))
        }
        if (any(is.na(model_value))) {
          mean_value <- mean(model_value, na.rm = TRUE) # Compute mean excluding NA values
          model_value[is.na(model_value)] <- mean_value # Replace NA values with the mean
        }
        return(sum((variogram$gamma - model_value)^2))
      }

      for (A in a_values) {
        if (model_name %in% c("Mat1", "Quad1", "Cau1", "Sta1", "Gam1")) {
          initial_params <- c(nugget_init, sill_init, range_init, A)
        } else if (model_name == "Pow1") {
          initial_params <- c(nugget_init, sill_init, A)
        } else if (model_name == "Card1") {
          initial_params <- c(sill_init)
        } else if (model_name == "Cos1") {
          initial_params <- c(sill_init, A)
        } else {
          initial_params <- c(nugget_init, sill_init, range_init)
        }
        #suppressMessages({
          tryCatch({
            result <- optim(
              initial_params,
              obj_func,
              method = "L-BFGS-B",
              lower = lower_bounds,
              upper = upper_bounds,
              control = list(maxit = 1200, factr = 1e-8)
            )

            if (!is.null(result$value) && result$convergence == 0 && result$value < best_result) {
              best_result <- result$value
              best_params <- result$par
              best_extra_params <- list(a = A)
            }
          }, warning = function(warning_condition) {
            message("Warning in iteration ", i, ": ", warning_condition)
          }, error = function(error_condition) {
            message("Error in iteration ", i, ": ", error_condition)
          })
        #})
      }
    }

    if (is.na(best_params[1])) {
      stop("Optimization did not converge with any random start")
    } else {
      return(list(best_result = best_result, nugget=best_params[1], sill=best_params[2], range=best_params[3], extra_params=best_extra_params))
    }
  }

  get_model_values <- function(h, nugget, sill, range_val, a, model_name) {
    if (model_name %in% c("Mat1", "Quad1", "Cau1", "Sta1", "Gam1")) {
      return(models_to_test[[model_name]](h, nugget, sill, range_val, a))
    } else if (model_name == "Pow1") {
      return(models_to_test[[model_name]](h, nugget, sill, a))
    } else if (model_name == "Card1") {
      return(models_to_test[[model_name]](h, sill))
    } else if (model_name == "Cos1") {
      return(models_to_test[[model_name]](h, sill, a))
    } else{
      return(models_to_test[[model_name]](h, nugget, sill, range_val))
    }
  }

  # Variogram optimization for each target factor
  kriging_results_df <- data.frame(value = numeric(), C_X = numeric(), C_Y = numeric(), FA = character())
  for (target_var in paste0("FAC", 1:num_fac)) {
    while (TRUE) {
      message(target_var)
      results <- list()
      for (model_name in names(models_to_test)) {
        try({
          opt_result <- optimize_variogram_random_starts(target_var = target_var, data_scores, models_to_test[[model_name]], model_name, n_starts = 1, nugget_bounds, sill_bounds, range_bounds)
          results[[model_name]] <- opt_result
        }, silent = TRUE)

      }

      # Extracting the best_result for each model
      best_results <- sapply(results, function(x) x$best_result)

      # Ordering the best_results in ascending order and taking the first three indices
      top_three_indices <- order(best_results)[1:3]

      # Finding the corresponding model names
      top_three_models <- results[top_three_indices]

      message("Optimizing")
      for (model_name in names(top_three_models)) {
        if (top_three_models[[model_name]]$best_result < .2) break
        opt_result <- optimize_variogram_random_starts(target_var, data_scores, models_to_test[[model_name]], model_name, n_starts = num_fin_test, nugget_bounds, sill_bounds, range_bounds)
        if (opt_result$best_result < top_three_models[[model_name]]$best_result)
          top_three_models[[model_name]] <- opt_result
      }
      # Find the model with the lowest best_result
      model_name <- names(which.min(sapply(top_three_models, function(x) x$best_result)))
      message(model_name)
      message(top_three_models[[model_name]])
      model_vals <- top_three_models[[model_name]]

      # Extract the values and assign them to variables
      nugget <- model_vals$nugget
      sill <- model_vals$sill
      range_val <- model_vals$range
      a <- model_vals$extra_params$a

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

      if (any(is.na(C))) {
        message("NA values detected in C. Imputing with mean...\n")

        global_mean <- mean(C, na.rm = TRUE)
        C[is.na(C)] <- global_mean
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
      # If any value in prediction_data$value is greater than 4, increase bin_change and recalculate
      if (any(prediction_data$value > 4)) {
        bin_change <- bin_change + 0.1
        message(("Kriging failed. Retrying"))
      } else {
        message(("Kriging success"))
        # If no value is greater than 4, exit the loop
        break
      }
    }


    # Append prediction_data to the kriging_results_df
    kriging_results_df <- rbind(kriging_results_df, prediction_data)
  }

  return(kriging_results_df)
}

