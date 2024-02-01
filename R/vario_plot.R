#' Create Variogram Plot
#'
#' This function calculates the empirical variogram for a given target factor (FAC)
#' and plots it along with the fitted variogram based on the specified variogram model.
#'
#' @param data A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis, see README for details
#' @param factor The target factor (FAC) to be used for variogram calculation and fitting.
#' @param nlags The number of lag bins for variogram calculation. Default is NA, which will
#' use Sturges' formula to determine the number of lags.
#' @param lags The lag width for variogram calculation. Default is NA, which will calculate
#' the lag width based on the range of distances.
#' @param nugget The nugget effect parameter for the variogram model.
#' @param sill The sill parameter for the variogram model.
#' @param range_val The range parameter for the variogram model.
#' @param a Additional parameter (depends on the variogram model) use NA if not needed.
#' @param model_name The name of the model to use for variogram fitting.
#' Available options include "Sph1", "Exp1", "Gau1", "Mat1", "Pow1", "Quad1", "Card1", "Gam1", "Cau1",
#' "Sta1", "Ord1", "Tri1", and "Cos1". Use function `print_model_names()` for more information.
#'
#' @return A plot displaying the empirical variogram and the fitted variogram based on the specified model.
#'
#' @examples
#' data(Klovan_Row80)
#' # Plot variogram for FAC1
#' vario_plot(Klovan_Row80, factor = 1, nlags = 10, nugget = 0.01, sill = 2.5,
#' range_val = 1000, a = NA, model_name = "Sph1")
#'
#' @importFrom dplyr select
#' @importFrom dplyr %>% select_if
#' @importFrom stats prcomp
#' @importFrom pracma pinv
#' @importFrom MASS ginv
#' @import ggplot2
#'
#' @export
vario_plot <- function(data, factor, nlags = NA, lags = NA, nugget, sill, range_val, a, model_name) {

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

  target_var <- paste0("FAC", factor)

  variogram_calculate <- function(spdf, target_var, nlags = NA, lags = NA) {
    dists <- as.matrix(dist(spdf[, c("C_X", "C_Y")]))

    # Use Sturges' formula
    if (is.na(nlags)) {
      nlags <- 1 + log2(length(dists))
    }

    if (is.na(lags)) {
      lags <- (max(dists) - min(dists)) / nlags
    }

    values <- spdf[[target_var]]
    pairs <- which(upper.tri(dists), arr.ind = TRUE)
    lag_seq <- seq(from = 0, to = max(dists), by = lags)
    variogram <- data.frame(dist = double(), gamma = double())

    for (i in 1:(length(lag_seq) - 1)) {
      in.bin <- which(dists[pairs] >= lag_seq[i] & dists[pairs] < lag_seq[i + 1])

      if (length(in.bin) > 0) {
        bin.pairs <- pairs[in.bin, , drop = FALSE]
        gamma <- sum((values[bin.pairs[, 1]] - values[bin.pairs[, 2]])^2) / (2 * length(in.bin))
        variogram <- rbind(variogram, data.frame(dist = lag_seq[i], gamma = gamma))
      }
    }

    return(variogram)
  }

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

  # Calculate variogram
  spdf_variogram <- variogram_calculate(spdf = data_scores, target_var = target_var, nlags = nlags, lags = lags)
  base::plot(spdf_variogram)

  if (is.na(nlags)){
    dist_seq <- seq(from = 0, to = max(spdf_variogram$dist), length.out = max(spdf_variogram$dist) - min(spdf_variogram$dist) / (1 + log2(length(spdf_variogram$dist)) * (.8)))
  }else{
    dist_seq <- seq(from = 0, to = max(spdf_variogram$dist), length.out = max(spdf_variogram$dist) - min(spdf_variogram$dist) / (nlags))
  }


    # Get model values based on the provided model_name
    model_values <- get_model_values(dist_seq, nugget = nugget, sill = sill, range_val = range_val, a = a, model_name = model_name)

    # Create a data frame for the model values
    model_df <- data.frame(dist = dist_seq, gamma = model_values)

    ggplot() +
      geom_point(data = spdf_variogram, aes(x = dist, y = gamma)) +
      geom_line(data = model_df, aes(x = dist, y = gamma), color = "red") +
      labs(x = "Distance", y = "Semivariance", title = "Variogram") +
      theme_minimal()
}
