#' Perform Factor Analysis
#'
#' @description
#' This function performs a Factor Analysis on a provided dataset using the "Varimax" orthogonal rotation method.
#' It also calculates the factor scores for each factor.
#'
#' @param data A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis, see README for details
#'
#' @param num_fac A numeric value indicating the number of factors to analyze.
#' It's recommended to use 3, which is also the default value.
#'
#' @return A data frame containing the calculated factors.
#' @export
#'
#' @examples
#' data("Klovan_Row80")
#' factor_analysis(Klovan_Row80)
#' factor_analysis(Klovan_Row80, 3)
#'
#' @importFrom stats varimax
factor_analysis <- function(data, num_fac = 3){
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

  RLoadings <- as.data.frame(rotatedLoadings[, 1:num_fac]) # Select the number of factors you want from 1:num_fac
  colnames(RLoadings) <- paste0("FAC", 1:num_fac)  # Generate column names dynamically
  RLoadings <- data.frame(VariableName = rownames(RLoadings), round(RLoadings[, 1:num_fac], 3))  # Add VariableName column and round the values
  return(RLoadings)
}
