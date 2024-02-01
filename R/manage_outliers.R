#' Outlier Management
#'
#' @description
#' This function appends a new column to the input data, marking potential outliers.
#' Once identified, these outliers can either be removed or imputed.
#'
#' @param data A dataset of class data.frame.
#' @param property A string representing the property on which the range transformation is based.
#'
#' @return The input dataset, supplemented with a new Boolean column. TRUE signifies a high likelihood of an outlier, while FALSE signifies a low likelihood.
#' @export
#'
#' @examples
#' data("Klovan_2D_all_outlier")
#' manage_outliers(Klovan_2D_all_outlier, "P_Mg")
#'
#' @importFrom stats IQR
#' @importFrom stats quantile
manage_outliers <- function(data, property){
  data <- data %>%
    dplyr::mutate(
      outlier =
        (!!sym(property)) > (quantile(!!sym(property), prob = 0.75, na.rm = TRUE) + ((IQR(!!sym(property))) * 1.5)) |
        (!!sym(property) < (quantile(!!sym(property), prob = .25, na.rm = TRUE) - ((IQR(!!sym(property))) * 1.5)))
    )
  return(data)
}
