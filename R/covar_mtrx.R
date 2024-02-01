#' Create Co-variance Matrix
#'
#' @description
#' This function creates a non-normalized co-variance matrix
#' from the given klovan dataset. For further details on klovan datasets,
#' refer to the README.
#'
#' @param data A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis.
#'
#' @return A non-normalized co-variance matrix of the klovan data.
#' @export
#'
#' @examples
#' data("Klovan_Row80")
#' covar_mtrx(Klovan_Row80) # view co-variance matrix
covar_mtrx <- function(data){
  Cov_Mtrx <- stats::cov(data %>% dplyr::select(-tidyselect::starts_with("C_")) %>%
                           select_if(~{
                             is.numeric(.) &&
                               length(unique(.)) > 1 &&
                               !all(diff(sort(.)) == 1)
                           }))
  return(Cov_Mtrx)
}
