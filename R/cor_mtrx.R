#' correlation matrix
#'
#' @param data A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis.
#'
#' @return a correlation matrix as matrix object
#' @export
#'
#' @examples
#' data("Klovan_Row80")
#' cor_mtrx(range_transform(Klovan_Row80)) # view correlation matrix
#' corMtrx <- cor_mtrx(Klovan_Row80) # save correlation matrix as object
cor_mtrx <- function(data){
  Cor_Mtrx <- stats::cor(data)
  return(Cor_Mtrx)
}
