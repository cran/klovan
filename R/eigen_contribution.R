#' Calculate Eigen Contribution
#'
#' @description
#' This function calculates the sum of all the eigenvalues from a provided covariance matrix or klovan dataset.
#' Each eigenvalue is divided by the sum of the eigenvalues to determine its proportional contribution.
#' This yields the percent contribution of each eigenvalue and provides an understanding of the proportion of total variance explained by each eigenvalue.
#'
#' @param Cov_Mtrx A covariance matrix used to compute the eigenvalues or
#' A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis, see README for details
#'
#' @return A data frame with columns: "EigenValues", "CumSum", "CumSumPct", "pc.names". Where:
#'         - "EigenValues": The eigenvalues
#'         - "CumSum": The cumulative sum of the eigenvalues
#'         - "CumSumPct": The proportional contribution of each eigenvalue
#'         - "pc.names": The principal component names
#' @export
#'
#' @examples
#' data("Klovan_Row80")
#' your_cov_Mtrx <- covar_mtrx(Klovan_Row80) # example covariance matrix
#' eigen_contribution(Klovan_Row80) # view the data frame
#' eigen_contribution(your_cov_Mtrx) # view the data frame
#' eigen_contribution(covar_mtrx(Klovan_Row80)) # view the data frame
#'
#' @importFrom tidyselect starts_with
eigen_contribution <- function(Cov_Mtrx){
  if (dim(Cov_Mtrx)[1] != dim(Cov_Mtrx)[2] && !all(Cov_Mtrx < 1)) { #check if imput is cov mtrx
    Cov_Mtrx <- stats::cov(Cov_Mtrx %>% dplyr::select(-tidyselect::starts_with("C_")) %>%
                             select_if(~{
                               is.numeric(.) &&
                                 length(unique(.)) > 1 &&
                                 !all(diff(sort(.)) == 1)
                             }))
  }
  pcnames1.vec <- paste0("PC", 0:(ncol(Cov_Mtrx) + 1))
  Cov_Mtrx.eigen <- eigen(Cov_Mtrx, symmetric = TRUE, only.values = FALSE)
  Evalues_COV <- as.data.frame(Cov_Mtrx.eigen$values)
  EigenFrame <-
    data.frame(
      EigenValues = c(NA, Cov_Mtrx.eigen$values/sum(Cov_Mtrx.eigen$values)),
      CumSum = c(0, cumsum(Cov_Mtrx.eigen$values/sum(Cov_Mtrx.eigen$values))),
      CumSumPct = c(0, cumsum(Cov_Mtrx.eigen$values/sum(Cov_Mtrx.eigen$values))) * 100,
      pc.names = factor(pcnames1.vec[1:(length(Evalues_COV[[1]]) + 1)], levels = pcnames1.vec))
  return(EigenFrame)
}

