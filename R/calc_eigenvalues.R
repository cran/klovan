#' Calculate Eigenvalues
#'
#' @description
#' This function calculates the eigenvalues of a given covariance matrix
#' or a klovan dataset. In case of a klovan dataset, it is first converted
#' into a covariance matrix. For further details on klovan datasets, refer to the README.
#'
#' @param data A covariance matrix made with `covar_mtrx()` function, or a
#' A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis.
#'
#' @return A data frame with two columns: "Evalues_COV" and "pc.names1".
#' "Evalues_COV" represents the eigenvalues for each principal component
#' listed in "pc.names1".
#' @export
#'
#' @examples
#' data("Klovan_Row80")
#' calc_eigenvalues(covar_mtrx(Klovan_Row80)) # view eigenvalues
#'
#' @importFrom magrittr %>%
calc_eigenvalues <- function(data) {
  if (dim(data)[1] != dim(data)[2] && !all(data < 1)) { #check if imput is cov mtrx
    data <- data %>% dplyr::select(-tidyselect::starts_with("C_")) %>%
      select_if(~{
        is.numeric(.) &&
          length(unique(.)) > 1 &&
          !all(diff(sort(.)) == 1)
      })
    data <- stats::cov(data)
  }
  pcnames1.vec <- paste0("PC", 1:(ncol(data) + 1))
  Cov_Mtrx.eigen <- base::eigen(data, symmetric = TRUE, only.values = FALSE)
  Evalues_COV <- as.data.frame(Cov_Mtrx.eigen$values)
  Evalues_COV <- base::data.frame(Evalues_COV, pc.names1 = base::factor(pcnames1.vec[1:length(Evalues_COV)], levels = pcnames1.vec))
  return(Evalues_COV)
}
