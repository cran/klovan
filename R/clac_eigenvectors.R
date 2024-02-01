#' Calculate Eigenvectors
#'
#' @description
#' This function calculates the Eigenvectors of a given covariance matrix
#' or a klovan dataset. In case of a klovan dataset, it is first converted
#' into a covariance matrix. For further details on klovan datasets, refer to the README.
#'
#' @param data A covariance matrix made with `covar_mtrx()` function, or a
#' A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis.
#'
#' @return A data frame with two columns: "Evalues_COV" and "pc.names1".
#' "Evalues_COV" represents the eigenvectors for each principal component
#' listed in "pc.names1".
#' @export
#'
#' @examples
#' data("Klovan_Row80")
#' calc_eigenvectors(covar_mtrx(Klovan_Row80)) # view eigenvectors
#'
#' @importFrom magrittr %>%
calc_eigenvectors <- function(data){
  if (dim(data)[1] != dim(data)[2] && !all(data < 1)) { #check if imput is cov mtrx
    data <- data %>% dplyr::select(-tidyselect::starts_with("C_")) %>%
      select_if(~{
        is.numeric(.) &&
          length(unique(.)) > 1 &&
          !all(diff(sort(.)) == 1)
      })
    data <- stats::cov(data)
  }
  pc.names <- paste0("PC", 1:(ncol(data) + 1))
  Cov_Mtrx.eigen <- base::eigen(data, symmetric = TRUE, only.values = FALSE)
  Evectors_COV <- as.data.frame(Cov_Mtrx.eigen$vectors)
  base::colnames(Evectors_COV) <- pc.names[1:length(Evectors_COV)]
  return(Evectors_COV)
}
