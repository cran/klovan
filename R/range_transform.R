#' range transform
#'
#' @description
#' Normalize the data using a 'Range' transform .
#' In the returned data table, note that in each column of the normalized Data Table, the variables will range from 0 to 1.
#'
#' @param data A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis, see README for details
#'
#' @return a range transformed version of a klovan dataset.
#' @export
#'
#' @examples
#' data("Klovan_Row80")
#' T_Klovan <- range_transform(Klovan_Row80)
#'
range_transform <- function(data){
  temp <- data %>% dplyr::select(-tidyselect::starts_with("C_")) %>%
    select_if(~{
      is.numeric(.) &&
        length(unique(.)) > 1 &&
        !all(diff(sort(.)) == 1)
    })

  T_data <- fields::transformx(temp, scale.type = "range")
  T_data <- as.data.frame(T_data)
  par_col_names <- colnames(temp)
  T_data <- cbind(data[,-match(par_col_names, colnames(data))], T_data)
  return(T_data)
}
