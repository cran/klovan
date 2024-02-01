#' factor correlation plot
#'
#' @param data A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis, see README for details. It will be that will be converted into
#' a plottable dataframe, see README for details or
#' a plottable data frame created from the `factor_analysis()` function
#'
#' @param FAC_1 the first factor to be compared. A string that can be chosen from FA1:FA3 or FA1:FAnum_fac e.g. "FA1"
#' @param FAC_2 the first factor to be compared. A string that can be chosen from FA1:FA3 or FA1:FAnum_fac e.g. "FA2"
#' @param num_fac a numeric value for how many factors to analyze. Recommended to use 3 and default to 3.
#' @param text_col an R color, the color of the text lables, defaults to "red"
#' @param line_col an R color, the color of the lines, defaults to "lightblue"
#'
#' @return a ggplot object of the correlation plot
#' @export
#'
#' @examples
#' data("Klovan_Row80")
#' fa_plot1 <- factor_cor_plot(Klovan_Row80, "FAC1", "FAC2", 2)
#' fa_plot1
#'
#' factor_cor_plot(Klovan_Row80, "FAC1", "FAC3")
#'
#' fa_plot2 <-factor_cor_plot(factor_analysis(Klovan_Row80), "FAC1", "FAC3", 4)
#' fa_plot2
#'
factor_cor_plot <- function(data, FAC_1, FAC_2, num_fac = 3, text_col = "red", line_col = "lightblue"){
  if (colnames(data)[1] != "VariableName") {
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
  }else{
    RLoadings <- data
  }


  #make correlation plot
  factor_plot <- ggplot(RLoadings, aes_string(FAC_1, FAC_2, label = "VariableName")) +
    ggtitle("Correlation Plot for Rotated Factors") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_segment(aes(x = 1, y = 0, xend = -1, yend = 0)) +
    geom_segment(aes(x = 0, y = 1, xend = 0, yend = -1)) +
    ggrepel::geom_text_repel(
      hjust = .5, vjust = -.5,  color = text_col, size = 3) +
    ggforce::geom_circle(
      aes(x0 = 0, y0 = 0, r = 1),
      color = "gray",
      inherit.aes = FALSE) +
    xlim(-1.2, 1.2) +
    ylim(-1, 1) +
    geom_segment(
      aes(xend = 0, yend = 0),
      arrow = arrow(ends = "first", length = unit(.15, "inches")),
      color = line_col) +
    geom_point() +
    coord_fixed(ratio = 1)

  return(factor_plot)
}
