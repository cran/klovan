#' Principal Component Correlation Plot
#'
#' @description
#' This function generates a correlation plot, also known as a "circle" plot, which compares the loadings from one principal component (PC) against another.
#' It visualizes the similarity among original variables and their correlation with each PC, revealing potential clusters. The function also adds annotations for understanding positive and negative values in different quadrants.
#'
#' @param data
#' A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis, see README for details
#' @param PC_1 A string specifying the first PC for comparison, can be chosen from "PC1" to "PC10". For example, "PC1".
#' @param PC_2 A string specifying the second PC for comparison, can be chosen from "PC1" to "PC10". For example, "PC2".
#' @param text_col An R color for the text labels. Defaults to "red".
#'
#' @return A ggplot object representing the correlation plot.
#' @export
#'
#' @examples
#' data("Klovan_Row80")
#' pc_cor_plot(Klovan_Row80, "PC1", "PC2")
#'
pc_cor_plot <- function(data, PC_1, PC_2, text_col = "red"){
  temp <- data %>% dplyr::select(-tidyselect::starts_with("C_")) %>%
    select_if(~{
      is.numeric(.) &&
        length(unique(.)) > 1 &&
        !all(diff(sort(.)) == 1)
    })

  pca <-
    stats::prcomp(temp, center = TRUE, scale. = TRUE)
  PCLoadings <- as.data.frame(pca$rotation[, 1:length(temp)]) %>%
    dplyr::rename() %>%
    tibble::rownames_to_column("VariableName")


  ggplot(PCLoadings, aes_string(PC_1, PC_2, label = "VariableName")) +
    ggtitle("Correlation Plot for Principal Components") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_segment(aes(x = 1, y = 0, xend = -1, yend = 0)) +
    geom_segment(aes(x = 0, y = 1, xend = 0, yend = -1)) +
    geom_point() +
    ggrepel::geom_text_repel(
      hjust = .5,
      vjust = -.5,
      color = text_col,
      size = 3) +
    ggforce::geom_circle(
      aes(x0 = 0, y0 = 0, r = 1),
      color = "gray",
      inherit.aes = FALSE) +
    xlim(-1.2, 1.2) +
    ylim(-1, 1) +
    coord_fixed(ratio = 1)
}
