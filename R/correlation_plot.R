#' Correlation Plot
#'
#' The correlation plot is a summary showing the relationship among variables.The plot below is a 10 x 10 table where each variable is plotted against every other variable.In the top half of the table, the correlation coefficients are displayed.  In the bottom half, the scatter plots are shown along with a regression line.  Down the diagonal axis, the variable histograms are show.
#'
#' @param data A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis, see README for details
#'
#' @return a correlation plot displaying correlation coefficients, the scatter plots with a regression line and, the variable histograms in a 10 x 10 table.
#' @export
#'
#' @examples
#' \donttest{
#' data("Klovan_Row80")
#' correlation_plot(Klovan_Row80)
#' }
correlation_plot <- function(data) {
  data <- data %>% dplyr::select(-tidyselect::starts_with("C_")) %>%
    select_if(~{
      is.numeric(.) &&
        length(unique(.)) > 1 &&
        !all(diff(sort(.)) == 1)
    })
  data %>% GGally::ggpairs(
      aes(alpha = .5),
      #        lower = list(continuous = wrap(lowerFn)),
      upper = list(
        continuous = GGally::wrap(
          "cor",
          size = 3,
          color = "dark red")),
      lower = list(
        continuous = GGally::wrap(
          "smooth",
          alpha = 0.3,
          size = 0.1))) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 7,
        color = "dark blue"),
      axis.text.y = element_text(color = "dark blue", size = 6))
}
