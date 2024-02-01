#' Map the Factor Scores
#'
#' @description
#' This function creates a faceted plot of each rotated factor score,
#' which could be interpreted as the elements of a "phantom" variable.
#' This function aids in defining the relationship between the phantom variables
#' and the known ore body by producing a contoured map for each variable.
#'
#' @param Interp_Data A plottable data frame produced by the `inv_dis_wt()` or `kriging()` functions.
#' @param overlay A Boolean input. If TRUE, the plot will overlay isolines; if FALSE, it will not.
#' @param data a klovan dataset (transformed, untransformed, outlier, etc), see README for details.
#' @param FA_colors A named vector of colors for different factors. Defaults
#' to a set color palette.
#'
#' @import ggplot2
#'
#' @return A ggplot object representing the Factor Scores plot.
#' @export
#'
#' @examples
#' \donttest{
#' data("Klovan_Row80")
#' factor_plot1 <- factor_score_plot(inv_dis_wt(Klovan_Row80), TRUE, data = Klovan_Row80)
#' factor_plot1
#'
#' your_interp_data_IDW <- inv_dis_wt(Klovan_Row80, 3)
#' factor_score_plot(your_interp_data_IDW, FALSE, data = Klovan_Row80)
#' }
#'
factor_score_plot <- function(Interp_Data, overlay, data, FA_colors = c(FA1 = "black", FA2 = "blue", FA3 = "darkred", FA4 = "green", FA5 = "purple", FA6 = "orange", FA7 = "yellow", FA8 = "pink", FA9 = "cyan", FA10 = "magenta")){
  C_X <- C_Y <- value <- FA <- NULL
  if (overlay == FALSE) {
    factor_plot <- ggplot(
      Interp_Data,
      aes(C_X, C_Y, z = value, color = FA)) +
      ggplot2::geom_contour_filled(
        binwidth = .25,
        size = .2,
        alpha = .5) +
      geom_point(
        data = data,
        mapping = aes(x = C_X, y = C_Y),
        inherit.aes = FALSE) +
      # metR::geom_text_contour(
      #   label.size = .5,
      #   stroke = 0) +
      facet_wrap(vars(FA), ncol = 2) +
      # geom_ellipse(
      #   aes(x0 = 3900, y0 = 1700, a = 600, b = 400, angle = pi/2.5),
      #   color = "white") +
      scale_colour_manual(values = FA_colors, guide = "none" )
  }else{
    factor_plot <- ggplot(Interp_Data, aes(C_X, C_Y, z = value, color = FA)) +
      #  geom_contour(binwidth = .5, size =.4) +
      geom_contour_filled(binwidth = .5, size = .4, alpha = .2) +
      metR::geom_text_contour(
        min.size = 5,
        binwidth = .5,
        label.size = .5,
        nudge_x = 0,
        nudge_y = 0) +
      # geom_label_contour(
      #   min.size = 5,
      #   label.size = .5,
      #   nudge_x = 0,
      #   nudge_y = 0) +
      geom_point(
        data = data,
        mapping = aes(x = C_X, y = C_Y),
        inherit.aes = FALSE) +
      # geom_ellipse(
      #   aes(x0 = 3900, y0 = 1700, a = 600, b = 400, angle = pi/2.5),
      #   color = "white") +
      # geom_circle(
      #   aes(x = NULL, y = NULL, x0 = 3300, y0 = 3500, r = 400),
      #   color = "white",
      #   inherit.aes = FALSE  ) +
      coord_fixed(ratio = 1) +
      scale_colour_manual(values = FA_colors, guide = "none" )

  }
  return(factor_plot)
}
