#' scree plot
#'
#' @param EigenPlot A dataset of class data.frame. The data should contain 'C_X' and 'C_Y' columns
#' representing the x and y coordinates of the data points and excludes any rank,
#' ID, or column not for analysis, see README for details or,
#'                  a covariance matrix that will be converted into a plottable data frame or,
#'                  a plottable data frame created by the `eigen_contribution()` function
#'
#' @param bar_fill an R color, The fill color for the bars, defaults to "lightblue"
#' @param outline an R color, the outline color of the bars, defaults to "darkblue"
#' @param eigen_line an R color, the color of the eigenvalues line, defaults to "red"
#' @param cum_eigen_line an R color, the color of the cummulative eigenvalues line, defaults to "blue"
#'
#' @return a ggplot object of the scree plot
#' @export
#'
#' @examples
#' data("Klovan_Row80")
#' scree_plot(eigen_contribution(covar_mtrx(Klovan_Row80)))
#'
#'
#' Scree1 <- scree_plot(Klovan_Row80)
#' Scree1
#'
#' your_eigen_data1 <- eigen_contribution(Klovan_Row80)
#' scree_plot(your_eigen_data1)
#'
scree_plot <- function(EigenPlot, bar_fill = "lightblue", outline = "darkblue", eigen_line = "red", cum_eigen_line = "blue") {
 CumSum <- pc.names <- EigenValues <- NULL
   if ((dim(EigenPlot)[1] != dim(EigenPlot)[2] && !all(EigenPlot < 1)) & colnames(EigenPlot)[2] != "CumSum") {
    EigenPlot <- stats::cov(EigenPlot %>% dplyr::select(-rank) %>% dplyr::select(-tidyselect::starts_with("C_")) %>%
                                          select_if(~{
                                            is.numeric(.) &&
                                              length(unique(.)) > 1 &&
                                              !all(diff(sort(.)) == 1)
                                          }))
  }
  if (colnames(EigenPlot)[2] != "CumSum") {
    pcnames1.vec <- paste0("PC", 0:(ncol(EigenPlot) + 1))
    Cov_Mtrx.eigen <- eigen(EigenPlot, symmetric = TRUE, only.values = FALSE)
    Evalues_COV <- as.data.frame(Cov_Mtrx.eigen$values)
    EigenPlot <-
      data.frame(
        EigenValues = c(NA, Cov_Mtrx.eigen$values/sum(Cov_Mtrx.eigen$values)),
        CumSum = c(0, cumsum(Cov_Mtrx.eigen$values/sum(Cov_Mtrx.eigen$values))),
        CumSumPct = c(0, cumsum(Cov_Mtrx.eigen$values/sum(Cov_Mtrx.eigen$values))) * 100,
        pc.names = factor(pcnames1.vec[1:(length(Evalues_COV[[1]])+1)], levels = pcnames1.vec))
  }

  ggplot(EigenPlot[2:length(EigenPlot$EigenValues),], aes(as.numeric(pc.names), EigenValues)) +
    ggtitle("Scree Plot") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_bar( fill = bar_fill, stat = "identity", color = outline) +
    geom_line(aes(color = "Eigenvalues"), size = 1) +
    geom_line(aes(y = CumSum, color = "Cummulative Eigenvalues"), size = 1) +
    geom_point() +
    geom_point(aes(y = CumSum)) +
    geom_text(aes(label = round(EigenValues, 2)), stat = "identity", vjust = 1.5, colour = "black") +
    scale_x_continuous(
      name = "Principal Components",
      breaks = 1:nrow(EigenPlot),
      labels = EigenPlot$pc.names) +
    scale_color_manual(
      name = "Legend",
      values = c("Eigenvalues" = eigen_line, "Cummulative Eigenvalues" = cum_eigen_line))
}

