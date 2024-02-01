## ----eval = FALSE-------------------------------------------------------------
#  # install.packages("Klovan_0.0.9.tar.gz", repos = NULL, type = "source")
#  # library(klovan)

## ----include=FALSE------------------------------------------------------------
# install.packages("klovan")
# library(klovan)
library(ggplot2)

## ----warning=FALSE, results='hide'--------------------------------------------
#loading data
data("Klovan_Row80", package = "klovan")
data("Klovan_2D_all_outlier", package = "klovan")

#apply a range transform to your data 
T_klovan <- klovan::range_transform(Klovan_Row80)

## ----warning = FALSE----------------------------------------------------------
#build a correlation matrix 
cov_mtrx <- klovan::covar_mtrx(T_klovan)
cov_mtrx

#calulate Eiegn values
klovan::calc_eigenvalues(cov_mtrx)

#calulate Eiegn vectors
klovan::calc_eigenvectors(cov_mtrx)

## ----warning = FALSE----------------------------------------------------------
eigen_data <- klovan::eigen_contribution(T_klovan)
eigen_data

## ----warning = FALSE, eval = FALSE--------------------------------------------
#  klovan::scree_plot(eigen_data)

## ----out.width="675px", dpi=1000, echo=FALSE, fig.cap="Scree Plot"------------
knitr::include_graphics("scree_plot.png")

## ----warning = FALSE, eval = FALSE--------------------------------------------
#  klovan::scree_plot(eigen_data, bar_fill = "green", outline = "darkgreen", eigen_line = "lightblue")

## ----out.width="675px", dpi=1000, echo=FALSE, fig.cap="Scree Plot"------------
knitr::include_graphics("scree_plot_color.png")

## ----warning = FALSE----------------------------------------------------------
#make a correlation plot 
klovan::cor_mtrx(Klovan_Row80)

## ----warning = FALSE----------------------------------------------------------
klovan::pc_cor_plot(Klovan_Row80, "PC1", "PC2")
#see function decimation for more information on how to interpret this plot

## ----warning = FALSE----------------------------------------------------------
#factor analysis 
klovan::factor_analysis(Klovan_Row80)

## ----warning = FALSE----------------------------------------------------------
#make correlation plot using factor data
klovan::factor_cor_plot(klovan::factor_analysis(Klovan_Row80), "FAC1", "FAC2")
#customize color choices 
klovan::factor_cor_plot(Klovan_Row80, "FAC1", "FAC3", text_col = "pink", line_col = "red")

## ----warning = FALSE----------------------------------------------------------
#use inverse distance weighted method for interpolation
inv_dis_data <- klovan::inv_dis_wt(Klovan_Row80, 3)
summary(inv_dis_data) #view data summary

## ----warning = FALSE, eval = FALSE--------------------------------------------
#  klovan::factor_score_plot(a, FALSE, data = Klovan_Row80) + ggforce::geom_ellipse(
#      aes(x0 = 3900, y0 = 1700, a = 600, b = 400, angle = pi/2.5),
#      color = "white")

## ----out.width="675px", dpi=1000, echo=FALSE, fig.cap="Inverse Distance Weighting Plot"----
knitr::include_graphics("invs_dis3.png")

## ----warning = FALSE, eval = FALSE--------------------------------------------
#  klovan::factor_score_plot(inv_dis_data, TRUE, data = Klovan_Row80) + ggforce::geom_ellipse(
#      aes(x0 = 3900, y0 = 1700, a = 600, b = 400, angle = pi/2.5),
#      color = "white") +
#    ggforce::geom_circle(
#      aes(x = NULL, y = NULL, x0 = 3300, y0 = 3500, r = 400),
#    color = "white",
#    inherit.aes = FALSE)

## ----out.width="675px", dpi=1000, echo=FALSE, fig.cap="Inverse Distance Weighting Plot"----
knitr::include_graphics("invs_dis.png")

## ----warning = FALSE----------------------------------------------------------
#plot variogram for use in kriging
klovan::vario_plot(Klovan_Row80, factor = 1, nugget = .214, nlags = 10, sill = 7.64507, range_val = 6271.83, model_name = "Gau1")

## ----warning = FALSE, eval = FALSE--------------------------------------------
#  #use kriging method for interpolation and
#  #plot with factors overlapped and separated
#  krig_data <- klovan::kriging.auto(Klovan_Row80, 3) #customize available for nugget, psill, range, and model see function documentation for more details
#  summary(krig_data) #view data summary

## ----warning = FALSE, eval = FALSE--------------------------------------------
#  klovan::factor_score_plot(krig_data, TRUE, data = Klovan_Row80) + ggforce::geom_ellipse(
#      aes(x0 = 3900, y0 = 1700, a = 600, b = 400, angle = pi/2.5),
#      color = "white") +
#    ggforce::geom_circle(
#      aes(x = NULL, y = NULL, x0 = 3300, y0 = 3500, r = 400),
#    color = "white",
#    inherit.aes = FALSE)

## ----out.width="675px", dpi=1000, echo=FALSE, fig.cap="Krige Plot"------------
knitr::include_graphics("krig_fin1.png")

