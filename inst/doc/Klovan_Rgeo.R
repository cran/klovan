## ----eval = FALSE-------------------------------------------------------------
#  # install.packages("Klovan_0.0.9.tar.gz", repos = NULL, type = "source")
#  # library(klovan)

## ----eval = FALSE-------------------------------------------------------------
#  # install.packages("klovan")
#  # library(klovan)

## ----warning = FALSE, eval = FALSE--------------------------------------------
#  #library(RGeostats)
#  #data("Klovan_Row80", package = "klovan")
#  Klovan_Row80 <- load(file = "~/CSE_MSE_RXF131/cradle-members/sdle/jeg165/git/klovan/packages/Klovan0.0.9/data/Klovan_Row80.rda")

## ----warning = FALSE, eval = FALSE--------------------------------------------
#  # Building a database based on RC1 factor
#  db <- Rgeo_database(Klovan_Row80, 3, "RC1")

## ----warning = FALSE, eval = FALSE--------------------------------------------
#  print(db)

## ----warning = FALSE, eval = FALSE--------------------------------------------
#  # Construct and plot the experimental variogram
#  Rgeo_vario_construct_plot(db, 3, "RC1", lag = 500)

## ----warning = FALSE, eval = FALSE--------------------------------------------
#  # Fit the variogram model based on experimental variogram
#  model <- Rgeo_vario_model(db, 3, "RC1", lag = 500, model = 13)

## ----warning = FALSE, eval = FALSE--------------------------------------------
#  print(model)

## ----warning = FALSE, eval = FALSE--------------------------------------------
#  krig <- Rgeo_kriging(db, model)

## ----warning = FALSE, eval = FALSE--------------------------------------------
#  print(krig)

## ----warning = FALSE, eval = FALSE, eval = FALSE------------------------------
#  # Plot the kriging estimation results
#  Rgeo_kriging_plot(krig, db, "RC1")

## ----out.width="675px", dpi=1000, echo=FALSE, fig.cap="Krige Plot"------------
knitr::include_graphics("Rgeo_krig.png")

