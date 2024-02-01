#' Spherical Model
#'
#' Calculate the spherical model based on the given parameters.
#'
#' @param h Distance.
#' @param nugget Nugget effect.
#' @param sill Sill value.
#' @param range Range value.
#'
#' @return Computed value based on the spherical model.
#' @export
spherical_model <- function(h, nugget, sill, range) {
  ifelse(h <= range,
         nugget + (sill - nugget) * (1.5 * h / range - 0.5 * (h / range) ^ 3),
         sill)
}

#' Exponential Model
#'
#' Calculate the exponential model based on the given parameters.
#'
#' @param h Distance.
#' @param nugget Nugget effect.
#' @param sill Sill value.
#' @param range Range value.
#'
#' @return Computed value based on the exponential model.
#' @export
exponential_model <- function(h, nugget, sill, range) {
  nugget + (sill - nugget) * (1 - exp(-h / range))
}

#' Gaussian Model
#'
#' Calculate the Gaussian model based on the given parameters.
#'
#' @param h Distance.
#' @param nugget Nugget effect.
#' @param sill Sill value.
#' @param range Range value.
#'
#' @return Computed value based on the Gaussian model.
#' @export
gaussian_model <- function(h, nugget, sill, range) {
  nugget + (sill - nugget) * (1 - exp(-(h / range)^2))
}

#' Matern Model
#'
#' Calculate the Matern model based on the given parameters.
#'
#' @param h Distance.
#' @param nugget Nugget effect.
#' @param sill Sill value.
#' @param range Range value.
#' @param a Shape parameter.
#'
#' @return Computed value based on the Matern model.
#' @export
matern_model <- function(h, nugget, sill, range, a) {
  gamma_nu = gamma(a)
  double_gamma_nu = gamma(2 * a)

  nugget + (sill - nugget) * ((2 ^ (1 - a)) / double_gamma_nu) * ((h / range) ^ a) * besselK(h / range, a, expon.scaled = FALSE)
}

#' Power Model
#'
#' Calculate the power model based on the given parameters.
#'
#' @param h Distance.
#' @param nugget Nugget effect.
#' @param sill Sill value.
#' @param a Power parameter.
#'
#' @return Computed value based on the power model.
#' @export
power_model <- function(h, nugget, sill, a) {
  nugget + sill * h^a
}

#' Quadratic Exponential Model
#'
#' Calculate the quadratic exponential model based on the given parameters.
#'
#' @param h Distance.
#' @param nugget Nugget effect.
#' @param sill Sill value.
#' @param range Range value.
#' @param a Additional parameter.
#'
#' @return Computed value based on the quadratic exponential model.
#' @export
quadratic_exponential_model <- function(h, nugget, sill, range, a) {
  nugget + a * h^2 + sill * (1 - exp(-h / range))
}

#' Cardinal Sine Model
#'
#' Calculate the cardinal sine model based on the given parameters.
#'
#' @param h Distance.
#' @param sill Sill value.
#'
#' @return Computed value based on the cardinal sine model.
#' @export
cardinal_sine_model <- function(h, sill) {
  sill * sin(h) / h
}

#' Gamma Model
#'
#' Calculate the gamma model based on the given parameters.
#'
#' @param h Distance.
#' @param nugget Nugget effect.
#' @param sill Sill value.
#' @param range Range value.
#' @param a Additional parameter.
#'
#' @return Computed value based on the gamma model.
#' @export
gamma_model <- function(h, nugget, sill, range, a) {
  nugget + sill * (1 - (1 + (h / range)^a)^(-1/a))
}

#' Cauchy Model
#'
#' Calculate the Cauchy model based on the given parameters.
#'
#' @param h Distance.
#' @param nugget Nugget effect.
#' @param sill Sill value.
#' @param range Range value.
#' @param a Additional parameter.
#'
#' @return Computed value based on the Cauchy model.
#' @export
cauchy_model <- function(h, nugget, sill, range, a) {
  nugget + sill * (1 / (1 + (h / range)^a))
}

#' Stable Model
#'
#' Calculate the stable model based on the given parameters.
#'
#' @param h Distance.
#' @param nugget Nugget effect.
#' @param sill Sill value.
#' @param range Range value.
#' @param a Additional parameter.
#'
#' @return Computed value based on the stable model.
#' @export
stable_model <- function(h, nugget, sill, range, a) {
  nugget + sill * exp(-1 * (h / range)^a)
}

#' Order-1 G.C. (General Covariance) Model
#'
#' Calculate the order-1 G.C. model based on the given parameters.
#'
#' @param h Distance.
#' @param nugget Nugget effect.
#' @param sill Sill value.
#' @param range Range value.
#'
#' @return Computed value based on the order-1 G.C. model.
#' @export
order_1_gc_model <- function(h, nugget, sill, range) {
  nugget + sill * exp(-h / range)
}

#' Cosinus Model
#'
#' Calculate the cosinus model based on the given parameters.
#'
#' @param h Distance.
#' @param sill Sill value.
#' @param a Additional parameter.
#'
#' @return Computed value based on the cosinus model.
#' @export
cosinus_model <- function(h, sill, a) {
  sill * cos(a * h)
}

#' Triangle Model
#'
#' Calculate the triangle model based on the given parameters.
#'
#' @param h Distance.
#' @param nugget Nugget effect.
#' @param sill Sill value.
#' @param range Range value.
#'
#' @return Computed value based on the triangle model.
#' @export
triangle_model <- function(h, nugget, sill, range) {
  ifelse(h <= range,
         nugget + sill * (1 - h / range),
         sill)
}

get_model <- function(name){
  # Define a list of models to test
  models_to_test <- list("Sph1" = spherical_model, "Exp1" = exponential_model, "Gau1" = gaussian_model, "Mat1" = matern_model, "Pow1" = power_model, "Quad1" = quadratic_exponential_model, "Card1" = cardinal_sine_model, "Gam1" = gamma_model, "Cau1" = cauchy_model, "Sta1" = stable_model, "Ord1" = order_1_gc_model, "Tri1" = triangle_model, "Cos1" = cosinus_model)

  return(models_to_test[[name]])
}


#' Print Model Names
#'
#' This function prints the names of the predefined model functions.
#'
#' @return NULL (This function is used for printing the model names only.)
#'
#' @examples
#' print_model_names()
#'
#' @importFrom utils head
#'
#' @export
print_model_names <- function() {
  message(list("Sph1" = "spherical_model", "Exp1" = "exponential_model", "Gau1" = "gaussian_model", "Mat1" = "matern_model", "Pow1" = "power_model", "Quad1" = "quadratic_exponential_model", "Card1" = "cardinal_sine_model", "Gam1" = "gamma_model", "Cau1" = "cauchy_model", "Sta1" = "stable_model", "Ord1" = "order_1_gc_model", "Tri1" = "triangle_model", "Cos1" = "cosinus_model"))
}



