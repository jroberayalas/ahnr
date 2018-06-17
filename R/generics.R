#' Summary Artificial Hydrocarbon Network
#'
#' @description Summary method for objects of class \code{ahn}.
#'
#' @param x an object of class "\code{ahn}" produced from the \link{fit} function.
#' @param ... further arguments passed to or from other methods.
#'
#' @return summary description of the AHN.
#' @export
#'
#' @examples
#' \dontrun{
#' # Create data
#' x <- 2 * runif(1000) - 1;
#' x <- sort(x)
#'
#' y <- (x < 0.1) * (0.05 * runif(100) + atan(pi*x)) +
#'     (x >= 0.1 & x < 0.6) * (0.05 * runif(1000) + sin(pi*x)) +
#'     (x >= 0.6) * (0.05 * runif(1000) + cos(pi*x))
#'
#' # Create Sigma list
#' Sigma <- list(X = data.frame(x = x), Y = data.frame(y = y))
#'
#' # Train AHN
#' ahn <- fit(Sigma, 5, 0.01, 500)
#'
#' # Summary AHN
#' summary(ahn)
#' }
#'
#summary <- function(x, ...) {
#    UseMethod("summary", x)
#}


#' predict
#'
#' @description Function to simulate a trained Artificial Hydrocarbon Network.
#'
#' @param x an object of class "\code{ahn}" produced from the \link{fit} function.
#' @param ... further arguments passed to or from other methods.
#'
#' @return predicted output values for inputs in \code{new_data}.
#' @export
#'
#' @examples
#' \dontrun{
#' # Create data
#' x <- 2 * runif(1000) - 1;
#' x <- sort(x)
#'
#' y <- (x < 0.1) * (0.05 * runif(100) + atan(pi*x)) +
#'     (x >= 0.1 & x < 0.6) * (0.05 * runif(1000) + sin(pi*x)) +
#'     (x >= 0.6) * (0.05 * runif(1000) + cos(pi*x))
#'
#' # Create Sigma list
#' Sigma <- list(X = data.frame(x = x), Y = data.frame(y = y))
#'
#' # Train AHN
#' ahn <- fit(Sigma, 5, 0.01, 500)
#'
#' # Test AHN
#' X <- data.frame(x = x)
#' ysim <- predict(ahn, X)
#' }
#'
#predict <- function(x, ...) {
#    UseMethod("predict", x)
#}


#' Plot Artificial Hydrocarbon Network
#'
#' @description Plot method for objects of class \code{ahn}.
#'
#' @param x an object of class "\code{ahn}" produced from the \link{fit} function.
#' @param ... further arguments passed to visNetwork functions.
#'
#' @return dynamic visualization of the AHN.
#' @export
#'
#' @examples
#' \dontrun{
#' # Create data
#' x <- 2 * runif(1000) - 1;
#' x <- sort(x)
#'
#' y <- (x < 0.1) * (0.05 * runif(100) + atan(pi*x)) +
#'     (x >= 0.1 & x < 0.6) * (0.05 * runif(1000) + sin(pi*x)) +
#'     (x >= 0.6) * (0.05 * runif(1000) + cos(pi*x))
#'
#' # Create Sigma list
#' Sigma <- list(X = data.frame(x = x), Y = data.frame(y = y))
#'
#' # Train AHN
#' ahn <- fit(Sigma, 5, 0.01, 500)
#'
#' # Plot AHN
#' plot(ahn)
#' }
#'
#plot <- function(x, ...) {
#    UseMethod("plot", x)
#}
