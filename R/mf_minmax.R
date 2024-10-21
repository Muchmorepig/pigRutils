#' Normalize Values to [0, 1] Range
#'
#' This function normalizes the input values to a range between 0 and 1 using min-max normalization.
#' It can handle both vectors and matrices. For matrices, normalization is applied row-wise.
#'
#' @param x A numeric vector or matrix to be normalized.
#'
#' @return 
#' If `x` is a vector, returns a normalized vector.
#' If `x` is a matrix, returns a matrix with each row normalized.
#'
#' @details
#' The function uses the formula: (x - min(x)) / (max(x) - min(x))
#' If all values in a vector or a row of a matrix are the same (i.e., min = max),
#' a warning is issued and the original values are returned for that vector or row.
#'
#' @examples
#' # Vector example
#' vec <- c(1, 2, 6, 9, 4, 4, 7.2)
#' mf_minmax(vec)
#'
#' # Matrix example
#' mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
#' mf_minmax(mat)
#'
#' @export

mf_minmax <- function(x) {
  minmax <- function(x) {
    min_val <- min(x, na.rm = TRUE)
    max_val <- max(x, na.rm = TRUE)
    
    if (min_val == max_val) {
      warning("min = max, returning original values")
      return(x)
    }
    
    (x - min_val) / (max_val - min_val)
  }
  
  if (is.vector(x)) {
    return(minmax(x))
  } else if (is.matrix(x)) {
    return(t(apply(x, 1, minmax)))
  } else {
    stop("Input must be a vector or matrix")
  }
}

#mf_minmax <- function(x) {
#  minmax <- function(x) {
#    min_val <- min(x, na.rm = TRUE)
#    max_val <- max(x, na.rm = TRUE)
#    if (min_val == max_val) {
#      stop("min = max")
#    }
#    (x - min_val) / (max_val - min_val)
#  }
#  t(apply(x, 1, minmax))
#}
