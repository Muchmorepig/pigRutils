#' Calculate PDF Dimensions for Multiple Plots
#'
#' This function calculates the appropriate width and height for a PDF document
#' that will contain multiple plots arranged in a grid layout.
#'
#' @param n_plots Integer. The total number of plots to be arranged.
#' @param ncol Integer. Number of columns in the plot grid layout. Default is 2.
#' @param square_size Numeric. The size (in inches) of each plot square. Default is 8.
#' @param margin Numeric. The margin factor as a proportion of square size. Default is 0.15.
#'
#' @return A list containing:
#'   \itemize{
#'     \item width: The calculated width of the PDF document in inches
#'     \item height: The calculated height of the PDF document in inches
#'   }
#'
#' @examples
#' # Calculate dimensions for 5 plots
#' dims <- calculate_pdf_dimensions(5)
#'
#' # Calculate dimensions with custom settings
#' dims <- calculate_pdf_dimensions(
#'   n_plots = 6,
#'   ncol = 3,
#'   square_size = 6,
#'   margin = 0.2
#' )
#'
#' @export
calculate_pdf_dimensions <- function(n_plots,
                                     ncol = 2,
                                     square_size = 8,
                                     margin = 0.15) {
  # Input validation
  if (!is.numeric(n_plots) || n_plots < 1 || n_plots != round(n_plots)) {
    stop("'n_plots' must be a positive integer")
  }
  if (!is.numeric(ncol) || ncol < 1 || ncol != round(ncol)) {
    stop("'ncol' must be a positive integer")
  }
  if (!is.numeric(square_size) || square_size <= 0) {
    stop("'square_size' must be a positive number")
  }
  if (!is.numeric(margin) || margin < 0) {
    stop("'margin' must be a non-negative number")
  }

  # Calculate number of rows needed
  nrow <- ceiling(n_plots / ncol)

  # Calculate dimensions with margins
  width <- square_size * ncol * (1 + margin)
  height <- square_size * nrow * (1 + margin)

  # Return dimensions as a list
  return(list(
    width = width,
    height = height
  ))
}
