#' Calculate Mean Values by Group for Matrix Data
#'
#' @description
#' This function calculates the mean values of a matrix for each specified group.
#' It is particularly useful for summarizing high-dimensional data (e.g., gene expression data)
#' across different experimental conditions or groups.
#'
#' @param mat A numeric matrix where rows represent features (e.g., genes) and columns represent samples.
#' @param group A vector specifying the group for each sample. Must be the same length as the number of columns in `mat`.
#' @param samples A vector of sample names. Default is `colnames(mat)`. If provided, must match the length of `group`.
#' @param na.rm Logical. Should NA values be removed when calculating means? Default is FALSE.
#' @param verbose Logical. If TRUE, print progress messages. Default is TRUE.
#'
#' @return A data frame where:
#'   \itemize{
#'     \item Each column represents a unique group
#'     \item Each row represents a feature (e.g., gene)
#'     \item Values are the mean of the features for each group
#'   }
#'
#' @examples
#' # Create a sample matrix (20 genes, 10 samples)
#' set.seed(123)
#' mat <- matrix(rnorm(200), nrow = 20, ncol = 10)
#' rownames(mat) <- paste0("Gene", 1:20)
#' colnames(mat) <- paste0("Sample", 1:10)
#'
#' # Define groups
#' group <- rep(c("Control", "Treatment"), each = 5)
#'
#' # Calculate mean values
#' result <- mf_mean(mat, group)
#' print(head(result))
#'
#' @export
#'
#' @importFrom stats setNames
#' @importFrom methods is

mf_mean <- function(mat, group, samples = colnames(mat), verbose = FALSE) {
  # Input validation
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("'mat' must be a numeric matrix")
  }
  if (length(samples) != length(group)) {
    stop("Length of 'group' must equal the number of columns in 'mat'")
  }
  if (length(samples) != ncol(mat)) {
    stop("Length of 'samples' must equal the number of columns in 'mat'")
  }

  # Assign group names to samples
  names(samples) <- group

  # Get unique group names
  unique_groups <- unique(group)

  # Calculate mean for each group
  group_means <- vapply(unique_groups, FUN = function(x) {
    if (verbose) message("Processing group: ", x)
    rowMeans(mat[, names(samples) == x, drop = FALSE])
  }, FUN.VALUE = numeric(nrow(mat)))

  # Create and return the result data frame
  result_df <- as.data.frame(group_means)
  rownames(result_df) <- rownames(mat)
  return(result_df)
}

