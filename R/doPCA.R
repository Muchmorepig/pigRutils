#' This function performs PCA on input data and extracts specified principal components
#' while preserving group information.
#'
#' @param indata A numeric matrix or data frame containing the input data 
#'               (typically normalized gene expression data)
#' @param coldata A data frame containing sample metadata/annotations
#' @param intgroup A character string specifying the column in coldata to use for grouping
#'                 (default is "groups")
#' @param ntop Number of top variable features to use for PCA calculation. 
#'             Defaults to the minimum of 2000 or 10% of total features.
#' @param pcsToUse A numeric vector specifying which principal components to extract 
#'                 (default is first 3 components)
#'
#' @return A data frame containing:
#'   - Principal component coordinates for selected components
#'   - Group information
#'   - Sample names
#'   With an attribute "percentVar" showing the percentage of variance explained 
#'   by the selected principal components
#'
#' @importFrom stats prcomp var
#'
#' @examples
#' # Assuming you have expression data and sample metadata
#' expr_data <- matrix(rnorm(1000), ncol=10)
#' sample_info <- data.frame(groups = c(rep("A", 5), rep("B", 5)))
#' pca_result <- doPCA(expr_data, sample_info)
#'
#' @export
doPCA <- function(
    indata,
    coldata,
    intgroup = "groups",
    ntop = ifelse(0.1 * nrow(indata) > 2000, 0.1 * nrow(indata), 2000),
    pcsToUse = 1:3
) {
    # Validate inputs
    if (!is.matrix(indata) && !is.data.frame(indata)) {
        stop("'indata' must be a matrix or data frame")
    }
    
    if (!is.data.frame(coldata)) {
        stop("'coldata' must be a data frame")
    }
    
    # Convert to matrix for consistent processing
    indata <- as.matrix(indata)
    
    # Calculate variance for each row (feature)
    rv <- apply(indata, 1, var)
    
    # Select top variable features
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    
    # Perform PCA on transposed matrix of selected features
    pca <- tryCatch(
        prcomp(t(indata[select, ]), rank. = length(pcsToUse)),
        error = function(e) {
            stop("PCA calculation failed. Check input data for potential issues.")
        }
    )
    
    # Calculate percentage of variance explained
    percentVar <- pca$sdev^2 / sum(pca$sdev^2)
    
    # Validate intgroup
    if (!all(intgroup %in% names(coldata))) {
        stop("The specified 'intgroup' does not match any column in 'coldata'")
    }
    
    # Prepare group information
    group <- tryCatch(
        {
            if (is.factor(coldata[[intgroup]]) || is.character(coldata[[intgroup]])) {
                factor(coldata[[intgroup]])
            } else {
                stop("'intgroup' must be a factor or character vector")
            }
        },
        error = function(e) {
            stop("Error processing group information: ", e$message)
        }
    )
    
    # Create result data frame
    d <- setNames(
        as.data.frame(pca$x[, pcsToUse, drop = FALSE]), 
        paste0("PC", pcsToUse)
    )
    
    # Add group and sample information
    d$group <- group
    d$samples <- colnames(indata)
    
    # Attach percentage of variance as an attribute
    attr(d, "percentVar") <- percentVar[pcsToUse]
    
    return(d)
}
