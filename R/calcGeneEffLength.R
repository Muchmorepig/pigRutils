#' Calculate Effective Gene Lengths from GTF File
#' 
#' @description
#' This function calculates the effective length of genes by summing the lengths
#' of their unique exonic regions, accounting for overlapping exons.
#' 
#' @param gtf_file Character string specifying the path to GTF/GFF file. Can be gzipped.
#' @param cores Integer specifying the number of CPU cores to use. If NULL, uses 50% of available cores.
#' @param skip_gene_check Logical indicating whether to skip gene ID validation (default: FALSE)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#' 
#' @return A data.frame with columns:
#'   \itemize{
#'     \item geneid: Character vector of gene identifiers
#'     \item efflen: Numeric vector of effective lengths
#'     \item n_exons: Integer vector of exon counts per gene (optional)
#'   }
#' 
#' @import GenomicFeatures
#' @import parallel
#' @import IRanges
#' @importFrom stats setNames
#' 
#' @examples
#' \dontrun{
#' # Basic usage
#' gene_lengths <- calcGeneEL("path/to/annotation.gtf.gz")
#' 
#' # Use specific number of cores
#' gene_lengths <- calcGeneEL("path/to/annotation.gtf.gz", cores = 4)
#' 
#' # Run quietly
#' gene_lengths <- calcGeneEL("path/to/annotation.gtf.gz", verbose = FALSE)
#' }
#' 
#' @export
calcGeneEL <- function(gtf_file, 
                      cores = NULL,
                      skip_gene_check = FALSE,
                      verbose = TRUE) {
    
    # Input validation
    if (!is.character(gtf_file) || length(gtf_file) != 1) {
        stop("gtf_file must be a single character string")
    }
    
    if (!file.exists(gtf_file)) {
        stop("GTF file does not exist: ", gtf_file)
    }
    
    # Function to print verbose messages
    msg <- function(...) {
        if (verbose) message(paste0("[calcGeneEL] ", ...))
    }
    
    # Set up parallel processing
    if (is.null(cores)) {
        cores <- max(1, floor(0.1 * parallel::detectCores()))
    }
    cores <- as.integer(cores)
    if (cores < 1) stop("cores must be a positive integer")
    
    msg("Using ", cores, " cores")
    
    # Create cluster with error handling
    cl <- tryCatch({
        parallel::makeCluster(cores)
    }, error = function(e) {
        stop("Failed to create cluster: ", e$message)
    })
    on.exit(parallel::stopCluster(cl))
    
    # Load required packages on worker nodes
    msg("Setting up worker nodes...")
    parallel::clusterEvalQ(cl, {
        suppressPackageStartupMessages({
            library(IRanges)
            library(GenomicRanges)
        })
    })
    
    # Main processing
    tryCatch({
        msg("Creating TxDb from GTF...")
        txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf")
        
        # Validate gene IDs if requested
        if (!skip_gene_check) {
            msg("Validating gene IDs...")
            genes <- genes(txdb)
            if (length(genes) == 0) {
                stop("No genes found in GTF file")
            }
        }
        
        msg("Extracting exons by gene...")
        exons_gene <- GenomicFeatures::exonsBy(txdb, by = "gene")
        
        # Progress tracking for parallel processing
        msg("Calculating effective lengths...")
        total_genes <- length(exons_gene)
        counter <- 0
        
        # Calculate effective lengths with progress tracking
        exons_gene_lens <- parallel::parLapply(cl, exons_gene, function(x) {
            # Calculate basic stats
            n_exons <- length(x)
            total_len <- sum(IRanges::width(GenomicRanges::reduce(x)))
            
            # Return both length and exon count
            list(
                efflen = total_len,
                n_exons = n_exons
            )
        })
        
        # Process results
        msg("Processing results...")
        result <- data.frame(
            geneid = names(exons_gene_lens),
            efflen = vapply(exons_gene_lens, function(x) x$efflen, numeric(1)),
            n_exons = vapply(exons_gene_lens, function(x) x$n_exons, numeric(1)),
            stringsAsFactors = FALSE
        )
        
        # Add basic statistics if verbose
        if (verbose) {
            msg(sprintf("Processed %d genes", nrow(result)))
            msg(sprintf("Mean effective length: %.2f", mean(result$efflen)))
            msg(sprintf("Mean exons per gene: %.2f", mean(result$n_exons)))
        }
        
        return(result)
        
    }, error = function(e) {
        stop("Error in calculating gene lengths: ", e$message)
    })
}

#' Validate Gene Length Results
#' 
#' @param gene_lengths Data frame output from calcGeneEL
#' @return Logical indicating if the results pass basic validation
#' 
#' @export
validateGeneLengths <- function(gene_lengths) {
    # Basic structure checks
    if (!is.data.frame(gene_lengths)) {
        stop("Input must be a data frame")
    }
    
    required_cols <- c("geneid", "efflen")
    if (!all(required_cols %in% colnames(gene_lengths))) {
        stop("Missing required columns: ", 
             paste(setdiff(required_cols, colnames(gene_lengths)), collapse = ", "))
    }
    
    # Value checks
    if (any(gene_lengths$efflen <= 0, na.rm = TRUE)) {
        warning("Found genes with zero or negative effective length")
        return(FALSE)
    }
    
    if (any(duplicated(gene_lengths$geneid))) {
        warning("Found duplicate gene IDs")
        return(FALSE)
    }
    
    return(TRUE)
}
