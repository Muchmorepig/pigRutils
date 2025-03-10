#' Calculate Effective Gene Lengths from GTF File (Alternative Method)
#'
#' @description
#' Calculates effective gene lengths using direct GTF parsing and IRanges
#' for efficient interval operations.
#'
#' @param gtf_file Character string specifying path to GTF file
#' @param cores Integer number of cores to use
#' @param chunk_size Integer size of chunks for parallel processing
#' @param verbose Logical for progress messages
#'
#' @return Data frame with gene IDs and effective lengths
#'
#' @import rtracklayer
#' @import IRanges
#' @import parallel
#' @import GenomicRanges
#'
#' @export
calcGeneEL_alt <- function(gtf_file, cores = NULL, chunk_size = 1000, verbose = TRUE) {
    # Input validation
    if (!file.exists(gtf_file)) {
        stop("GTF file does not exist: ", gtf_file)
    }

    # Set up parallel processing
    if (is.null(cores)) {
        cores <- max(1, floor(0.1 * parallel::detectCores()))
    }

    msg <- function(...) {
        if (verbose) message(paste0("[calcGeneEL_alt] ", ...))
    }

    # Create and register cluster
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))

    # Load required packages on worker nodes
    parallel::clusterEvalQ(cl, {
        library(IRanges)
        library(GenomicRanges)
    })

    tryCatch({
        msg("Importing GTF file...")
        gtf <- rtracklayer::import(gtf_file)

        # Convert to data frame and filter exons
        msg("Processing exon data...")
        exon_data <- as.data.frame(gtf[gtf$type == "exon"])

        # Create GRanges object for efficient interval operations
        exon_ranges <- GRanges(
            seqnames = exon_data$seqnames,
            ranges = IRanges(
                start = exon_data$start,
                end = exon_data$end
            ),
            strand = exon_data$strand,
            gene_id = exon_data$gene_id
        )

        # Split by gene ID
        msg("Splitting data by gene...")
        gene_groups <- split(exon_ranges, exon_ranges$gene_id)

        # Process in chunks for better memory management
        chunk_indices <- split(seq_along(gene_groups),
                             ceiling(seq_along(gene_groups) / chunk_size))

        msg(sprintf("Processing %d genes in %d chunks...",
                   length(gene_groups), length(chunk_indices)))

        # Process chunks in parallel
        results <- lapply(chunk_indices, function(chunk_idx) {
            current_genes <- gene_groups[chunk_idx]

            # Calculate effective lengths for current chunk
            chunk_lengths <- parallel::parLapply(cl, current_genes, function(gene_exons) {
                # Merge overlapping ranges and calculate total width
                reduced_ranges <- reduce(ranges(gene_exons))
                sum(width(reduced_ranges))
            })

            # Create data frame for current chunk
            data.frame(
                geneid = names(current_genes),
                efflen = as.numeric(chunk_lengths),
                stringsAsFactors = FALSE
            )
        })

        # Combine results
        msg("Combining results...")
        final_results <- do.call(rbind, results)

        # Add summary statistics
        if (verbose) {
            msg(sprintf("Processed %d genes", nrow(final_results)))
            msg(sprintf("Mean effective length: %.2f", mean(final_results$efflen)))
            msg(sprintf("Median effective length: %.2f", median(final_results$efflen)))
        }

        return(final_results)

    }, error = function(e) {
        stop("Error in calculating gene lengths: ", e$message)
    })
}

#' Compare Gene Length Calculation Methods
#'
#' @param gtf_file Path to GTF file
#' @return Data frame comparing results from both methods
#'
#' @export
compareGeneLengthMethods <- function(gtf_file) {
    # Calculate lengths using both methods
    msg <- message
    msg("Calculating lengths using method 1...")
    lengths1 <- calcGeneEL(gtf_file)

    msg("Calculating lengths using method 2...")
    lengths2 <- calcGeneEL_alt(gtf_file)

    # Merge results
    comparison <- dplyr::inner_join(
        lengths1, lengths2,
        by = "geneid",
        suffixes = c("_method1", "_method2"),
        all = TRUE
    )

    # Calculate differences
    comparison$diff <- abs(comparison$efflen_method1 - comparison$efflen_method2)

    # Summary statistics
    msg("\nComparison Summary:")
    msg("Total genes: ", nrow(comparison))
    msg("Genes with differences: ", sum(comparison$diff > 0))
    msg("Max difference: ", max(comparison$diff))
    msg("Mean difference: ", mean(comparison$diff))

    return(comparison)
}
