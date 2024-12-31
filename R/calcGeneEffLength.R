#' Calculate Effective Gene Lengths from GTF File
#' 
#' @description 
#' This function calculates the effective length of genes by summing up the 
#' non-overlapping exon lengths for each gene from a GTF annotation file.
#' 
#' @param gtf_file Character string specifying the path to the GTF file
#' @param cores Numeric value indicating the number of CPU cores to use (default is half of available cores)
#' 
#' @return A data frame with two columns:
#'   \item{geneid}{Gene identifiers}
#'   \item{efflen}{Effective lengths of genes}
#'   
#' @import GenomicFeatures
#' @import parallel
#' 
#' @examples
#' \dontrun{
#' gene_lengths <- calcGeneEL("path/to/your/annotation.gtf.gz")
#' head(gene_lengths)
#' }
#' 
#' @export
calcGeneEL <- function(gtf_file, cores = NULL) {
  # Input validation
  if (!file.exists(gtf_file)) {
    stop("GTF file does not exist: ", gtf_file)
  }
  
  # Set up parallel processing
  if (is.null(cores)) {
    cores <- max(1, floor(0.5 * parallel::detectCores()))
  }
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl)) # Ensure cluster is stopped when function exits
  
  tryCatch({
    # Create TxDb object from GTF
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf")
    
    # Extract exons by gene
    exons_gene <- GenomicFeatures::exonsBy(txdb, by = "gene")
    
    # Calculate effective lengths using parallel processing
    exons_gene_lens <- parallel::parLapply(cl, exons_gene, function(x) {
      sum(width(reduce(x)))
    })
    
    # Create output data frame
    result <- data.frame(
      geneid = names(exons_gene_lens),
      efflen = as.numeric(exons_gene_lens),
      stringsAsFactors = FALSE
    )
    
    return(result)
    
  }, error = function(e) {
    stop("Error in calculating gene lengths: ", e$message)
  })
}
