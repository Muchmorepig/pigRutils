#' @title GO Enrichment Analysis Tools
#' @name go_analysis
#' @description A collection of functions for GO enrichment analysis and visualization
#' @import ggplot2
#' @import clusterProfiler
#' @importFrom stats reorder
#' @importFrom methods is
#' @importFrom writexl write_xlsx
#' @importFrom grDevices pdf
#' @docType package
NULL

#' Perform Gene Ontology (GO) Enrichment Analysis
#'
#' @description This function conducts GO enrichment analysis on a set of gene identifiers
#' using clusterProfiler and simplifies the results.
#'
#' @param gene_ids A character vector of gene identifiers
#' @param orgdb An OrgDb object for the specific organism (e.g., org.Hs.eg.db for human)
#' @param ont A character string specifying the ontology type ("BP", "CC", or "MF")
#' @param keyType The key type of input gene IDs
#' @param pvalueCutoff Adjusted p-value cutoff (default: 0.05)
#' @param qvalueCutoff Q-value cutoff (default: 0.2)
#' @param simplify_cutoff Similarity cutoff for simplifying GO terms (default: 0.7)
#' @param use_cache Logical, whether to use cached results (default: TRUE)
#' @param cache_dir Directory for storing cached results (default: "cache")
#'
#' @return A simplified enrichGO result object or NULL if no enrichment is found
#' @examples
#' \dontrun{
#' library(org.At.tair.db)
#' orgdb <- org.At.tair.db
#' genes <- c("AT1G75240", "AT2G17950", ...)
#' result <- pGO(genes, orgdb, keyType = "TAIR")
#' }
#' @export
pGO <- function(
    gene_ids,
    orgdb,
    ont = "BP",
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    simplify_cutoff = 0.7
) {
    # Input validation
    if (length(gene_ids) == 0) {
        warning("No gene IDs provided.")
        return(NULL)
    }
    
    if (!is(orgdb, "OrgDb")) {
        stop("orgdb must be an OrgDb object")
    }
    
    if (!keyType %in% keytypes(orgdb)) {
        stop(sprintf("Invalid keyType. Must be one of: %s",
                    paste(keytypes(orgdb), collapse = ", ")))
    }
    
    valid_onts <- c("BP", "CC", "MF")
    if (!ont %in% valid_onts) {
        stop(sprintf("Invalid ontology type. Must be one of: %s",
                    paste(valid_onts, collapse = ", ")))
    }
    
    # Remove duplicates
    gene_ids <- unique(gene_ids)
    
    # Perform GO enrichment analysis with error handling
    result <- tryCatch({
        # Perform GO enrichment
        message("Running ...")
        enrich_res <- clusterProfiler::enrichGO(
            gene = gene_ids,
            OrgDb = orgdb,
            keyType = keyType,
            ont = ont,
            pvalueCutoff = pvalueCutoff,
            qvalueCutoff = qvalueCutoff
        )
        
        # Check if enrichment results are empty
        if (is.null(enrich_res) || nrow(enrich_res@result) == 0) {
            message("No significant GO terms found.")
            return(NULL)
        }
        
        # Simplify GO terms
        message("Simplifying GO terms...")
        simplified_res <- clusterProfiler::simplify(
            enrich_res,
            cutoff = simplify_cutoff,
            by = "p.adjust",
            select_fun = min
        )
        
        # Check if simplification resulted in any terms
        if (is.null(simplified_res) || nrow(simplified_res@result) == 0) {
            message("No GO terms remained after simplification.")
            return(enrich_res)
        }
        
        simplified_res
        
    }, error = function(e) {
        warning(sprintf("Error in GO analysis: %s", e$message))
        return(NULL)
    })
    
    # Cache results if successful
    if (!is.null(result) && use_cache) {
        saveRDS(result, cache_file)
    }
    
    return(result)
}

#' Create a Bar Plot for GO Enrichment Results
#'
#' @description Generates a customized bar plot for GO enrichment analysis results
#'
#' @param data A data frame containing GO enrichment results
#' @param mapping Aesthetic mapping
#' @param barfill Color of bars
#' @param baralpha Bar transparency
#' @param title_name Plot title
#' @param fontsize Base font size
#' @param text_angle Text rotation angle
#' @param text_size Text size
#' @param show_count Show gene counts (logical)
#'
#' @return A ggplot object
#' @export
gobar <- function(
    data,
    mapping = aes(
        x = reorder(Description, -log10(p.adjust)),
        y = -log10(p.adjust),
        label = Description
    ),
    barfill = "#ff4d79",
    baralpha = 0.9,
    title_name = "",
    fontsize = 12,
    text_angle = 0,
    text_size = 5,
    show_count = TRUE
) {
    # Check required packages
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package is required but not installed.")
    }
    
    # Create base plot
    p <- ggplot(data, mapping) +
        geom_bar(
            stat = "identity",
            position = position_stack(reverse = TRUE),
            fill = barfill,
            alpha = baralpha
        ) +
        theme_minimal() +
        coord_flip() +
        labs(
            x = "",
            y = "-log10(padj)",
            title = title_name
        ) +
        theme(
            text = element_text(size = fontsize),
            legend.position = "none",
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.text.y = element_text(size = fontsize)
        )
    
    # Add gene counts if requested
    if (show_count) {
        p <- p + geom_text(
            aes(label = sprintf("%s\n(n=%d)", Description, Count)),
            hjust = 0,
            position = position_stack(vjust = 0.5),
            size = text_size,
            angle = text_angle
        )
    } else {
        p <- p + geom_text(
            aes(label = Description),
            hjust = 0,
            position = position_stack(vjust = 0.5),
	    size = text_size,
            angle = text_angle
        )
    }
    
    return(p)
}

#' Plot GO Enrichment Results with Multiple Visualization Options
#'
#' @param go_result GO enrichment result object
#' @param plot_type Type of plot ("bar", "dot", "network", "heatmap")
#' @param top_n Number of top terms to plot
#' @param output_file Output file path (optional)
#' @param ... Additional parameters passed to plotting functions
#'
#' @return A plot object or NULL if no results
#' @export
plotGO <- function(
    go_result,
    plot_type = "bar",
    top_n = 15,
    output_file = NULL,
    ...
) {
    if (is.null(go_result) || nrow(go_result@result) == 0) {
        warning("No GO enrichment results to plot.")
        return(NULL)
    }
    
    # Extract top N results
    plot_data <- go_result@result[
        order(go_result@result$p.adjust)[1:min(top_n, nrow(go_result@result))],
    ]
    
    # Create plot based on type
    p <- switch(
        plot_type,
        "bar" = gobar(plot_data, ...),
        # "dot" = dotplot(go_result, showCategory = top_n, ...),
        # "network" = emapplot(go_result, showCategory = top_n, ...),
        # "heatmap" = heatplot(go_result, showCategory = top_n, ...),
        stop("Invalid plot type")
    )
    
    # Save plot if output file specified
    if (!is.null(output_file)) {
        ggsave(output_file, p, ...)
    }
    
    return(p)
}

#' Export GO Enrichment Results
#'
#' @param go_result GO enrichment result object or a list of GO result objects
#' @param file Output file path (must end in .xlsx or .csv)
#' @param sheets Named list of data frames for Excel export (optional)
#'
#' @return Invisible NULL
#' @export
export_GO <- function(go_result, file, sheets = NULL) {
    if (is.null(go_result)) {
        warning("No GO enrichment results to export.")
        return(invisible(NULL))
    }
    
    # Function to extract results from a single GO result
    extract_results <- function(go_obj) {
        if (is.null(go_obj)) return(NULL)
        go_obj@result
    }
    
    # Handle different input types
    if (grepl("\\.xlsx$", file)) {
        if (is.list(go_result) && !is(go_result, "enrichResult")) {
            # Handle list of GO results
            sheets <- lapply(go_result, extract_results)
            
            # If the list is named, use those names for sheets
            if (!is.null(names(go_result))) {
                names(sheets) <- names(go_result)
            } else {
                # Create default sheet names if list is unnamed
                names(sheets) <- paste0("GO_Results_", seq_along(sheets))
            }
        } else {
            # Handle single GO result
            if (is.null(sheets)) {
                sheets <- list(GO_Results = extract_results(go_result))
            }
        }
        writexl::write_xlsx(sheets, path = file)
        
    } else if (grepl("\\.csv$", file)) {
        if (is.list(go_result) && !is(go_result, "enrichResult")) {
            warning("Multiple GO results cannot be exported to a single CSV file. Only the first result will be exported.")
            results <- extract_results(go_result[[1]])
        } else {
            results <- extract_results(go_result)
        }
        write.csv(results, file = file, row.names = FALSE)
    } else {
        stop("Unsupported file format. Use .xlsx or .csv")
    }
    
    message(sprintf("Results exported to %s", file))
    return(invisible(NULL))
}

#' Run Batch GO Analysis
#'
#' @param gene_lists Named list of gene vectors
#' @param orgdb OrgDb object
#' @param ... Additional parameters passed to pGO
#'
#' @return List of GO enrichment results
#' @export
batch_GO <- function(gene_lists, orgdb, ...) {
    if (!is.list(gene_lists) || is.null(names(gene_lists))) {
        stop("gene_lists must be a named list of gene vectors")
    }
    
    results <- lapply(names(gene_lists), function(name) {
        message(sprintf("Processing: %s", name))
        result <- pGO(gene_lists[[name]], orgdb, ...)
        return(result)
    })
    
    names(results) <- names(gene_lists)
    return(results)
}
