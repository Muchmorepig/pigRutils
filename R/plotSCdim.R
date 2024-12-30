#' Plot Single Cell Dimensional Reduction
#'
#' This function creates dimensional reduction plots (UMAP, tSNE, etc.) for single cell data
#' using a consistent theme and styling.
#'
#' @param seurat_obj A Seurat object containing the single cell data
#' @param reduction Character string specifying which dimensional reduction to use.
#'                 Default: "umap"
#' @param group_by Character string specifying the metadata column to use for coloring.
#'                Must be present in seurat_obj@meta.data
#' @param split_by Optional character string specifying a metadata column for faceting.
#'                Must be present in seurat_obj@meta.data
#' @param colors Optional vector of colors for groups. If NULL, uses default color scheme
#' @param point_size Numeric value specifying point size. Default: 0.1
#' @param save Logical indicating whether to save the plot. Default: FALSE
#' @param save_path Character string specifying save location if save=TRUE
#' @param width Numeric value specifying plot width in inches. Default: 8
#' @param height Numeric value specifying plot height in inches. Default: 6
#'
#' @return A ggplot object
#' 
#' @import Seurat
#' @import ggplot2
#' @import dplyr
#' @import ggrastr
#' @importFrom tidydr theme_dr
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' plot_sc_dim(seurat_obj, group_by = "seurat_clusters")
#'
#' # With splitting
#' plot_sc_dim(seurat_obj, 
#'             group_by = "celltype", 
#'             split_by = "orig.ident",
#'             point_size = 0.2)
#'
#' # Custom colors
#' plot_sc_dim(seurat_obj,
#'             group_by = "orig.ident",
#'             colors = c("#1dc1ab", "#56acfc", "#fc877f"))
#' }
#'
#' @export
plot_sc_dim <- function(seurat_obj,
                       reduction = "umap",
                       group_by,
                       split_by = NULL,
                       colors = NULL,
                       point_size = 0.1,
                       save = FALSE,
                       save_path = NULL,
                       width = 8,
                       height = 6) {
  
  # Input validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }
  
  if (!reduction %in% Reductions(seurat_obj)) {
    stop(sprintf("Reduction '%s' not found in Seurat object", reduction))
  }
  
  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("group_by variable '%s' not found in metadata", group_by))
  }
  
  if (!is.null(split_by) && !split_by %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("split_by variable '%s' not found in metadata", split_by))
  }
  
  if (save && is.null(save_path)) {
    stop("save_path must be specified when save = TRUE")
  }
  
  # Get dimensional reduction coordinates
  dim_coords <- Embeddings(seurat_obj, reduction = reduction) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")
  
  # Get metadata
  metadata <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, !!sym(group_by))
  
  # Add split information if needed
  if (!is.null(split_by)) {
    metadata <- metadata %>%
      dplyr::mutate(split_group = seurat_obj@meta.data[[split_by]])
  }
  
  # Merge data
  plot_data <- dim_coords %>%
    dplyr::left_join(metadata, by = "cell")
  
  # Create base plot
  dim_names <- colnames(dim_coords)[-1]
  p <- ggplot(plot_data, 
              aes_string(x = dim_names[1], 
                        y = dim_names[2], 
                        color = group_by)) +
    ggrastr::geom_point_rast(size = point_size) +
    get_theme_settings()
  
  # Add colors
  if (is.null(colors)) {
    n_groups <- length(unique(plot_data[[group_by]]))
    if (is.numeric(plot_data[[group_by]])) {
      p <- p + scale_color_gradientn(
        colours = c("#c8c6c3", "#f6ee6c", "#f9a432", "#eb212f", "#88181d")
      )
    } else {
      p <- p + scale_color_manual(values = pigRutils::select_colors("col_33")[1:n_groups])
    }
  } else {
    p <- p + scale_color_manual(values = colors)
  }
  
  # Add faceting if specified
  if (!is.null(split_by)) {
    p <- p + facet_wrap(~split_group)
  }
  
  # Add legend guides
  p <- p + guides(
    color = guide_legend(
      keywidth = 0.4,
      ncol = ifelse(length(unique(plot_data[[group_by]])) > 12, 2, 1),
      override.aes = list(size = 2, stroke = 2)
    )
  )
  
  # Save plot if requested
  if (save) {
    ggsave(
      filename = save_path,
      plot = p,
      width = width,
      height = height,
      dpi = 300
    )
  }
  
  return(p)
}

#' Get theme settings for single cell plots
#'
#' @return A ggplot2 theme object
#' @keywords internal
#'
#' @importFrom tidydr theme_dr
get_theme_settings <- function() {
  tidydr::theme_dr() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    )
}
