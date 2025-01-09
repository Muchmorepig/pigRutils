#' Generate UMAP plot for ArchR object
#'
#' This function creates a UMAP plot for an ArchR project, allowing visualization
#' of cells colored by clusters or gene expression.
#'
#' @param ArchR.proj An ArchR project object.
#' @param embed.umap Character string specifying the name of the UMAP embedding in the ArchR project. Default is "UMAP_Tile_8".
#' @param colorBy Character string specifying either a cluster name in cellColData or a gene name. Default is "Cluster_Tile_8".
#' @param mtx Optional SummarizedExperiment object of PeakMatrix or GeneScoreMatrix. Default is NULL.
#' @param marker.color Vector of three colors for low, mid, and high values when coloring by gene expression. Default is c(low = "#eafff2", mid = "#fad8bf", high = "#ff0000").
#' @param palette Color palette for clusters. Default is paletteer::paletteer_d("ggsci::category20c_d3").
#' @param multiple Logical, whether to return a single plot with all clusters or a faceted plot with single clusters. Default is FALSE.
#' @param label.size Numeric, size of cluster labels. Default is 3.6.
#' @param label.shadow Numeric, size of label shadow. Default is 0.25.
#' @param label.color Character, color of cluster labels. Default is "#000103".
#' @param shadow.color Character, color of label shadow. Default is "#f5f0f0".
#' @param point.size Numeric, size of points in the plot. Default is 0.8.
#' @param point.alpha Numeric, alpha transparency of points. Default is 0.8.
#' @param theme Character, either "mini" or "clear", specifying the plot theme. Default is "mini".
#' @param title Character, optional title for the plot. Default is NULL.
#' @param title.size Numeric, size of the plot title. Default is 13.
#' @param legend.size Numeric, size of the legend text. Default is 10.
#' @param log2Norm Logical, whether to apply log2 normalization to gene expression values. Default is TRUE.
#'
#' @return A ggplot object representing the UMAP plot.
#'
#' @import dplyr
#' @import ggplot2
#' @import gtools
#' @import SummarizedExperiment
#' @import ggrepel
#' @importFrom tidydr theme_dr
#'
#' @examples
#' \dontrun{
#' # Plot UMAP colored by clusters
#' umapForArchR(ArchR.proj = myArchRProject, colorBy = "Clusters")
#'
#' # Plot UMAP colored by gene expression
#' umapForArchR(ArchR.proj = myArchRProject, colorBy = "GENE1", mtx = geneScoreMatrix)
#' }
#'
#' @export
umapForArchR <- function(ArchR.proj,
                         embed.umap = "UMAP_Tile_8",
                         colorBy = "Cluster_Tile_8",
                         mtx = NULL,
                         marker.color = c(low = "#eafff2", mid = "#fad8bf", high = "#ff0000"),
                         palette = paletteer::paletteer_d("ggsci::category20c_d3"),
                         multiple = FALSE,
                         label.size = 3.6,
                         label.shadow = 0.25,
                         label.color = "#000103",
                         shadow.color = "#f5f0f0",
                         point.size = 0.8,
                         point.alpha = 0.8,
                         theme = c("mini", "clear"),
                         title = NULL,
                         title.size = 13,
                         legend.size = 10,
                         log2Norm = TRUE) {
  
  # Input validation
  if (!embed.umap %in% names(ArchR.proj@embeddings)) {
    stop("`embed.umap` must be included in `names(ArchR.proj@embeddings)`")
  }

  # Extract UMAP coordinates
  df <- ArchR.proj@embeddings[[embed.umap]]$df
  colnames(df) <- c("UMAP_Dim1", "UMAP_Dim2")

  if (is.null(mtx)) {
    # Color by cluster
    celcol <- colnames(ArchR.proj@cellColData)
    colorBy <- celcol[match(colorBy, celcol)]
    if (is.na(colorBy)) {
      stop("`colorBy` must be included in `colnames(ArchR.proj@cellColData)`")
    }
    df$Clusters <- factor(ArchR.proj@cellColData[, colorBy], levels = gtools::mixedsort(unique(ArchR.proj@cellColData[, colorBy])))
    
    if (multiple) {
      p <- create_multiple_plot(df, palette, point.size, legend.size)
    } else {
      p <- create_single_plot(df, palette, point.size, point.alpha, label.size, label.shadow, label.color, shadow.color, legend.size)
    }
  } else {
    # Color by gene expression
    p <- create_gene_expression_plot(df, mtx, colorBy, marker.color, point.size, point.alpha, log2Norm)
    title <- colorBy
  }

  # Apply theme
  p <- apply_theme(p, match.arg(theme), title, title.size, legend.size)

  return(p)
}

# Helper functions

create_multiple_plot <- function(df, palette, point.size, legend.size) {
  ggplot(df, aes(UMAP_Dim1, UMAP_Dim2)) +
    geom_point(data = df[, 1:2], color = "#e7e7e7", size = point.size) +
    geom_point(size = point.size, aes(color = Clusters)) +
    facet_wrap(~Clusters) +
    scale_color_manual(values = palette, breaks = gtools::mixedsort(unique(df$Clusters))) +
    labs(x = "UMAP Dim1", y = "UMAP Dim2") +
    .mytheme +
    theme(legend.text = element_text(size = legend.size)) +
    guides(color = guide_legend(keywidth = 0.4, nrow = 1, override.aes = list(size = 2, stroke = 2)))
}

create_single_plot <- function(df, palette, point.size, point.alpha, label.size, label.shadow, label.color, shadow.color, legend.size) {
  df2 <- df %>%
    dplyr::group_by(Clusters) %>%
    dplyr::summarize(UMAP_Dim1 = mean(UMAP_Dim1), UMAP_Dim2 = mean(UMAP_Dim2))

  ggplot(df, aes(UMAP_Dim1, UMAP_Dim2)) +
    geom_point(size = point.size, alpha = point.alpha, aes(color = Clusters)) +
    ggrepel::geom_text_repel(
      data = df2,
      aes(label = Clusters),
      color = label.color,
      bg.color = shadow.color,
      bg.r = label.shadow,
      size = label.size,
      box.padding = 0,
      seed = 1
    ) +
    scale_color_manual(values = palette, breaks = gtools::mixedsort(df2$Clusters)) +
    guides(color = guide_legend(keywidth = 0.4, override.aes = list(size = 2, stroke = 2)))
}

create_gene_expression_plot <- function(df, mtx, colorBy, marker.color, point.size, point.alpha, log2Norm) {
  mx <- SummarizedExperiment::assay(mtx)
  rownames(mx) <- SummarizedExperiment::rowData(mtx)$name
  df$count <- mx[colorBy, rownames(df)]
  df <- df[order(df$count), ]
  if (log2Norm) df$count <- log2(df$count + 1)
  
  ggplot(df, aes(UMAP_Dim1, UMAP_Dim2)) +
    geom_point(size = point.size, alpha = point.alpha, aes(color = count)) +
    scale_color_gradient2(
      midpoint = max(df$count) * 0.5,
      low = marker.color["low"],
      mid = marker.color["mid"],
      high = marker.color["high"]
    )
}

apply_theme <- function(p, theme, title, title.size, legend.size) {
  if (theme == "mini") {
    p + labs(title = title, caption = "X-axis: UMAP Dim1\nY-axis: UMAP Dim2") +
      theme_minimal() +
      theme(
        legend.text = element_text(size = legend.size),
        legend.title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.background = element_rect(color = "#f3eaea"),
        plot.title = element_text(size = title.size)
      )
  } else if (theme == "clear") {
    p + labs(title = title, x = "UMAP Dim1", y = "UMAP Dim2") +
      tidydr::theme_dr() +
      theme(
        legend.text = element_text(size = legend.size),
        legend.title = element_blank(),
        legend.background = element_rect(color = "#f3eaea"),
        plot.title = element_text(size = title.size),
        panel.grid = element_blank()
      )
  } else {
    stop("theme is required: 'mini' or 'clear'")
  }
}

.mytheme <- theme_minimal() + theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  strip.text = element_blank(),
  panel.spacing = unit(0, "lines"),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  legend.title = element_blank(),
  legend.position = "top",
  legend.background = element_rect(color = "#f3eaea")
)
