#' Single Cell Visualization Function
#' @param seurat_obj Seurat对象
#' @param reduction 降维方法 (例如: "umap", "tsne", "harmony_umap")
#' @param group_by 分组变量 (例如: "seurat_clusters", "celltype", "orig.ident")
#' @param split_by 是否需要分面显示的变量(可选)
#' @param colors 自定义颜色向量(可选)
#' @param point_size 点的大小
#' @param save 是否保存图片
#' @param save_path 保存路径
#' @param width 图片宽度
#' @param height 图片高度
#' @return ggplot对象

# 主题设置函数
get_theme_settings <- function() {
  tidydr::theme_dr() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    )
}

plotSCdim <- function(seurat_obj,
                       reduction = "umap",
                       group_by,
                       split_by = NULL,
                       colors = NULL,
                       point_size = 0.1,
                       save = FALSE,
                       save_path = NULL,
                       width = 8,
                       height = 6) {
  
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(ggrastr)
  
  # 获取降维坐标
  dim_coords <- Embeddings(seurat_obj, reduction = reduction) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell")
  
  # 获取分组信息
  metadata <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    select(cell, !!sym(group_by))
  
  # 如果需要split，添加split信息
  if (!is.null(split_by)) {
    metadata <- metadata %>%
      mutate(split_group = seurat_obj@meta.data[[split_by]])
  }
  
  # 合并数据
  plot_data <- dim_coords %>%
    left_join(metadata, by = "cell")
  
  # 设置坐标轴名称
  dim_names <- colnames(dim_coords)[-1]
  
  # 基础图层
  p <- ggplot(plot_data, 
              aes_string(x = dim_names[1], 
                        y = dim_names[2], 
                        color = group_by)) +
    ggrastr::geom_point_rast(size = point_size) +
    get_theme_settings()  # 使用自定义主题
  
  # 添加颜色
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
  
  # 添加分面
  if (!is.null(split_by)) {
    p <- p + facet_wrap(~split_group)
  }
  
  # 设置图例
  p <- p + guides(
    color = guide_legend(
      keywidth = 0.4,
      ncol = ifelse(length(unique(plot_data[[group_by]])) > 12, 2, 1),
      override.aes = list(size = 2, stroke = 2)
    )
  )
  
  # 保存图片
  if (save && !is.null(save_path)) {
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
