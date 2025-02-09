#' Calculate Cell Type Distribution by Condition
#' 
#' @description
#' Calculates the distribution of cell types across different experimental conditions
#' using various statistical methods including chi-square test and Fisher's exact test.
#' 
#' @param dat.tb A data.table object containing the cell metadata
#' @param bySample Logical, whether to perform analysis by sample. Default is FALSE
#' @param colname.cluster Character, column name for cell clusters. Default is "seurat_clusters"
#' @param colname.sample Character, column name for sample IDs. Required if bySample = TRUE
#' @param colname.condition Character, column name for experimental conditions. Default is "orig.ident"
#' @param method Character, statistical method to use. Options are "chisq", "fisher", or "freq". Default is "chisq"
#' @param min.rowSum Numeric, minimum row sum threshold for filtering in fisher test. Default is 0
#' 
#' @return 
#' Depending on the method:
#' \itemize{
#'   \item For "chisq": Returns a matrix of observed/expected ratios
#'   \item For "fisher": Returns a list containing count distribution, p-values, adjusted p-values, and odds ratios
#'   \item For "freq": Returns a list containing frequency analysis results with logFC and significance tests
#' }
#' 
#' @examples
#' \dontrun{
#' # Basic usage with chi-square test
#' result <- calcConDist(dat.tb)
#' 
#' # Using Fisher's exact test
#' fisher_result <- calcConDist(dat.tb, method = "fisher")
#' 
#' # Analysis by sample
#' results <- calcConDist_v2(
#'   dat.tb = your_data,
#'   method = "mixed_model",
#'   min.cells.per.sample = 100
#' )}
#'
#' @importFrom data.table data.table as.data.table melt dcast :=
#' @importFrom stats fisher.test chisq.test t.test wilcox.test p.adjust
#' 
#' @export
#' @name calcConDist_v2

library(data.table)

calcConDist_v2 <- function(dat.tb,
                          colname.cluster = "seurat_clusters",
                          colname.sample = "sample_id",
                          colname.condition = "orig.ident",
                          method = "mixed_model",
                          min.cells.per.sample = 50) {
  # 参数检查
  required_cols <- c(colname.cluster, colname.sample, colname.condition)
  missing_cols <- setdiff(required_cols, colnames(dat.tb))
  if (length(missing_cols) > 0) {
    stop(sprintf("缺少必要列: %s", paste(missing_cols, collapse = ", ")))
  }

  # 过滤低质量样本
  dat.tb <- dat.tb[, sample_cells := .N, by = colname.sample][sample_cells >= min.cells.per.sample]

  # 核心统计方法
  if (method == "mixed_model") {
    # 使用混合效应模型
    library(lme4)
    library(broom.mixed)

    # 生成样本-簇-条件计数矩阵
    sample_counts <- dat.tb[, .(.N), by = c(colname.sample, colname.cluster, colname.condition)]

    # 拟合模型
    model <- glmer(
      formula = N ~ condition + (1|sample) + (1|cluster),
      data = sample_counts,
      family = poisson()
    )

    # 提取结果
    res <- tidy(model, effects = "fixed") %>%
      filter(term == colname.condition) %>%
      mutate(
        FDR = p.adjust(p.value, "BH"),
        effect_size = exp(estimate)
      )

    return(list(
      model = model,
      results = res,
      message = "使用混合效应模型（包含样本和簇的随机效应）"
    ))

  } else if (method == "sample_level_chisq") {
    # 样本层级的卡方检验

    # 生成三维列联表（样本 x 簇 x 条件）
    crosstab <- dcast(dat.tb,
                     formula = paste(colname.sample, "~", colname.cluster, "+", colname.condition),
                     value.var = colname.sample,
                     fun.aggregate = length)

    # 执行卡方检验
    chisq_res <- chisq.test(crosstab[, -1])  # 排除样本名列

    # 计算标准化残差
    standardized_resid <- chisq_res$residuals / sqrt(sum(chisq_res$residuals^2))

    return(list(
      observed = chisq_res$observed,
      expected = chisq_res$expected,
      residuals = standardized_resid,
      p.value = chisq_res$p.value
    ))

  } else if (method == "fisher_exact") {
    # 样本层级的Fisher精确检验
    library(plyr)

    res <- dlply(dat.tb, colname.cluster, function(cluster_data) {
      # 对每个簇进行样本级别的检验
      sample_counts <- cluster_data[, .N, by = c(colname.sample, colname.condition)]

      # 构建2x2表（condition x counts）
      contingency_table <- table(sample_counts[[colname.condition]], sample_counts$N)

      if (all(dim(contingency_table) == c(2,2))) {
        fisher.test(contingency_table)$p.value
      } else {
        NA_real_
      }
    })

    # 多重检验校正
    res_df <- data.frame(
      cluster = names(res),
      p.value = unlist(res),
      FDR = p.adjust(unlist(res), "BH")
    )

    return(res_df)
  }

  stop("不支持的统计方法，可用选项: mixed_model, sample_level_chisq, fisher_exact")
}
