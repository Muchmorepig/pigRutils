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
                                  colname.sample = "orig.ident",
                                  colname.condition = "group",
                                  method = c("mixed_model", "sample_level_chisq", "fisher_exact"),
                                  min.cells.per.sample = 50,
                                  effect_threshold = 1.2) {
  
  # 检查必要包安装
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("需要安装data.table包：install.packages('data.table')")
  }
  if (method == "mixed_model" && !requireNamespace("lme4", quietly = TRUE)) {
    stop("需要安装lme4包：install.packages('lme4')")
  }
  
  # 数据格式转换和验证
  dt <- data.table::as.data.table(dat.tb)
  required_cols <- c(colname.cluster, colname.sample, colname.condition)
  missing_cols <- setdiff(required_cols, colnames(dt))
  if (length(missing_cols) > 0) {
    stop("缺少必要列：", paste(missing_cols, collapse = ", "))
  }
  
  # 检查条件列与样本列是否相同
  if (colname.condition == colname.sample) {
    stop("条件列不能与样本列相同，请指定独立的实验条件列")
  }
  
  # 数据预处理
  message("正在预处理数据...")
  dt[, sample_cell_count := .N, by = c(colname.sample)]
  dt_filtered <- dt[sample_cell_count >= min.cells.per.sample]
  
  # 检查数据是否有效
  if (nrow(dt_filtered) == 0) {
    stop("过滤后无有效数据，请降低min.cells.per.sample阈值")
  }
  
  # 选择统计方法
  method <- match.arg(method)
  message("使用分析方法：", method)
  
  # 核心分析逻辑
  if (method == "mixed_model") {
    # 混合效应模型方法
    message("正在拟合混合效应模型...")
    
    # 生成汇总数据
    count_data <- dt_filtered[, .(N = .N), 
      by = c(colname.sample, colname.cluster, colname.condition)]
    
    # 重命名列以适应模型公式
    data.table::setnames(count_data, 
      old = c(colname.sample, colname.cluster, colname.condition),
      new = c("sample", "cluster", "condition"))
    
    # 拟合模型
    model <- lme4::glmer(
      formula = N ~ condition + (1|sample) + (1|cluster),
      data = count_data,
      family = poisson()
    )
    
    # 提取结果
    res <- data.table::as.data.table(broom.mixed::tidy(model)) %>% 
      .[term == "condition", .(estimate, std.error, p.value)] %>% 
      .[, `:=`(
        effect_size = exp(estimate),
        FDR = p.adjust(p.value, "BH"),
        significance = data.table::fcase(
          FDR < 0.01 & abs(effect_size) > effect_threshold, "****",
          FDR < 0.05 & abs(effect_size) > effect_threshold, "***",
          FDR < 0.1 & abs(effect_size) > effect_threshold, "*",
          default = "ns"
        )
      )]
    
    return(list(
      model = model,
      results = res,
      filtered_data = dt_filtered,
      model_summary = summary(model)
    ))
    
  } else if (method == "sample_level_chisq") {
    # 样本层级的卡方检验
    message("执行样本层级卡方检验...")
    
    # 生成三维列联表
    crosstab <- data.table::dcast(
      dt_filtered,
      formula = reformulate(
        termlabels = c(colname.sample, colname.cluster),
        response = colname.condition),
      fun.aggregate = length
    )
    
    # 执行卡方检验
    chisq_test <- stats::chisq.test(crosstab[, -1])
    
    return(list(
      observed = chisq_test$observed,
      expected = chisq_test$expected,
      residuals = chisq_test$residuals,
      p.value = chisq_test$p.value,
      standardized_residuals = chisq_test$residuals / sqrt(sum(chisq_test$residuals^2))
    ))
    
  } else if (method == "fisher_exact") {
    # 簇特异的Fisher精确检验
    message("执行簇特异性Fisher检验...")
    
    results <- dt_filtered[, {
      # 对每个cluster进行分析
      cont_table <- table(get(colname.condition))
      if (all(dim(cont_table) == c(2,2))) {
        ft <- fisher.test(cont_table)
        .(p.value = ft$p.value, estimate = ft$estimate)
      } else {
        .(p.value = NA_real_, estimate = NA_real_)
      }
    }, by = c(colname.cluster)]
    
    # 多重检验校正
    results[, FDR := p.adjust(p.value, "BH")]
    
    return(results)
  }
}
