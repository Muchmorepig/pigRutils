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
#' sample_result <- calcConDist(dat.tb, 
#'                             bySample = TRUE,
#'                             colname.sample = "sample_id")
#' }
#' 
#' @importFrom data.table data.table as.data.table melt dcast :=
#' @importFrom stats fisher.test chisq.test t.test wilcox.test p.adjust
#' 
#' @export
#' @name calcConDist

library(data.table)

calcConDist <- function(dat.tb, bySample = FALSE,
                        colname.cluster = "seurat_clusters",
                        colname.sample = NULL,
                        colname.condition = "orig.ident",
                        method = "chisq", min.rowSum = 0) {
  # 
  required_cols <- c(colname.cluster, colname.condition)
  if (bySample) {
    stopifnot(!is.null(colname.sample))
    required_cols <- c(required_cols, colname.sample)
  }
  
  missing_cols <- setdiff(required_cols, colnames(dat.tb))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # 
  get_counts <- function(data) {
    table(data[[colname.cluster]], data[[colname.condition]])
  }
  
  # 
  .table.fisher <- function(count.dist) {
    count.dist.melt.tb <- melt(as.data.table(count.dist),
                               variable.name = "cid", value.name = "count",
                               id.vars = c("V1", "V2")
    )
    
    sum.row <- rowSums(count.dist)
    sum.col <- colSums(count.dist)
    
    count.dist.melt.ext.tb <- count.dist.melt.tb[,
                                                 {
                                                   rid <- V1
                                                   cid <- V2
                                                   m <- matrix(c(
                                                     count, sum.row[rid] - count, sum.col[cid] - count,
                                                     sum(sum.col) - sum.row[rid] - (sum.col[cid] - count)
                                                   ), ncol = 2)
                                                   res.test <- fisher.test(m)
                                                   .(p.value = res.test$p.value, OR = res.test$estimate)
                                                 },
                                                 by = .(V1, V2)
    ] 
    
    count.dist.melt.ext.tb[, p.adj := p.adjust(p.value, "BH")]
    return(count.dist.melt.ext.tb)
  }
  
  # 
  if (method == "freq") {
    freq_table <- dat.tb[, .N, by = .(sampleID, loc, majorCluster)][
      dcast(.SD, sampleID + loc ~ majorCluster, value.var = "N", fill = 0),
      on = .(sampleID, loc)
    ]
    
    freq_table[, freq.mcls := N / sum(N), by = .(sampleID, loc)]
    res <- freq_table[, lapply(unique(loc), function(loc_val) {
      freq.x <- freq.mcls[loc == loc_val]
      freq.y <- freq.mcls[loc != loc_val]
      if (length(freq.x) >= 3) {
        logFC <- log2(mean(freq.x) / mean(freq.y))
        p.value.t <- t.test(freq.x, freq.y)$p.value
        p.value.w <- wilcox.test(freq.x, freq.y)$p.value
        
        data.table(loc = loc_val, logFC = logFC, p.value.t = p.value.t, p.value.w = p.value.w)
      } else {
        NULL
      }
    }), by = majorCluster][order(majorCluster)]
    
    res[, `:=`(
      FDR.t = p.adjust(p.value.t, "BH"), FDR.w = p.adjust(p.value.w, "BH"),
      char.sig = fcase(FDR.t < 0.01, "\U2731\U2731", FDR.t < 0.05, "\U2731", default = "")
    )]
    
    return(list(
      res = res,
      dist.logFC.tb = dcast(res, majorCluster ~ loc, value.var = "logFC"),
      dist.FDR.tb = dcast(res, majorCluster ~ loc, value.var = "FDR.t"),
      dist.charSig.tb = dcast(res, majorCluster ~ loc, value.var = "char.sig")
    ))
  }
  
  # 
  count.dist <- get_counts(dat.tb)
  
  if (method == "chisq") {
    res.chisq <- chisq.test(count.dist)
    return(res.chisq$observed / res.chisq$expected)
  } else if (method == "fisher") {
    count.dist_filtered <- count.dist[rowSums(count.dist) > min.rowSum, , drop = FALSE]
    count.dist.melt.ext.tb <- .table.fisher(count.dist_filtered)
    return(list(
      count.dist = count.dist.melt.ext.tb,
      p.tb = dcast(count.dist.melt.ext.tb, V1 ~ V2, value.var = "p.value"),
      padj.tb = dcast(count.dist.melt.ext.tb, V1 ~ V2, value.var = "p.adj"),
      OR.tb = dcast(count.dist.melt.ext.tb, V1 ~ V2, value.var = "OR")
    ))
  }
  
  # 
  if (bySample) {
    N.o.bySample <- table(dat.tb[[colname.sample]], dat.tb[[colname.cluster]], dat.tb[[colname.condition]])
    results <- apply(N.o.bySample, 1, function(x) {
      if (method == "chisq") {
        res.chisq <- chisq.test(x)
        return(res.chisq$observed / res.chisq$expected)
      } else {
        res.fisher <- .table.fisher(matrix(x,
                                           nrow = length(unique(dat.tb[[colname.cluster]])),
                                           ncol = length(unique(dat.tb[[colname.condition]]))
        ))
        return(dcast(res.fisher, V1 ~ V2, value.var = "OR"))
      }
    })
    return(results)
  }
  
  stop("Unsupported method")
}
