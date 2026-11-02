#' Calculate Cell Type Distribution by Condition
#'
#' @description Describe cell composition or compare cell counts. Sample-level
#' tests require explicit opt-in and independent biological samples.
#' @param dat.tb A data.frame or data.table; one row per retained cell.
#' @param bySample Logical. For chisq/fisher, return a named list of results
#'   computed separately within each sample. Does not create biological replicates.
#'   For freq, proportions are always computed per sample.
#' @param colname.cluster Cluster column name. Unused factor levels are retained.
#' @param colname.sample Sample column name; required for bySample or freq.
#' @param colname.condition Condition column name.
#' @param method One of "chisq", "fisher", or "freq" (exact spelling).
#' @param min.rowSum Nonnegative threshold. Fisher reports clusters whose pooled
#'   count is strictly greater than this value; all cells remain in the background.
#' @param sample.inference Logical, default FALSE. For freq only, explicitly
#'   assert that sample IDs identify independent biological replicates and enable
#'   unpaired Welch and Wilcoxon tests. At least two samples on each side required.
#' @return chisq: observed/expected matrix, with counts, expected, p.value,
#'   status and test.warning attributes. fisher: list with count.dist, p.tb,
#'   padj.tb, OR.tb; long-table labels remain V1/V2. bySample returns a named
#'   list of these same structures. freq: list with res, dist.logFC.tb,
#'   dist.FDR.tb, dist.charSig.tb and sample.freq.tb. Legacy freq output labels
#'   majorCluster, loc, sampleID are fixed output names, regardless of input names.
#'   Undefined tests are NA with a status; no pseudocount is added.
#' @details
#' chisq and fisher use cells as observations, not biological replicates. They
#' cannot establish population-level condition effects. Chisq needs at least
#' two positive row and column margins; zero margins are retained as NA.
#' Fisher compares each cluster/condition against all other clusters/conditions;
#' zero 2x2 margins yield NA. BH adjustment covers all reported cluster-condition
#' tests (separately within each sample when bySample=TRUE).
#'
#' freq requires exactly one condition per sample (paired designs unsupported).
#' Missing clusters in observed samples are zero-filled; absent sample-condition
#' combinations are not invented. Means give equal weight to each sample. Each
#' condition is compared with all remaining samples pooled; with >2 conditions,
#' the reference is sample-weighted, not condition-weighted. logFC is log2 of the
#' mean proportion ratio: positive/zero gives Inf, zero/positive -Inf, 0/0 NA.
#' All-zero clusters are not tested. Each test's BH family includes all res rows;
#' for two conditions reciprocal contrasts are retained for compatibility.
#' Test errors/warnings are recorded in status.t/status.w; Wilcoxon uses a
#' tie-corrected normal approximation. Sample independence cannot be verified
#' from metadata. Technical libraries or random cell splits are not replicates.
#' @examples
#' d <- data.frame(cluster = c("A", "B", "A", "A"),
#'                 sample = c("s1", "s1", "s2", "s2"),
#'                 condition = c("M", "M", "F", "F"))
#' r <- calcConDist(d, colname.cluster = "cluster", colname.sample = "sample",
#'                  colname.condition = "condition", method = "freq")
#' r$sample.freq.tb
#' @export
calcConDist <- function(dat.tb, bySample = FALSE,
                        colname.cluster = "seurat_clusters",
                        colname.sample = NULL,
                        colname.condition = "orig.ident",
                        method = "chisq", min.rowSum = 0,
                        sample.inference = FALSE) {
  flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  if (!flag(bySample) || !flag(sample.inference))
    stop("bySample and sample.inference must be single non-missing logicals")
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% c("chisq", "fisher", "freq")) stop("Unsupported method")
  if (!is.numeric(min.rowSum) || length(min.rowSum) != 1L ||
      !is.finite(min.rowSum) || min.rowSum < 0) stop("Invalid min.rowSum")
  if (sample.inference && method != "freq")
    stop("sample.inference is only supported for method='freq'")
  if (!is.data.frame(dat.tb) || !nrow(dat.tb))
    stop("dat.tb must be a nonempty data.frame or data.table")
  if (anyDuplicated(names(dat.tb))) stop("Duplicate input column names")
  cols <- list(colname.cluster, colname.condition)
  if (bySample || method == "freq") cols <- c(cols, list(colname.sample))
  if (!all(vapply(cols, function(x) is.character(x) && length(x) == 1L &&
                  !is.na(x) && nzchar(x), logical(1))))
    stop("Column names must be nonempty strings; colname.sample is required for bySample or freq")
  required <- unique(unlist(cols))
  missing <- setdiff(required, names(dat.tb))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  for (nm in required) {
    x <- dat.tb[[nm]]
    if (!is.atomic(x) || !is.null(dim(x)) || anyNA(x) ||
        any(!nzchar(trimws(as.character(x)))))
      stop("Missing, blank or unsupported metadata in column: ", nm)
  }
  # Copy labels, never modify the caller's data.table by reference.
  lev <- function(x) if (is.factor(x)) levels(x) else sort(unique(as.character(x)))
  clusters <- lev(dat.tb[[colname.cluster]])
  conditions <- lev(dat.tb[[colname.condition]])
  d <- data.frame(cluster = as.character(dat.tb[[colname.cluster]]),
                  condition = as.character(dat.tb[[colname.condition]]),
                  stringsAsFactors = FALSE)
  if (bySample || method == "freq") d$sample <- as.character(dat.tb[[colname.sample]])
  counts <- function(x) table(factor(x$cluster, levels = clusters),
                              factor(x$condition, levels = conditions),
                              dnn = c("V1", "V2"))
  # Capture exceptions and numerical warnings without turning failure into p=1.
  safe <- function(fun) {
    warnings <- character()
    value <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
    if (inherits(value, "error"))
      return(list(value = NULL, status = paste0("error: ", conditionMessage(value))))
    if (is.null(value$p.value) || !is.finite(value$p.value))
      return(list(value = NULL, status = "undefined_test"))
    list(value = value, status = if (length(warnings))
      paste0("warning: ", paste(unique(warnings), collapse = "; ")) else "ok")
  }
  wide <- function(x, value, rows, columns) {
    # Use fixed internal names to support arbitrary (including non-syntactic) labels.
    mat <- matrix(x[[value]][NA_integer_], nrow = length(rows),
                  ncol = length(columns), dimnames = list(NULL, columns))
    if (nrow(x)) {
      ij <- cbind(match(x[[1L]], rows), match(x[[2L]], columns))
      mat[ij] <- x[[value]]
    }
    ans <- cbind(data.frame(rows, stringsAsFactors = FALSE), as.data.frame(mat,
                  optional = TRUE, stringsAsFactors = FALSE))
    names(ans) <- c(names(x)[1L], columns)
    data.table::as.data.table(ans, check.names = FALSE)
  }
  cell_result <- function(ct) {
    if (method == "chisq") {
      expected <- outer(rowSums(ct), colSums(ct)) / sum(ct)
      dimnames(expected) <- dimnames(ct)
      ratio <- matrix(NA_real_, nrow(ct), ncol(ct), dimnames = dimnames(ct))
      rr <- rowSums(ct) > 0; cc <- colSums(ct) > 0
      status <- "insufficient_positive_margins"; pv <- NA_real_
      if (sum(rr) >= 2L && sum(cc) >= 2L) {
        z <- safe(function() stats::chisq.test(ct[rr, cc, drop = FALSE]))
        status <- z$status
        if (!is.null(z$value)) {
          ratio[rr, cc] <- ct[rr, cc, drop = FALSE] / expected[rr, cc, drop = FALSE]
          pv <- z$value$p.value
        }
      }
      attr(ratio, "counts") <- ct; attr(ratio, "expected") <- expected
      attr(ratio, "p.value") <- pv; attr(ratio, "status") <- status
      attr(ratio, "test.warning") <- if (startsWith(status, "warning:")) status else ""
      return(ratio)
    }
    # Filtering selects reported targets only; background remains all cells.
    targets <- which(rowSums(ct) > min.rowSum)
    grid <- expand.grid(V1 = rownames(ct)[targets], V2 = colnames(ct),
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    grid$count <- numeric(nrow(grid)); grid$p.value <- grid$OR <- rep(NA_real_, nrow(grid))
    grid$status <- rep("zero_margin", nrow(grid))
    for (k in seq_len(nrow(grid))) {
      i <- match(grid$V1[k], rownames(ct)); j <- match(grid$V2[k], colnames(ct))
      a <- as.double(ct[i, j]); b <- sum(ct[i, ]) - a
      c <- sum(ct[, j]) - a; dd <- sum(ct) - a - b - c
      m <- matrix(c(a, b, c, dd), nrow = 2L)
      grid$count[k] <- a
      if (all(rowSums(m) > 0) && all(colSums(m) > 0)) {
        z <- safe(function() stats::fisher.test(m))
        grid$status[k] <- z$status
        if (!is.null(z$value)) {
          grid$p.value[k] <- z$value$p.value; grid$OR[k] <- unname(z$value$estimate)
        }
      }
    }
    grid$p.adj <- stats::p.adjust(grid$p.value, "BH")
    list(count.dist = data.table::as.data.table(grid),
         p.tb = wide(grid, "p.value", rownames(ct)[targets], colnames(ct)),
         padj.tb = wide(grid, "p.adj", rownames(ct)[targets], colnames(ct)),
         OR.tb = wide(grid, "OR", rownames(ct)[targets], colnames(ct)))
  }
  if (method != "freq") {
    warning("Cell-count analysis treats cells as observations; it does not establish population-level effects from biological replicates.", call. = FALSE)
    if (bySample) {
      samples <- sort(unique(d$sample))
      ans <- lapply(samples, function(s) cell_result(counts(d[d$sample == s, , drop = FALSE])))
      names(ans) <- samples
      if (any(vapply(samples, function(s) length(unique(d$condition[d$sample == s])) < 2L,
                     logical(1)))) warning("Some samples contain only one condition; within-sample tests are undefined (NA).", call. = FALSE)
      return(ans)
    }
    return(cell_result(counts(d)))
  }
  samples <- sort(unique(d$sample))
  nconditions <- vapply(samples, function(s) length(unique(d$condition[d$sample == s])), integer(1))
  if (any(nconditions != 1L)) stop("freq requires exactly one condition per sample; paired/repeated designs are unsupported")
  sc <- vapply(samples, function(s) d$condition[match(s, d$sample)], character(1))
  # Only observed samples; zeros mean a cluster was not recovered in that sample.
  ct <- table(factor(d$sample, levels = samples), factor(d$cluster, levels = clusters))
  sf <- expand.grid(sampleID = samples, majorCluster = clusters,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  sf$loc <- sc[match(sf$sampleID, samples)]
  sf$N <- as.vector(ct); sf$total <- as.numeric(rowSums(ct)[match(sf$sampleID, samples)])
  sf$freq.mcls <- sf$N / sf$total
  observed_conditions <- conditions[conditions %in% d$condition]
  res <- expand.grid(majorCluster = clusters, loc = observed_conditions,
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  for (nm in c("mean.x", "mean.y", "logFC", "p.value.t", "p.value.w")) res[[nm]] <- rep(NA_real_, nrow(res))
  res$n.x <- res$n.y <- integer(nrow(res))
  res$status.t <- res$status.w <- rep("inference_disabled", nrow(res))
  for (k in seq_len(nrow(res))) {
    f <- sf[sf$majorCluster == res$majorCluster[k], , drop = FALSE]
    x <- f$freq.mcls[f$loc == res$loc[k]]; y <- f$freq.mcls[f$loc != res$loc[k]]
    res$n.x[k] <- length(x); res$n.y[k] <- length(y)
    mx <- mean(x); my <- if (length(y)) mean(y) else NA_real_
    res$mean.x[k] <- mx; res$mean.y[k] <- my
    res$logFC[k] <- if (is.na(my) || (mx == 0 && my == 0)) NA_real_ else log2(mx / my)
    if (!sample.inference) next
    if (length(x) < 2L || length(y) < 2L) {
      res$status.t[k] <- res$status.w[k] <- "insufficient_replicates"; next
    }
    if (all(c(x, y) == 0)) {
      res$status.t[k] <- res$status.w[k] <- "all_zero"; next
    }
    for (test in c("t", "w")) {
      z <- safe(function() if (test == "t") stats::t.test(x, y) else
        stats::wilcox.test(x, y, exact = FALSE))
      res[[paste0("status.", test)]][k] <- z$status
      if (!is.null(z$value)) res[[paste0("p.value.", test)]][k] <- z$value$p.value
    }
  }
  if (sample.inference && any(res$status.t == "insufficient_replicates"))
    warning("At least two independent biological samples per contrast side are required; affected p-values are NA.", call. = FALSE)
  res$FDR.t <- stats::p.adjust(res$p.value.t, "BH")
  res$FDR.w <- stats::p.adjust(res$p.value.w, "BH")
  res$char.sig <- ifelse(is.na(res$FDR.t), "", ifelse(res$FDR.t < .01, "\u2731\u2731", ifelse(res$FDR.t < .05, "\u2731", "")))
  list(res = data.table::as.data.table(res),
       dist.logFC.tb = wide(res, "logFC", clusters, observed_conditions),
       dist.FDR.tb = wide(res, "FDR.t", clusters, observed_conditions),
       dist.charSig.tb = wide(res, "char.sig", clusters, observed_conditions),
       sample.freq.tb = data.table::as.data.table(sf))
}
