# Can also be run without installing the package: Rscript tools/test-calcConDist.R
make_cells <- function(counts, sample, condition) {
  data.frame(cell = rep(names(counts), counts), lib = sample, sex = condition,
             stringsAsFactors = FALSE)
}
run_dist <- function(d, ...) calcConDist(d, colname.cluster = "cell",
  colname.sample = "lib", colname.condition = "sex", ...)
quiet_dist <- function(...) suppressWarnings(run_dist(...))
fixture <- function() rbind(
  make_cells(c(A = 2, B = 8), "m1", "M"),
  make_cells(c(A = 8, B = 2), "m2", "M"),
  make_cells(c(A = 1, B = 9), "f1", "F"),
  make_cells(c(A = 3, B = 7), "f2", "F"),
  make_cells(c(A = 4, B = 6), "f3", "F"))

test_that("pooled chi-square agrees with stats and keeps matrix API", {
  d <- fixture(); ct <- table(d$cell, d$sex)
  expect_warning(r <- run_dist(d), "Cell-count")
  ref <- suppressWarnings(chisq.test(ct))
  expect_true(is.matrix(r))
  expect_equal(as.vector(r), as.vector(ref$observed / ref$expected))
  expect_equal(attr(r, "p.value"), ref$p.value)
  expect_equal(unname(attr(r, "counts")), unname(ct))
})

test_that("Fisher 2x2 targets, orientation, BH and nonnumeric labels are correct", {
  d <- fixture(); d$cell[d$cell == "A"] <- "oil body"
  r <- quiet_dist(d, method = "fisher")
  ct <- table(d$cell, d$sex)
  for (k in seq_len(nrow(r$count.dist))) {
    row <- as.data.frame(r$count.dist)[k, ]; i <- row$V1; j <- row$V2
    a <- ct[i, j]; m <- matrix(c(a, sum(ct[i, ]) - a,
      sum(ct[, j]) - a, sum(ct) - sum(ct[i, ]) - sum(ct[, j]) + a), 2)
    ref <- fisher.test(m)
    expect_equal(row$p.value, ref$p.value)
    expect_equal(row$OR, unname(ref$estimate))
  }
  expect_equal(r$count.dist$p.adj, p.adjust(r$count.dist$p.value, "BH"))
  expect_named(r, c("count.dist", "p.tb", "padj.tb", "OR.tb"))
  expect_s3_class(r$OR.tb, "data.table")
})

test_that("bySample dispatches before pooled returns and preserves named lists", {
  d <- fixture()
  # A sample spanning conditions allows a within-sample association test.
  d$lib <- ifelse(d$lib %in% c("m1", "f1"), "pair 1", "pair 2")
  for (m in c("chisq", "fisher")) {
    ans <- quiet_dist(d, method = m, bySample = TRUE)
    expect_named(ans, c("pair 1", "pair 2"))
    for (s in names(ans)) expect_equal(ans[[s]], quiet_dist(d[d$lib == s, ], method = m))
  }
})

test_that("single-condition samples are undefined, not pooled or goodness-of-fit tests", {
  d <- fixture()
  expect_warning(expect_warning(ans <- run_dist(d, bySample = TRUE),
                                "within-sample"), "Cell-count")
  expect_length(ans, 5L)
  expect_true(all(vapply(ans, function(z) all(is.na(z)), logical(1))))
  expect_true(all(vapply(ans, function(z) attr(z, "status") ==
    "insufficient_positive_margins", logical(1))))
  fish <- quiet_dist(d, method = "fisher", bySample = TRUE)
  expect_true(all(vapply(fish, function(z) all(is.na(z$count.dist$p.value)), logical(1))))
  expect_true(all(vapply(fish, function(z) all(z$count.dist$status == "zero_margin"), logical(1))))
})

test_that("zero factor margins and filtering retain coherent shapes and backgrounds", {
  d <- fixture(); d$cell <- factor(d$cell, levels = c("A", "B", "absent"))
  d$sex <- factor(d$sex, levels = c("F", "M", "unused"))
  z <- quiet_dist(d)
  expect_equal(dim(z), c(3L, 3L))
  expect_true(all(is.na(z["absent", ])))
  expect_true(all(is.na(z[, "unused"])))
  all <- quiet_dist(d, method = "fisher")
  selected <- quiet_dist(d, method = "fisher", min.rowSum = 20)
  ref <- as.data.frame(all$count.dist)
  ref <- ref[ref$V1 == "B", ]; rownames(ref) <- NULL
  expect_equal(as.data.frame(selected$count.dist), ref)
  empty <- quiet_dist(d, method = "fisher", min.rowSum = 1e6)
  expect_equal(nrow(empty$count.dist), 0L)
  expect_equal(dim(empty$p.tb), c(0L, 4L))
  expect_equal(quiet_dist(d[1, ], method = "fisher")$count.dist$status[1], "zero_margin")
})

test_that("freq handles custom non-syntactic input names and does not mutate inputs", {
  d <- data.table::as.data.table(fixture()); before <- data.table::copy(d)
  r <- run_dist(d, method = "freq")
  expect_identical(d, before)
  dd <- as.data.frame(d); names(dd) <- c("cell type", "sample id", "group name")
  rr <- calcConDist(dd, colname.cluster = "cell type", colname.sample = "sample id",
    colname.condition = "group name", method = "freq", bySample = TRUE)
  expect_equal(rr, r)
  expect_named(r, c("res", "dist.logFC.tb", "dist.FDR.tb", "dist.charSig.tb", "sample.freq.tb"))
  expect_true(all(is.na(r$res$p.value.t)))
  expect_true(all(r$res$status.t == "inference_disabled"))
})

test_that("freq tests use equal-weight independent samples with unequal group sizes", {
  d <- fixture(); r <- run_dist(d, method = "freq", sample.inference = TRUE)
  a <- as.data.frame(r$res); a <- a[a$majorCluster == "A" & a$loc == "M", ]
  x <- c(.2, .8); y <- c(.1, .3, .4)
  expect_equal(a$mean.x, mean(x)); expect_equal(a$mean.y, mean(y))
  expect_equal(a$logFC, log2(mean(x) / mean(y)))
  expect_equal(a$n.x, 2L); expect_equal(a$n.y, 3L)
  expect_equal(a$p.value.t, t.test(x, y)$p.value)
  expect_equal(a$p.value.w, wilcox.test(x, y, exact = FALSE)$p.value)
  expect_equal(r$res$FDR.t, p.adjust(r$res$p.value.t, "BH"))
  expect_equal(r$res$FDR.w, p.adjust(r$res$p.value.w, "BH"))
  sums <- tapply(r$sample.freq.tb$freq.mcls, r$sample.freq.tb$sampleID, sum)
  expect_equal(as.numeric(sums), rep(1, 5))
  # Unequal cell yields must not weight the condition mean.
  d2 <- rbind(d, d[d$lib == "m1", ], d[d$lib == "m1", ])
  expect_equal(run_dist(d2, method = "freq")$res$mean.x,
               run_dist(d, method = "freq")$res$mean.x)
})

test_that("freq fills missing clusters but not unobserved sample-condition combinations", {
  d <- fixture(); d <- d[!(d$lib == "m1" & d$cell == "A"), ]
  d$cell <- factor(d$cell, levels = c("B", "A", "zero"))
  d$lib <- factor(d$lib, levels = c(unique(d$lib), "absent sample"))
  d$sex <- factor(d$sex, levels = c("M", "F", "absent condition"))
  r <- run_dist(d, method = "freq")
  sf <- as.data.frame(r$sample.freq.tb)
  expect_equal(nrow(sf), 5L * 3L)
  expect_equal(sf$N[sf$sampleID == "m1" & sf$majorCluster == "A"], 0L)
  expect_equal(sf$freq.mcls[sf$sampleID == "m1" & sf$majorCluster == "A"], 0)
  expect_false("absent condition" %in% r$res$loc)
  expect_equal(nrow(r$res), 6L)
  expect_true(all(is.na(r$res$logFC[r$res$majorCluster == "zero"])))
})

test_that("zero denominators and constant proportions are explicit without pseudocounts", {
  d <- rbind(make_cells(c(A = 10), "m1", "M"), make_cells(c(A = 10), "m2", "M"),
             make_cells(c(B = 10), "f1", "F"), make_cells(c(B = 10), "f2", "F"))
  d$cell <- factor(d$cell, levels = c("A", "B", "zero"))
  r <- run_dist(d, method = "freq", sample.inference = TRUE)
  a <- as.data.frame(r$res); am <- a[a$majorCluster == "A" & a$loc == "M", ]
  af <- a[a$majorCluster == "A" & a$loc == "F", ]
  expect_equal(am$logFC, Inf); expect_equal(af$logFC, -Inf)
  expect_true(is.na(am$p.value.t)); expect_match(am$status.t, "error:")
  expect_true(is.finite(am$p.value.w))
  expect_true(all(a$status.t[a$majorCluster == "zero"] == "all_zero"))
  d$cell <- "A"
  r <- run_dist(d, method = "freq", sample.inference = TRUE)
  expect_true(all(is.na(r$res$p.value.w)))
  expect_true(all(r$res$status.w == "undefined_test"))
})

test_that("no replicates (Tak1/Tak2 pattern) gives descriptive output only", {
  d <- rbind(make_cells(c(A = 2, B = 8), "Tak1_11d", "male"),
             make_cells(c(A = 6, B = 4), "Tak2_11d", "female"))
  r <- run_dist(d, method = "freq")
  expect_true(all(is.na(r$res$p.value.t)))
  expect_true(all(is.finite(r$res$logFC)))
  expect_warning(rr <- run_dist(d, method = "freq", sample.inference = TRUE), "two independent")
  expect_true(all(rr$res$status.t == "insufficient_replicates"))
  expect_true(all(is.na(rr$res$FDR.w)))
  # Having 3 vs 1 also fails: both sides are checked.
  d <- fixture(); d <- d[d$lib != "m2", ]
  r <- quiet_dist(d, method = "freq", sample.inference = TRUE)
  expect_true(all(r$res$status.t == "insufficient_replicates"))
})

test_that("one condition, one cluster and more than two conditions are stable", {
  d <- fixture(); d$sex <- "M"
  r <- run_dist(d, method = "freq")
  expect_true(all(r$res$n.y == 0L)); expect_true(all(is.na(r$res$logFC)))
  d <- fixture(); d$sex[d$lib == "f3"] <- "X"
  r <- quiet_dist(d, method = "freq", sample.inference = TRUE)
  expect_equal(nrow(r$res), 6L)
  expect_true(all(r$res$status.t[r$res$loc == "X"] == "insufficient_replicates"))
  expect_equal(r$res$mean.y[r$res$majorCluster == "A" & r$res$loc == "M"], mean(c(.1,.3,.4)))
  d$cell <- "A"
  expect_true(all(is.na(quiet_dist(d))))
})

test_that("invalid metadata and options fail early", {
  d <- fixture()
  expect_error(run_dist(d, method = "other"), "Unsupported")
  expect_error(run_dist(d, method = c("freq", "chisq")), "Unsupported")
  expect_error(run_dist(d, bySample = NA), "logicals")
  expect_error(run_dist(d, min.rowSum = -1), "min.rowSum")
  expect_error(run_dist(d, sample.inference = TRUE), "only supported")
  expect_error(run_dist(d[0, ]), "nonempty")
  expect_error(run_dist(d[, -1]), "Missing columns")
  expect_error(calcConDist(d, colname.cluster = "cell", colname.condition = "sex",
                          method = "freq"), "colname.sample")
  for (bad in list(NA_character_, "", "  ")) {
    dd <- d; dd$cell[1] <- bad
    expect_error(run_dist(dd), "metadata")
  }
  dd <- d; dd$sex[dd$lib == "m1"][1] <- "F"
  expect_error(run_dist(dd, method = "freq"), "one condition per sample")
  dd <- d; names(dd)[2] <- "cell"
  expect_error(run_dist(dd), "Duplicate")
})

test_that("labels colliding with output id names or separators retain all values", {
  d <- fixture(); d$sex[d$sex == "M"] <- "majorCluster"
  d$cell[d$cell == "A"] <- "a\rb"
  r <- run_dist(d, method = "freq")
  expect_equal(names(r$dist.logFC.tb), c("majorCluster", "F", "majorCluster"))
  expect_equal(ncol(r$dist.logFC.tb), 3L)
  d$sex[d$sex == "majorCluster"] <- "V1"
  f <- quiet_dist(d, method = "fisher")
  expect_equal(names(f$p.tb), c("V1", "F", "V1"))
})
