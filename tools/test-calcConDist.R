# Run from repository root; isolates this function from unrelated Seurat/GO imports.
stopifnot(requireNamespace("data.table", quietly = TRUE),
          requireNamespace("testthat", quietly = TRUE))
source("R/mf_calcConDist.R")
testthat::test_file("tests/testthat/test-calcConDist.R", reporter = "summary",
                    stop_on_failure = TRUE, stop_on_warning = TRUE)
