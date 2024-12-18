% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mf_mean.R
\name{mf_mean}
\alias{mf_mean}
\title{Calculate Mean Values by Group for Matrix Data}
\usage{
mf_mean(mat, group, samples = colnames(mat), verbose = FALSE)
}
\arguments{
\item{mat}{A numeric matrix where rows represent features (e.g., genes) and columns represent samples.}

\item{group}{A vector specifying the group for each sample. Must be the same length as the number of columns in \code{mat}.}

\item{samples}{A vector of sample names. Default is \code{colnames(mat)}. If provided, must match the length of \code{group}.}

\item{verbose}{Logical. If TRUE, print progress messages. Default is TRUE.}

\item{na.rm}{Logical. Should NA values be removed when calculating means? Default is FALSE.}
}
\value{
A data frame where:
\itemize{
\item Each column represents a unique group
\item Each row represents a feature (e.g., gene)
\item Values are the mean of the features for each group
}
}
\description{
This function calculates the mean values of a matrix for each specified group.
It is particularly useful for summarizing high-dimensional data (e.g., gene expression data)
across different experimental conditions or groups.
}
\examples{
# Create a sample matrix (20 genes, 10 samples)
set.seed(123)
mat <- matrix(rnorm(200), nrow = 20, ncol = 10)
rownames(mat) <- paste0("Gene", 1:20)
colnames(mat) <- paste0("Sample", 1:10)

# Define groups
group <- rep(c("Control", "Treatment"), each = 5)

# Calculate mean values
result <- mf_mean(mat, group)
print(head(result))

}
