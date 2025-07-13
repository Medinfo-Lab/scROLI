# scMethyCA: Region-Based Integrated Epigenome-Transcriptome Analysis of Single-Cell

Unifying DNA Methylation, Chromatin Accessibility, and Transcriptome Profiling in Single-Cell Multi-Omics via Region-Centric Integration

# Installation instructions:

```R
#install devtools if you don't have it already for easy installation
install.packages("devtools")
library(devtools)
install_github("Medinfo-Lab/scMethyCA")
```

If you prefer to build the package by hand, follow these steps:

- Make sure that you have the dependencies from CRAN ("dply","reticulate","utils","png")

- Download and build from source:

```R
git clone git@github.com:Medinfo-Lab/scMethyCA.git
R CMD build scMethyCA
R CMD INSTALL scMethyCA-0.1.0.tar.gz
```

# Usage

![workflow](https://imgur.com/dYhHZyB.png)

```R
library(scMethyCA)

#First provide a coverage data, bed data and chromosome data.
merge_coverage <- list.files(
  coverage_path,
  full.names = TRUE,
  pattern = "\\.cov.gz$"
)

for (i in 1:length(merge_coverage)) {
  lines <- readLines(merge_coverage[i], warn = FALSE)
  # lines
  if (length(lines)<1) {
    message("跳过空文件: ", merge_coverage[i])
    next
  }
  
  cov_data <- read.table(merge_coverage[i])
  b <- data.frame()
  # 使用mclapply并行处理每个染色体
  merge_list <- mclapply(chromosome_data, function(chr_tmp) {
    merge_chr <- cov_to_data(merge_coverage[i], cov_data, chr_tmp, bed_data, suffixname)
    merge_name <- paste0("merge_", chr_tmp)
    return(list(merge_name = merge_chr))
  }, mc.cores = 10)  # 使用的CPU核心
  
  b <- do.call(rbind, lapply(merge_list, function(x) x$merge_name))
  sample_name <- colnames(b)
  gene_data_paste[, sample_name] <- b[, sample_name]
  
  cat("bed data ",merge_coverage[i],'\n')
  print(paste(i,Sys.time(),seq=""))
}
```









