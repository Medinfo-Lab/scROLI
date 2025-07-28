# scROLI

# Installation instructions

*scROLI: An R Package for* ***R**egion-**O**rdered Multi-Omics **L**inkage and **I**ntegration for **S**ingle **C**ells*

*scROLI* is an R software package for the joint analysis of transcriptome, DNA methylation and chromatin open data. The package is designed to process sequencing data for NOMe-seq and can also process BS-seq data. In addition, scROLI can convert sequencing data into region-based data for easy storage and subsequent analysis.

```R
#install devtools if you don't have it already for easy installation
install.packages("devtools")
library(devtools)
install_github("Medinfo-Lab/scROLI")
```

If you prefer to build the package by hand, follow these steps:

- Make sure that you have the dependencies from CRAN ("dply","reticulate","utils","png")

- Download and build from source:

```R
git clone git@github.com:Medinfo-Lab/scROLI.git
R CMD build scMethyCA
R CMD INSTALL scMethyCA-0.1.0.tar.gz
```

# Usage

**The Workflow:**

![](https://imgur.com/CO92ngt.png)

**The Epigenomic Processing Flow:**

![](https://imgur.com/GQ0WaVf.png)

**The Transcriptomic Processing Flow:**

![](https://imgur.com/gVY0VJ0.png)

```R
library(scMethyCA)
library(dply)
library(data.table)

#First provide a coverage data, bed data and chromosome data
merge_coverage <- list.files(
  coverage_path,
  full.names = TRUE,
  pattern = "\\.cov.gz$"
)

load("data/list.RData")
bed_data <- bed_data %>% 
	filter(chr %in% chromosome_data)
bed_data_paste <- as.data.frame(sprintf("%s:%s-%s", bed_data$chr, bed_data$start,
                                        bed_data$end))
colnames(bed_data_paste) <- "chr"
bed_data_paste_methlevel <- bed_data_paste
bed_data_paste_meth <- bed_data_paste

#site methlevel
for (i in 1:length(merge_coverage)) {
  start_time <- Sys.time()
  lines <- readLines(merge_coverage[i], warn = FALSE)
  # lines
  if (length(lines)<1) {
    message("Skip empty files: ", merge_coverage[i])
    next
  }
  
  cov_data <- fread(merge_coverage[i])
  b <- data.frame()
  # Parallel processing of each chromosome using mclapply
  merge_list <- mclapply(chromosome_data, function(chr_tmp) {
    merge_chr <- cov_to_data(merge_coverage[i], cov_data, chr_tmp, bed_data, suffixname, "methlevel")
    merge_name <- paste0("merge_", chr_tmp)
    return(list(merge_name = merge_chr))
  }, mc.cores = 10)  # Number of CPU cores used
  
  b <- do.call(rbind, lapply(merge_list, function(x) x$merge_name))
  sample_name <- colnames(b)
  bed_data_paste_methlevel[, sample_name] <- b[, sample_name]
  
  cat(i,"bed data ",merge_coverage[i],'\n')
  cat("文件处理耗时:", round(Sys.time() - start_time, 1), "秒\n")
}
                             
#meth UNmeth
for (i in 1:length(merge_coverage)) {
  start_time <- Sys.time()
  lines <- readLines(merge_coverage[i], warn = FALSE)
  # lines
  if (length(lines)<1) {
    message("Skip empty files: ", merge_coverage[i])
    next
  }
  
  cov_data <- fread(merge_coverage[i])
  b <- data.frame()
  # Parallel processing of each chromosome using mclapply
  merge_list <- mclapply(chromosome_data, function(chr_tmp) {
    merge_chr <- cov_to_data(merge_coverage[i], cov_data, chr_tmp, bed_data, suffixname, "meth")
    merge_name <- paste0("merge_", chr_tmp)
    return(list(merge_name = merge_chr))
  }, mc.cores = 10)  # Number of CPU cores used
  
  b <- do.call(rbind, lapply(merge_list, function(x) x$merge_name))
  sample_name <- colnames(b)
  bed_data_paste_meth[, sample_name] <- b[, sample_name]
  
  cat(i,"bed data ",merge_coverage[i],'\n')
  cat("文件处理耗时:", round(Sys.time() - start_time, 1), "秒\n")
}
```







