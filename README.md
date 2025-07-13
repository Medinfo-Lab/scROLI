# scMethyCA

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
R CMD INSTALL scMethyCA_0.1.0.tar.gz
```

# Usage

![](data\workflow.jpg)
