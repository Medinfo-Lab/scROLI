#' distinguish chrXX:xxx-xxx as chrXX, xxx, xxx
#'
#' @param file_tmp Chromosome physical segment data, chrXX:xxx-xxx format or chrXX, xxx, xxx format
#'
#' @return Processed physical fragment data, chrXX, xxx, xxx format or chrXX:xxx-xxx format
#' @export
#'
#' @examples
Chr_region_process <- function(file_tmp,method){
  if(method=="split"){
    chr_data <- file_tmp
    chr_data_frame <- data.frame(region=chr_data)
    chr_split <- chr_data_frame %>%
      separate(region, into = c("chr", "start", "end"), sep = ":|-")
    chr_split$start <- as.numeric(chr_split$start)
    chr_split$end <- as.numeric(chr_split$end)
    return(chr_split)
  }else if(method=="paste"){
    chr_data <- file_tmp
    chr_data_paste <- sprintf("%s:%d-%d",chr_data$chr,chr_data$start,chr_data$end)
    chr_data_paste_frame <- as.data.frame(chr_data_paste)
    colnames(chr_data_paste_frame) <- "chr"
    return(chr_data_paste_frame)
  }
}


#' Group analysis
#'
#' @param file_data DNA methylation or chromatin accessibility level data
#' @param method Statistical methods
#' @param file_group_tmp Group Information
#'
#' @return Difference analysis results
#' @export
#'
#' @examples
Meth_group_variance_analysis <- function(file_data,file_group_tmp,method="t.test"){
  # file_data <- SRP151137_CpG_methlevel_choose_mean
  # file_group_tmp <- SRP151137_Epis_sample_filter$group1
  # method <- "var.test"

  file_group <- file_group_tmp
  file_group <- factor(file_group)

  file_group_factor <- levels(file_group)
  file_sample <- data.table()
  file_sample$sample <- colnames(file_data)
  file_sample$group <- file_group
  file_sample <- as.data.frame(file_sample)

  DEG_data <- data.table()

  for (i in 1:length(file_group_factor)) {
    file_group_tmp <- file_sample %>%
      filter(group == file_group_factor[i])
    file_group_data <- file_data[,file_group_tmp$sample]
    file_group_data <- as.matrix(file_group_data)

    file_nogroup_tmp <- file_sample %>%
      filter(group != file_group_factor[i])
    file_nogroup_data <- file_data[,file_nogroup_tmp$sample]
    file_nogroup_data <- as.matrix(file_nogroup_data)

    if (method=="t.test") {
      t.test_data <- data.frame(chr = rownames(file_group_data))
      for (j in 1:nrow(file_group_data)) {
        DEG_tmp <- t.test(file_group_data[j,],file_nogroup_data[j,])
        t.test_data$logFC[j] <- log2(mean(file_group_data[j,],na.rm = T)/mean(file_nogroup_data[j,],na.rm = T))
        t.test_data$P.value[j] <- DEG_tmp$p.value
        t.test_data$adj_P.value[j] <- p.adjust(DEG_tmp$p.value,method = "BH")
        t.test_data$group[j] <- file_group_factor[i]
      }
      DEG_data <- rbind(DEG_data,t.test_data)
    }

    if(method=="var.test"){
      var.test_data <- data.frame(chr = rownames(file_group_data))
      for (j in 1:nrow(file_group_data)) {
        DEG_tmp <- var.test(file_group_data[j,],file_nogroup_data[j,])
        var.test_data$logFC[j] <- log2(mean(file_group_data[j,],na.rm = T)/mean(file_nogroup_data[j,],na.rm = T))
        var.test_data$P.value[j] <- DEG_tmp$p.value
        var.test_data$adj_P.value[j] <- p.adjust(DEG_tmp$p.value,method = "BH")
        var.test_data$group[j] <- file_group_factor[i]
      }
      DEG_data <- rbind(DEG_data,var.test_data)
    }

    # if(method=="fisher.test"){
    #   fisher.test_data <- data.frame(chr = rownames(file_group_data))
    #   for (j in 1:nrow(file_group_data)) {
    #     DEG_tmp <- var.test(file_group_data[j,],file_nogroup_data[j,])
    #     fisher.test_data$logFC[j] <- log2(mean(file_group_data[j,],na.rm = T)/mean(file_nogroup_data[j,],na.rm = T))
    #     fisher.test_data$P.value[j] <- DEG_tmp$p.value
    #     fisher.test_data$adj_P.value[j] <- p.adjust(DEG_tmp$p.value,method = "BH")
    #     fisher.test_data$group[j] <- file_group_factor[i]
    #   }
    #   DEG_data <- rbind(DEG_data,fisher.test_data)
    # }
  }
  return(DEG_data)
}



