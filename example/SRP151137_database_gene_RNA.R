library(Seurat)
library(dplyr)
library(scMethyCA)


#chr data----
GRCm38_data <- read.csv("data/list/GRCm38_Genes.csv")
GRCm38_data_genes <- read.csv("data/list/GRCm38_Genes_genename.csv")
GRCm38_data_paste <- Chr_region_process(GRCm38_data,"paste")
GRCm38_data_paste$geneid <- GRCm38_data_genes$gene_id


#RNA----
SRP151137_RNA_data <- read.csv("data/SRP151137/SRP151137_RNA_counts_genename.csv")
SRP151137_RNA_sample <- read.csv("data/SRP151137/sample/SRP151137_RNA_CpGGpC_sample_data.csv")
# SRP151137_CpG_GpC_sample <- read.csv("data/SRP151137/sample/CpG_GpC_sample_data.csv")

# CpG_GpC_RNA_sample <- intersect(SRP151137_CpG_GpC_sample$Title,SRP151137_RNA_sample$Title)
# SRP151137_RNA_sample_choose <- SRP151137_RNA_sample %>%
#   filter(Title %in% CpG_GpC_RNA_sample)
# write.csv(SRP151137_RNA_sample_choose,"data/SRP151137/sample/SRP151137_RNA_CpGGpC_sample_data.csv")
# SRP151137_CpG_GpC_sample_choose <- SRP151137_CpG_GpC_sample %>%
#   filter(Title %in% CpG_GpC_RNA_sample)
# write.csv(SRP151137_CpG_GpC_sample_choose,"data/SRP151137/sample/SRP151137_CpGGpC_RNA_sample_data.csv")


SRP151137_RNA_data_length <- SRP151137_RNA_data[,c("Geneid","Length")]
SRP151137_RNA_data_length <-  unique(SRP151137_RNA_data_length, by = c("Geneid", "Length"))

SRP151137_RNA_data_counts <- SRP151137_RNA_data[,c(-3,-2)]
# any(duplicated(SRP151137_RNA_data_counts$Geneid))
# any(is.na(SRP151137_RNA_data$Genename))

dedup_dt <- SRP151137_RNA_data_counts[!duplicated(SRP151137_RNA_data_counts$Geneid), ]
sample_cols <- 2:ncol(dedup_dt)

# 计算每个基因在多少样本中表达量>3
n_samples_above_threshold <- rowSums(dedup_dt[, sample_cols] > 3)
# 保留至少满足3个样本的基因
SRP151137_RNA_data_counts_filter <- dedup_dt[n_samples_above_threshold >= 3, ]

rownames(SRP151137_RNA_data_counts_filter) <- SRP151137_RNA_data_counts_filter$Geneid
SRP151137_RNA_data_counts_filter <- SRP151137_RNA_data_counts_filter[,-1]


SRP151137_RNA_data_counts_filter_choose <- SRP151137_RNA_data_counts_filter[,SRP151137_RNA_sample$Run]
colnames(SRP151137_RNA_data_counts_filter_choose) <- SRP151137_RNA_sample$Title

remove_after_last_dot <- function(x) {
  sub("\\.[^.]*$", "", x)  # 正则替换最后一次出现的点及后续内容
}


#seurat----
SRP151137_pbmc <- CreateSeuratObject(counts = SRP151137_RNA_data_counts_filter_choose,
                                     min.cells = 100, min.features = 0)

SRP151137_RNA_sample_seurat_group <- SRP151137_RNA_sample %>%
  filter(Title %in% colnames(SRP151137_pbmc))
SRP151137_pbmc <- AddMetaData(object = SRP151137_pbmc,
                              metadata = SRP151137_RNA_sample_seurat_group$group1,
                              col.name = "group1")

SRP151137_RNA_data_relevance <- GetAssayData(SRP151137_pbmc,slot = "count")
SRP151137_RNA_data_relevance <- as.data.frame(SRP151137_RNA_data_relevance)

rownames(SRP151137_RNA_data_relevance) <- remove_after_last_dot(rownames(SRP151137_RNA_data_relevance))




#CpG GpC gene----
SRP151137_CpG_meth_level_data <- read.csv("data/SRP151137/SRP151137_CpGgene_meth_level.csv",row.names = 1)
SRP151137_GpC_meth_level_data <- read.csv("data/SRP151137/SRP151137_GpCgene_meth_level.csv",row.names = 1)



GRCm38_data_paste$geneid <- remove_after_last_dot(GRCm38_data_paste$geneid)

ensemblid_common <- intersect(rownames(SRP151137_RNA_data_relevance),GRCm38_data_paste$geneid)
# a <- remove_after_last_dot(SRP151137_RNA_data_counts$Geneid)
# ensemblid_common <- intersect(a,GRCm38_data_paste$geneid)

rownames(SRP151137_CpG_meth_level_data) <- GRCm38_data_paste$geneid
rownames(SRP151137_GpC_meth_level_data) <- GRCm38_data_paste$geneid

SRP151137_Epis_sample <- read.csv("data/SRP151137/sample/SRP151137_CpGGpC_RNA_sample_data.csv")
# SRP151137_Epis_sample$Run_methlevel <- paste0(SRP151137_Epis_sample$Run, ".methlevel")

SRP151137_CpG_methlevel <- Read_file_meth_colname(SRP151137_CpG_meth_level_data,"methlevel")
SRP151137_GpC_methlevel <- Read_file_meth_colname(SRP151137_GpC_meth_level_data,"methlevel")

# CpG_GpC_sample <- intersect(colnames(SRP151137_CpG_methlevel),colnames(SRP151137_GpC_methlevel))

SRP151137_CpG_methlevel_filter <- SRP151137_CpG_methlevel[,factor(SRP151137_Epis_sample$Run_methlevel)]
SRP151137_GpC_methlevel_filter <- SRP151137_GpC_methlevel[,factor(SRP151137_Epis_sample$Run_methlevel)]

# any(colnames(SRP151137_CpG_methlevel_filter)==colnames(SRP151137_GpC_methlevel_filter))
# SRP151137_Epis_sample_filter <- SRP151137_Epis_sample %>%
#   filter(Run_methlevel %in% CpG_GpC_sample)

colnames(SRP151137_CpG_methlevel_filter) <- SRP151137_Epis_sample$Title
colnames(SRP151137_GpC_methlevel_filter) <- SRP151137_Epis_sample$Title

SRP151137_CpG_methlevel_filter_ensembl_commonid <- SRP151137_CpG_methlevel_filter[ensemblid_common,]
SRP151137_GpC_methlevel_filter_ensembl_commonid <- SRP151137_GpC_methlevel_filter[ensemblid_common,]
SRP151137_RNA_data_relevance_ensembl_commonid <- SRP151137_RNA_data_relevance[ensemblid_common,]

# write.csv(SRP151137_CpG_methlevel_filter_ensembl_commonid,"a.csv")
# write.csv(SRP151137_RNA_data_relevance_ensembl_commonid,"b.csv")

SRP151137_CpG_methlevel_na_count <- NA_count_order(SRP151137_CpG_methlevel_filter_ensembl_commonid)
SRP151137_GpC_methlevel_na_count <- NA_count_order(SRP151137_GpC_methlevel_filter_ensembl_commonid)


SRP151137_CpG_methlevel_na_count_nona <- SRP151137_CpG_methlevel_na_count %>%
  filter(NA_percentage<=0.65) #0.45

SRP151137_CpG_methlevel_choose <- SRP151137_CpG_methlevel_filter_ensembl_commonid[rownames(SRP151137_CpG_methlevel_na_count_nona),]

SRP151137_GpC_methlevel_na_count_nona <- SRP151137_GpC_methlevel_na_count %>%
  filter(NA_percentage<=0.9) #0.85

SRP151137_GpC_methlevel_choose <- SRP151137_GpC_methlevel_filter[rownames(SRP151137_GpC_methlevel_na_count_nona),]

# SRP151137_CpG_methlevel_filter_ensembl_commonid_t <- t(SRP151137_CpG_methlevel_filter_ensembl_commonid)
# SRP151137_RNA_data_relevance_ensembl_commonid_t <- t(SRP151137_RNA_data_relevance_ensembl_commonid)

SRP151137_CpG_methlevel_choose_mean <- NA_padding_mean(SRP151137_CpG_methlevel_choose)

# all(is.numeric(SRP151137_CpG_methlevel_filter_ensembl_commonid[1,]))


result <- cor.test(SRP151137_CpG_methlevel_choose[1,],
                   SRP151137_RNA_data_relevance_ensembl_commonid[1,],
                   method = "pearson",use = "complete.obs")





