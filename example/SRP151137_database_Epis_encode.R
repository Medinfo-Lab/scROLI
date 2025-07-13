library(readxl)
library(dplyr)
library(ggplot2)
library(scMethyCA)
library(dplyr)
library(MOFA2)
library(data.table)
library(org.Mm.eg.db)
library(ggVennDiagram)
library(VennDiagram)
library(clusterProfiler)
library(enrichplot)






#epis data----
SRP151137_CpG_meth_level_data <- read.csv("data/SRP151137/SRP151137_CpGencode_meth_level.csv",row.names = 1)
SRP151137_GpC_meth_level_data <- read.csv("data/SRP151137/SRP151137_GpCencode_meth_level.csv",row.names = 1)

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



#encode----
mm10_encode <- read.csv("data/list//mm10_encode.csv")

# 处理 gene 列：保留每个记录的第一个基因名
# mm10_encode$gene <- sapply(mm10_encode$gene, function(g) {
#   # 处理 "NONE" 或空值
#   if (is.na(g) || g == "NONE") return("NONE")
#
#   # 分割字符串：按逗号分隔
#   genes <- unlist(strsplit(g, "\\s*,\\s*"))
#
#   # 提取第一个基因名（移除基因距离信息）
#   first_gene <- sub("\\s*\\(.*", "", genes[1])
#
#   return(first_gene)
# })

# mm10_encode_paste <- Chr_region_process(mm10_encode,"paste")
# mm10_encode_paste$gene <- mm10_encode$gene


#NA process----
SRP151137_CpG_methlevel_na_count <- NA_count_order(SRP151137_CpG_methlevel_filter)
SRP151137_GpC_methlevel_na_count <- NA_count_order(SRP151137_GpC_methlevel_filter)

SRP151137_CpG_methlevel_na_count_nona <- SRP151137_CpG_methlevel_na_count %>%
  filter(NA_percentage<=0.8) #0.7

SRP151137_CpG_methlevel_choose <- SRP151137_CpG_methlevel_filter[rownames(SRP151137_CpG_methlevel_na_count_nona),]

SRP151137_GpC_methlevel_na_count_nona <- SRP151137_GpC_methlevel_na_count %>%
  filter(NA_percentage<=0.98)

SRP151137_GpC_methlevel_choose <- SRP151137_GpC_methlevel_filter[rownames(SRP151137_GpC_methlevel_na_count_nona),]


# save.image("SRP151137_plot/ENCODE_ALL_Means.RData")
load("SRP151137_plot/ENCODE_ALL_Means.RData")


#CpG mofa----
SRP151137_CpG_methlevel_choose_mean <- NA_padding_mean(SRP151137_CpG_methlevel_choose)

SRP151137_CpG_methlevel_choose_mean_matrix <- as.matrix(SRP151137_CpG_methlevel_choose_mean)
mofa_list <- list()
mofa_list$file1 <- SRP151137_CpG_methlevel_choose_mean_matrix
# SRP151137_CpG_mofa <- MOFA_est(mofa_list,group = SRP151137_Epis_sample_filter$group1)

SRP151137_CpG_mofa <- create_mofa(mofa_list,group = SRP151137_Epis_sample$group1)

ModelOptions <- get_default_model_options(SRP151137_CpG_mofa)
TrainOptions <- get_default_training_options(SRP151137_CpG_mofa)
DataOptions <- get_default_data_options(SRP151137_CpG_mofa)

TrainOptions$convergence_mode <- "slow"
DataOptions$scale_views <- TRUE
# DataOptions$center_groups <- TRUE



SRP151137_CpG_mofa <- prepare_mofa(SRP151137_CpG_mofa,
                                   model_options = ModelOptions,
                                   training_options = TrainOptions,
                                   data_options = DataOptions)

SRP151137_CpG_mofa <- run_mofa(SRP151137_CpG_mofa)

# SRP151137_CpG_mofa <- run_umap(SRP151137_CpG_mofa,min_dist = 1,spread = 2)
SRP151137_CpG_mofa <- run_umap(SRP151137_CpG_mofa,
                               min_dist = 2,spread = 3)

plot_dimred(SRP151137_CpG_mofa, method = "UMAP",color_by = "group",dot_size = 2,label = F)+
  theme(legend.title=element_blank())

SRP151137_CpG_umap_data <- SRP151137_CpG_mofa@dim_red$UMAP
SRP151137_CpG_umap_data <- cbind(SRP151137_CpG_umap_data,SRP151137_CpG_mofa@samples_metadata$group)
colnames(SRP151137_CpG_umap_data)[4] <- "group"
rownames(SRP151137_CpG_umap_data) <- 1:nrow(SRP151137_CpG_umap_data)

# write.csv(SRP151137_CpG_umap_data,"SRP151137_plot/encode_data/SRP151137_CpGencode_umap_data.csv")

ggplot(SRP151137_CpG_umap_data)+
  geom_point(SRP151137_CpG_umap_data,mapping=aes(UMAP1,UMAP2,color=group),size = 2.5,alpha = 0.8)+
  theme_classic()+
  labs(x = "UMAP1",y = "UMAP2",title = "ENCODE CpG UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5,size = 18,face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14))+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5"),
                     values = c("#761C78","#5B6C2F","#4682B4"))


#GpC mofa----
SRP151137_GpC_methlevel_choose_mean <- NA_padding_mean(SRP151137_GpC_methlevel_choose)

SRP151137_GpC_methlevel_choose_mean_matrix <- as.matrix(SRP151137_GpC_methlevel_choose_mean)
mofa_list <- list()
mofa_list$file1 <- SRP151137_GpC_methlevel_choose_mean_matrix
# SRP151137_GpC_mofa <- MOFA_est(mofa_list,group = SRP151137_Epis_sample_filter$group1)

SRP151137_GpC_mofa <- create_mofa(mofa_list,group = SRP151137_Epis_sample$group1)

ModelOptions <- get_default_model_options(SRP151137_GpC_mofa)
TrainOptions <- get_default_training_options(SRP151137_GpC_mofa)
DataOptions <- get_default_data_options(SRP151137_GpC_mofa)

TrainOptions$convergence_mode <- "slow"
DataOptions$scale_views <- TRUE
# DataOptions$center_groups <- TRUE



SRP151137_GpC_mofa <- prepare_mofa(SRP151137_GpC_mofa,
                                   model_options = ModelOptions,
                                   training_options = TrainOptions,
                                   data_options = DataOptions)

SRP151137_GpC_mofa <- run_mofa(SRP151137_GpC_mofa)

# SRP151137_GpC_mofa <- run_umap(SRP151137_GpC_mofa,min_dist = 1,spread = 2)
SRP151137_GpC_mofa <- run_umap(SRP151137_GpC_mofa,
                               min_dist = 2,spread = 3)

plot_dimred(SRP151137_GpC_mofa, method = "UMAP",color_by = "group",dot_size = 2,label = F)+
  theme(legend.title=element_blank())

SRP151137_GpC_umap_data <- SRP151137_GpC_mofa@dim_red$UMAP
SRP151137_GpC_umap_data <- cbind(SRP151137_GpC_umap_data,SRP151137_GpC_mofa@samples_metadata$group)
colnames(SRP151137_GpC_umap_data)[4] <- "group"
rownames(SRP151137_GpC_umap_data) <- 1:nrow(SRP151137_GpC_umap_data)

# write.csv(SRP151137_GpC_umap_data,"SRP151137_plot/encode_data/SRP151137_GpCencode_umap_data.csv")

ggplot(SRP151137_GpC_umap_data)+
  geom_point(SRP151137_GpC_umap_data,mapping=aes(UMAP1,UMAP2,color=group),size = 2.5,alpha = 0.8)+
  theme_classic()+
  labs(x = "UMAP1",y = "UMAP2",title = "ENCODE GpC UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5,size = 18,face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14))+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5"),
                     values = c("#761C78","#5B6C2F","#4682B4"))


#CpG DEG----
# a <- SRP151137_Epis_sample_filter$group1
SRP151137_CpG_DEG <- Meth_group_variance_analysis(SRP151137_CpG_methlevel_choose_mean,SRP151137_Epis_sample$group1)
# write.csv(SRP151137_CpG_DEG,"SRP151137_plot/encode_data/SRP151137_CpG_DEG.csv")

SRP151137_CpG_DEG_P <- SRP151137_CpG_DEG %>%
  filter(P.value < 0.05)
# write.csv(SRP151137_CpG_DEG_P,"SRP151137_plot/encode_data/SRP151137_CpG_DEG_P0.05.csv")


# SRP151137_CpG_DEG_choose <- SRP151137_CpG_DEG_P %>%
#   group_by(group) %>%
#   top_n(n = 500, wt = logFC)

SRP151137_CpG_DEG_choose_chr <- unique(SRP151137_CpG_DEG_P$chr)
SRP151137_CpG_DEG_choose_chr <- as.data.frame(SRP151137_CpG_DEG_choose_chr)
colnames(SRP151137_CpG_DEG_choose_chr) <- "chr"

# SRP151137_CpG_DEG_choose_chr_gene <- inner_join(SRP151137_CpG_DEG_choose_chr, mm10_enhancers_paste, by = "chr") %>%
#   dplyr::select(chr, gene)  # 选择所需列


volcano_data <- SRP151137_CpG_DEG %>%
  mutate(
    # chr = rownames(SRP151137_CpG_DEG_P$chr),  # 确保有基因名列
    log_pval = -log10(P.value),  # 转换校正p值
    direction = case_when(         # 标记上下调
      logFC > 0.5 & P.value < 0.05 ~ "Up",
      logFC < -0.5 & P.value < 0.05 ~ "Down",
      TRUE ~ "Not sig"
    )
  ) %>%
  # 过滤无限值（避免log(0)错误）
  filter(!is.infinite(log_pval))

ggplot(volcano_data, aes(x = logFC, y = log_pval, colour = direction)) +
  geom_point(alpha=0.8, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2", "#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(Fold Change)",y="-log10 (P-value)")+
  theme_bw()+
  ggtitle("CpG ENCODE DEG Volcano")+
  # 图例
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5,size = 18,face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14))
  # xlim(-12, 12) + ylim(0, 35)





#CpG DEG enrichment----
SRP151137_CpG_DEG_gene <- bitr(SRP151137_CpG_DEG_choose_chr_gene$gene,
                               fromType = "SYMBOL",
                               toType = "ENTREZID",
                               OrgDb = org.Mm.eg.db)

ego <- enrichGO(gene = SRP151137_CpG_DEG_gene$ENTREZID,
                OrgDb = org.Mm.eg.db,
                ont = "ALL",        # "BP","MF","CC"或"ALL"
                pAdjustMethod = "BH",         # 校正方法：BH, bonferroni等
                pvalueCutoff = 0.05,         # 显著阈值
                qvalueCutoff = 0.05,
                readable = TRUE)         # 转换ID为基因名

enrichplot::dotplot(ego,
                    showCategory=10,         # 显示top15条目
                    # font.size=10,
                    title="GO Enrichment") +
  scale_color_gradient(low="#546de5", high="#ff4757")+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))


#GpC DEG----
SRP151137_GpC_DEG <- Meth_group_variance_analysis(SRP151137_GpC_methlevel_choose_mean,SRP151137_Epis_sample$group1)
# write.csv(SRP151137_GpC_DEG,"SRP151137_plot/encode_data/SRP151137_GpC_DEG.csv")

SRP151137_GpC_DEG_P <- SRP151137_GpC_DEG %>%
  filter(P.value < 0.05)

SRP151137_GpC_DEG_choose <- SRP151137_GpC_DEG_P %>%
  group_by(group) %>%
  top_n(n = 5, wt = logFC)

SRP151137_GpC_DEG_choose_chr <- unique(SRP151137_GpC_DEG_P$chr)
SRP151137_GpC_DEG_choose_chr <- as.data.frame(SRP151137_GpC_DEG_choose_chr)
colnames(SRP151137_GpC_DEG_choose_chr) <- "chr"

SRP151137_GpC_DEG_choose_chr_gene <- inner_join(SRP151137_GpC_DEG_choose_chr, GRCm38_data_paste, by = "chr") %>%
  select(chr, gene)


volcano_data <- SRP151137_GpC_DEG %>%
  mutate(
    # chr = rownames(SRP151137_CpG_DEG_P$chr),  # 确保有基因名列
    log_pval = -log10(P.value),  # 转换校正p值
    direction = case_when(         # 标记上下调
      logFC > 0.5 & P.value < 0.05 ~ "Up",
      logFC < -0.5 & P.value < 0.05 ~ "Down",
      TRUE ~ "Not sig"
    )
  ) %>%
  # 过滤无限值（避免log(0)错误）
  filter(!is.infinite(log_pval))

ggplot(volcano_data, aes(x = logFC, y = log_pval, colour = direction)) +
  geom_point(alpha=0.8, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2", "#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(Fold Change)",y="-log10 (P-value)")+
  theme_bw()+
  ggtitle("GpC ENCODE DEG Volcano")+
  # 图例
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5,size = 18,face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14))
# xlim(-12, 12) + ylim(0, 35)


#GpC DEG enrichment----
# valid_keys <- keys(org.Mm.eg.db, keytype="SYMBOL")
# invalid_genes <- setdiff(SRP151137_GpC_DEG_choose_chr_gene$gene, valid_keys)
#
# SRP151137_GpC_DEG_gene <- bitr(SRP151137_GpC_DEG_choose_chr_gene$gene,
#                                fromType = "SYMBOL",
#                                toType = "ENTREZID",
#                                OrgDb = org.Mm.eg.db)
#
# ego <- enrichGO(gene = SRP151137_GpC_DEG_gene$ENTREZID,
#                 OrgDb = org.Mm.eg.db,
#                 ont = "ALL",        # "BP","MF","CC"或"ALL"
#                 pAdjustMethod = "BH",         # 校正方法：BH, bonferroni等
#                 pvalueCutoff = 0.05,         # 显著阈值
#                 qvalueCutoff = 0.05,
#                 readable = TRUE)         # 转换ID为基因名




#CpG mean meth group----
# SRP151137_CpG_methlevel_DEG <- SRP151137_CpG_methlevel_filter[SRP151137_CpG_DEG_choose_chr$chr,]

SRP151137_CpG_methlevel_choose_colmean <- colMeans(SRP151137_CpG_methlevel_choose,na.rm = T)
SRP151137_CpG_methlevel_choose_colmean <- as.data.frame(SRP151137_CpG_methlevel_choose_colmean)
SRP151137_CpG_methlevel_choose_colmean$group <- SRP151137_Epis_sample$group1

SRP151137_CpG_methlevel_choose_colmean$group <- factor(SRP151137_CpG_methlevel_choose_colmean$group)

SRP151137_CpG_methlevel_choose_colmean_mean_data <- SRP151137_CpG_methlevel_choose_colmean %>%
  group_by(group) %>%
  summarise(mean_value = mean(SRP151137_CpG_methlevel_choose_colmean, na.rm = TRUE))

ggplot(SRP151137_CpG_methlevel_choose_colmean,
       aes(x = group, y = SRP151137_CpG_methlevel_choose_colmean, fill = group)) +
  geom_violin(alpha = 0.8,scale = "width")+
  geom_segment(data = SRP151137_CpG_methlevel_choose_colmean_mean_data,
               aes(x = as.numeric(group) - 0.1,
                   xend = as.numeric(group) + 0.1,
                   y = mean_value, yend = mean_value),
               color = "red", size = 0.5, linetype = "solid")+
  theme_bw()+
  labs(y="Percentage (%)",
       x="",
       title = "ENCODE CpG ColMean")+
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 12,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 14,face = "bold"),
        plot.title = element_text(size = 18,hjust = 0.5,face = "bold"))+
  # scale_fill_manual(values = group_colors)+
  guides(fill="none")



#GpC mean meth group----
# SRP151137_GpC_methlevel_DEG <- SRP151137_GpC_methlevel_filter[SRP151137_GpC_DEG_choose_chr$chr,]

SRP151137_GpC_methlevel_choose_colmean <- colMeans(SRP151137_GpC_methlevel_choose,na.rm = T)
SRP151137_GpC_methlevel_choose_colmean <- as.data.frame(SRP151137_GpC_methlevel_choose_colmean)
SRP151137_GpC_methlevel_choose_colmean$group <- SRP151137_Epis_sample$group1

SRP151137_GpC_methlevel_choose_colmean$group <- factor(SRP151137_GpC_methlevel_choose_colmean$group)

SRP151137_GpC_methlevel_choose_colmean_mean_data <- SRP151137_GpC_methlevel_choose_colmean %>%
  group_by(group) %>%
  summarise(mean_value = mean(SRP151137_GpC_methlevel_choose_colmean, na.rm = TRUE))

ggplot(SRP151137_GpC_methlevel_choose_colmean,
       aes(x = group, y = SRP151137_GpC_methlevel_choose_colmean, fill = group)) +
  geom_violin(alpha = 0.8,scale = "width")+
  geom_segment(data = SRP151137_GpC_methlevel_choose_colmean_mean_data,
               aes(x = as.numeric(group) - 0.1,
                   xend = as.numeric(group) + 0.1,
                   y = mean_value, yend = mean_value),
               color = "red", size = 0.5, linetype = "solid")+
  theme_bw()+
  labs(y="Percentage (%)",
       x="",
       title = "ENCODE GpC ColMean")+
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 12,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 14,face = "bold"),
        plot.title = element_text(size = 18,hjust = 0.5,face = "bold"))+
  # scale_fill_manual(values = group_colors)+
  guides(fill="none")


#Venn genes----
genes <- list(
  RNA_Gene = SRP151137_pbmc_markers_choose$gene,
  CpG_Gene = SRP151137_CpG_DEG_choose_chr_gene$gene,
  GpC_Gene = SRP151137_GpC_DEG_choose_chr_gene$gene
)

ggVennDiagram(genes,
              label_alpha = 0,  # 标签透明背景
              edge_size = 0.5) +
  ggplot2::scale_fill_gradient(low = "white", high = "firebrick") +
  ggplot2::labs(title = "Gene Overlap")+
  theme(
    plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=16,face="bold")
    # axis.title.x = element_text(size = 14,face = "bold"),
    # axis.title.y = element_text(size = 14,face = "bold")
  )

inter <- get.venn.partitions(genes)
for (i in 1:nrow(inter)) {
  inter[i,'Genes'] <- paste(inter[[i,'..values..']], collapse=", ")
}

venn_genes <- inter$Genes[5]
venn_genes <- unlist(strsplit(venn_genes, ","))


venn_DEG_genes <- bitr(venn_genes,
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Mm.eg.db)

ego <- enrichGO(gene = venn_DEG_genes$ENTREZID,
                OrgDb = org.Mm.eg.db,
                ont = "ALL",        # "BP","MF","CC"或"ALL"
                pAdjustMethod = "BH",         # 校正方法：BH, bonferroni等
                pvalueCutoff = 0.05,         # 显著阈值
                qvalueCutoff = 0.05,
                readable = TRUE)         # 转换ID为基因名

dotplot(ego,
        showCategory=10,         # 显示top15条目
        # font.size=10,
        title="Venn Gene GO Enrichment") +
  scale_color_gradient(low="#546de5", high="#ff4757")+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
