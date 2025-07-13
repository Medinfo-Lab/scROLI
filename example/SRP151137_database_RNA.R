library(readxl)
library(dplyr)
library(Seurat)
library(ggplot2)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggrepel)
library(enrichplot)



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
any(duplicated(SRP151137_RNA_data_counts$Geneid))
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


#seurat----
SRP151137_pbmc <- CreateSeuratObject(counts = SRP151137_RNA_data_counts_filter_choose,
                                min.cells = 100, min.features = 200)

SRP151137_RNA_sample_seurat_group <- SRP151137_RNA_sample %>%
  filter(Title %in% colnames(SRP151137_pbmc))
SRP151137_pbmc <- AddMetaData(object = SRP151137_pbmc,
                         metadata = SRP151137_RNA_sample_seurat_group$group1,
                         col.name = "group1")

SRP151137_pbmc <- NormalizeData(SRP151137_pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
SRP151137_pbmc <- FindVariableFeatures(SRP151137_pbmc, selection.method = "vst", nfeatures = 2000)
SRP151137_pbmc_variable_data <- VariableFeatures(SRP151137_pbmc)
SRP151137_pbmc <- ScaleData(SRP151137_pbmc)
SRP151137_pbmc <- RunPCA(SRP151137_pbmc, features = VariableFeatures(object = SRP151137_pbmc))
SRP151137_pbmc <- FindNeighbors(SRP151137_pbmc)
SRP151137_pbmc <- FindClusters(SRP151137_pbmc,resolution = 0.5)
SRP151137_pbmc <- RunUMAP(SRP151137_pbmc, dims = 1:50,
                          min.dist = 3,spread = 5)

DimPlot(SRP151137_pbmc,reduction = "umap",pt.size = 2.5,group.by = "group1")+
  ggtitle("Group UMAP")+
  theme(plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
        axis.title.x = element_text(face = "bold",size = 15),
        axis.title.y = element_text(face = "bold",size = 15))+
  scale_color_manual(values = c("#761C78","#5B6C2F","#4682B4"))

# DimPlot(SRP151137_pbmc,reduction = "umap",pt.size = 2.5,group.by = "seurat_clusters")+
#   ggtitle("Clusters UMAP")+
#   theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),
#         axis.title.x = element_text(face = "bold",size = 14),
#         axis.title.y = element_text(face = "bold",size = 14))
  # scale_color_manual(values=c("#264C58", "#30B3E9", "#4682B4"))


#DEG----
Idents(SRP151137_pbmc) <- SRP151137_pbmc$group1
# Idents(SRP151137_pbmc) <- "seurat_clusters"
SRP151137_pbmc_markers <- FindAllMarkers(SRP151137_pbmc)
# write.csv(SRP151137_pbmc_markers,"SRP151137_plot/RNA/SRP151137_RNA_DEG_markers.csv")

remove_after_last_dot <- function(x) {
  sub("\\.[^.]*$", "", x)  # 正则替换最后一次出现的点及后续内容
}

SRP151137_pbmc_markers$geneid <- remove_after_last_dot(SRP151137_pbmc_markers$gene)
# save.image("SRP151137_plot/RNA_markers.RData")

SRP151137_pbmc_markers_P <- SRP151137_pbmc_markers %>%
  filter(p_val < 0.05)

SRP151137_pbmc_markers_choose <- SRP151137_pbmc_markers_P %>%
  group_by(cluster) %>%
  top_n(n = 4000, wt = avg_log2FC)

# write.csv(SRP151137_pbmc_markers_choose,"SRP151137_plot/RNA/SRP151137_RNA_DEG_markers.csv")
# SRP151137_pbmc_markers_choose <- SRP151137_pbmc_markers %>%
#   group_by(cluster) %>%
#   top_n(n = 350, wt = avg_log2FC)

DoHeatmap(SRP151137_pbmc,features = SRP151137_pbmc_markers_choose$gene, angle = 20,
          group.colors = c("E4.5"="#30B3E9","E5.5"="#4682B4","E6.5"="#264C58"),
          label = T,slot = "scale.data")+
  scale_fill_gradient2(
    low = "#0571B0",
    mid = "#F7F7F7",
    high = "#CA0020",
    midpoint = 0)+
  theme(
    axis.text.y = element_blank(),  # 隐藏Y轴文本（基因名）
    axis.ticks.y = element_blank()   # 可选：隐藏Y轴刻度线
  )
  # scale_color_manual(values=c("#30B3E9", "#4682B4", "#264C58"))

# DotPlot(SRP151137_pbmc,features=SRP151137_pbmc_markers_choose$gene[10], split.by='group1')


volcano_data <- SRP151137_pbmc_markers %>%
  mutate(
    gene = rownames(SRP151137_pbmc_markers),  # 确保有基因名列
    log_pval = -log10(p_val),  # 转换校正p值
    direction = case_when(         # 标记上下调
      avg_log2FC > 0.5 & p_val < 0.05 ~ "Up",
      avg_log2FC < -0.5 & p_val < 0.05 ~ "Down",
      TRUE ~ "Not sig"
    )
  ) %>%
  # 过滤无限值（避免log(0)错误）
  filter(!is.infinite(log_pval))

top_genes_up <- volcano_data %>%
  filter(direction == "Up") %>%
  arrange(p_val) %>%
  head(5)

top_genes_down <- volcano_data %>%
  filter(direction == "Down") %>%
  arrange(p_val) %>%
  head(5)



ggplot(volcano_data, aes(x = avg_log2FC, y = log_pval, colour = direction)) +
  geom_point(alpha=0.5, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2", "#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(Fold Change)",y="-log10 (P-value)")+
  theme_bw()+
  ggtitle("RNA-seq DEG Volcano")+
  # 图例
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14))+
  xlim(-8, 8) + ylim(0, 50)
  geom_label_repel(data = top_genes_up, aes(label = gene),
                   size = 4,                           # 设置标签大小
                   box.padding = unit(0.8, "lines"),   # 设置标签内边距
                   point.padding = unit(0.8, "lines"), # 设置标签与点的距离
                   segment.color = "black",            # 设置标签边界线颜色
                   show.legend = FALSE,                # 不显示图例
                   max.overlaps = 10000)+             # 设置标签重叠的最大次数
  geom_label_repel(data = top_genes_down, aes(label = gene),
                   size = 4,                           # 设置标签大小
                   box.padding = unit(0.8, "lines"),   # 设置标签内边距
                   point.padding = unit(0.8, "lines"), # 设置标签与点的距离
                   segment.color = "black",            # 设置标签边界线颜色
                   show.legend = FALSE,                # 不显示图例
                   max.overlaps = 10000)               # 设置标签重叠的最大次数


#KEGG GO----
# enrichment_genes <- SRP151137_pbmc_markers_choose %>%
#   filter(cluster == "E4.5")

# enrichment_genes <- SRP151137_pbmc_markers_choose %>%
#   filter(cluster == "E5.5")

# enrichment_genes <- SRP151137_pbmc_markers_P %>%
#   filter(cluster == "E5.5")

enrichment_genes_choose <- SRP151137_pbmc_markers_P %>%
  # group_by(cluster) %>%
  top_n(n = 5000, wt = avg_log2FC)

# gene_ids <- bitr(enrichment_genes_choose$geneid,
#                  fromType="SYMBOL",
#                  toType="ENTREZID",
#                  OrgDb="org.Mm.eg.db")
# entrez_ids <- gene_ids$ENTREZID


#GO
# ego <- enrichGO(gene = entrez_ids,
#                 OrgDb = org.Mm.eg.db,
#                 ont = "ALL",        # "BP","MF","CC"或"ALL"
#                 pAdjustMethod = "BH",         # 校正方法：BH, bonferroni等
#                 pvalueCutoff = 0.05,         # 显著阈值
#                 qvalueCutoff = 0.05,
#                 readable = TRUE)         # 转换ID为基因名

ego <- enrichGO(gene = enrichment_genes_choose$geneid,
                keyType = "ENSEMBL",
                OrgDb = org.Mm.eg.db,
                ont = "ALL",        # "BP","MF","CC"或"ALL"
                pAdjustMethod = "BH",         # 校正方法：BH, bonferroni等
                pvalueCutoff = 0.05,         # 显著阈值
                qvalueCutoff = 0.05,
                readable = TRUE)         # 转换ID为基因名



#提取并简化结果（避免重复）
# go_simplified <- simplify(
#   ego,
#   cutoff = 0.1,           # 相似性阈值
#   by = "p.adjust",        # 按调整后p值筛选
#   select_fun = min
# )

ego_data <- ego@result

ego_data_group <- ego_data %>%
  group_by(GeneRatio)


ego_description <- c("mRNA processing",
                     "regulation of protein stability",
                     "histone modification",
                     "ncRNA processing",
                     "RNA splicing",
                     "protein localization to nucleus",
                     "mitotic cell cycle phase transition",
                     "autophagy",
                     "protein stabilization",
                     "non-membrane-bounded organelle assembly",
                     "double-strand break repair",
                     "DNA replication",
                     "RNA splicing, via transesterification reactions",
                     "nuclear division",
                     "gastrulation")


dotplot(ego,
        showCategory=ego_description,         # 显示top15条目
        color = "pvalue",
        title="RNA-seq DEG GO Enrichment") +
  # scale_color_gradient(low="#546de5", high="#ff4757")+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))+
  scale_color_gradientn(
    colours = c("#E64B35", "#4DBBD5", "#00A087"))



#KEGG
kk <- enrichKEGG(gene = entrez_ids,
                 organism = "mmu",       # 人类'hsa'，小鼠'mmu'
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)

dotplot(kk,
        showCategory=10,         # 显示top15条目
        # font.size=10,
        title="GO Enrichment") +
  scale_color_gradient(low="#546de5", high="#ff4757")+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))



#predit model----
SRP151137_pbmc_markers_choose_unique <- unique(SRP151137_pbmc_markers_choose$gene)
SRP151137_pbmc_NormalizeData <- GetAssayData(SRP151137_pbmc,slot = "data")
SRP151137_pbmc_NormalizeData <- as.data.frame(SRP151137_pbmc_NormalizeData)


SRP151137_RNA_LASSO_data <- SRP151137_pbmc_NormalizeData[SRP151137_pbmc_markers_choose_unique,]
SRP151137_RNA_LASSO_data_group <- SRP151137_pbmc$group1
SRP151137_RNA_LASSO_data_group <- as.data.frame(SRP151137_RNA_LASSO_data_group)


# RNA_LASSO_genes <- LASSO_feature(SRP151137_RNA_LASSO_data,SRP151137_RNA_LASSO_data_group)




