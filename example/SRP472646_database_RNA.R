library(readxl)
library(dplyr)
library(Seurat)
library(ggplot2)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggrepel)
library(enrichplot)
library(clusterProfiler)


# SMARCA5


#RNA----
SRP472646_RNA_data <- read.csv("data/SRP472646/SRP472646_RNA_counts_genename.csv")
SRP472646_RNA_sample <- read.csv("data/SRP472646/sample_SRP472646_scRNAseq.csv")


SRP472646_RNA_data_length <- SRP472646_RNA_data[,c("Genename","Length")]
SRP472646_RNA_data_length <-  unique(SRP472646_RNA_data_length, by = c("Genename", "Length"))

SRP472646_RNA_data_counts <- SRP472646_RNA_data[,c(-3,-1)]
any(duplicated(SRP472646_RNA_data_counts$Genename))

dedup_dt <- SRP472646_RNA_data_counts[!duplicated(SRP472646_RNA_data_counts$Genename), ]
sample_cols <- 2:ncol(dedup_dt)


# 计算每个基因在多少样本中表达量>3
n_samples_above_threshold <- rowSums(dedup_dt[, sample_cols] > 3)
# 保留至少满足3个样本的基因
SRP472646_RNA_data_counts_filter <- dedup_dt[n_samples_above_threshold >= 3, ]

rownames(SRP472646_RNA_data_counts_filter) <- SRP472646_RNA_data_counts_filter$Genename
SRP472646_RNA_data_counts_filter <- SRP472646_RNA_data_counts_filter[,-1]


SRP472646_RNA_data_counts_filter_choose <- SRP472646_RNA_data_counts_filter[,SRP472646_RNA_sample$Run]
colnames(SRP472646_RNA_data_counts_filter_choose) <- SRP472646_RNA_sample$Title

#seurat----
SRP472646_pbmc <- CreateSeuratObject(counts = SRP472646_RNA_data_counts_filter_choose,
                                     min.cells = 100, min.features = 200)

SRP472646_RNA_sample_seurat_group <- SRP472646_RNA_sample %>%
  filter(Title %in% colnames(SRP472646_pbmc))
SRP472646_pbmc <- AddMetaData(object = SRP472646_pbmc,
                              metadata = SRP472646_RNA_sample_seurat_group$group1,
                              col.name = "group1")

SRP472646_pbmc <- NormalizeData(SRP472646_pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
SRP472646_pbmc <- FindVariableFeatures(SRP472646_pbmc, selection.method = "vst", nfeatures = 2000)
SRP472646_pbmc_variable_data <- VariableFeatures(SRP472646_pbmc)
SRP472646_pbmc <- ScaleData(SRP472646_pbmc)
SRP472646_pbmc <- RunPCA(SRP472646_pbmc, features = VariableFeatures(object = SRP472646_pbmc))
SRP472646_pbmc <- FindNeighbors(SRP472646_pbmc)
SRP472646_pbmc <- FindClusters(SRP472646_pbmc,resolution = 0.5)
SRP472646_pbmc <- RunUMAP(SRP472646_pbmc, dims = 1:15,
                          min.dist = 2,spread = 3)

DimPlot(SRP472646_pbmc,reduction = "umap",pt.size = 2.5,group.by = "group1")+
  ggtitle("Group UMAP")+
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14))
  # scale_color_manual(values=c("#264C58", "#30B3E9", "#4682B4"))

DimPlot(SRP472646_pbmc,reduction = "umap",pt.size = 2.5,group.by = "seurat_clusters")+
  ggtitle("Clusters UMAP")+
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14))
# scale_color_manual(values=c("#264C58", "#30B3E9", "#4682B4"))


Idents(SRP472646_pbmc) <- SRP472646_pbmc$group1


VlnPlot(SRP472646_pbmc, features = c("Smarca2","Smarca4","Smarca5"))
FeaturePlot(SRP472646_pbmc, features = c("Smarca2","Smarca4","Smarca5"), ncol = 3)
RidgePlot(SRP472646_pbmc, features = c("Smarca2","Smarca4","Smarca5"), ncol = 3)
DotPlot(SRP472646_pbmc, features = c("Smarca2","Smarca4","Smarca5"))



#DEG----
Idents(SRP472646_pbmc) <- SRP472646_pbmc$group1
# Idents(SRP472646_pbmc) <- "seurat_clusters"
SRP472646_pbmc_markers <- FindAllMarkers(SRP472646_pbmc)

SRP472646_pbmc_markers_P <- SRP472646_pbmc_markers %>%
  filter(p_val < 0.05)

SRP472646_pbmc_markers_choose <- SRP472646_pbmc_markers_P %>%
  group_by(cluster) %>%
  top_n(n = 400, wt = avg_log2FC)

# SRP472646_pbmc_markers_choose <- SRP472646_pbmc_markers %>%
#   group_by(cluster) %>%
#   top_n(n = 350, wt = avg_log2FC)

DoHeatmap(SRP472646_pbmc,features = SRP472646_pbmc_markers_choose$gene, angle = 15,
          # group.by = SRP472646_pbmc$seurat_clusters,
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

# DotPlot(SRP472646_pbmc,features=SRP472646_pbmc_markers_choose$gene[10], split.by='group1')


volcano_data <- SRP472646_pbmc_markers %>%
  mutate(
    gene = rownames(SRP472646_pbmc_markers),  # 确保有基因名列
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
  head(4)

top_genes_down <- volcano_data %>%
  filter(direction == "Down") %>%
  arrange(p_val) %>%
  head(4)



ggplot(volcano_data, aes(x = avg_log2FC, y = log_pval, colour = direction)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2", "#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(Fold Change)",y="-log10 (P-value)")+
  theme_bw()+
  ggtitle("Volcano")+
  # 图例
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14))+
  xlim(-5, 5) + ylim(0, 20)+
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
IgG_genes <- SRP472646_pbmc_markers_choose %>%
  filter(cluster == "IgG")

enrichment_genes_choose <- IgG_genes %>%
  # group_by(cluster) %>%
  top_n(n = 1000, wt = avg_log2FC)

gene_ids <- bitr(enrichment_genes_choose$gene,
                 fromType="SYMBOL",
                 toType="ENTREZID",
                 OrgDb="org.Mm.eg.db")
entrez_ids <- gene_ids$ENTREZID


#GO
ego <- enrichGO(gene = entrez_ids,
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

showCategory_data <- c("chromosome segregation","condensed chromosome","condensed chromosome, centromeric region",
                       "iron ion homeostasis","meiosis I cell cycle process",
                       "DNA-templated transcription initiation","spermatid nucleus differentiation")

dotplot(ego,
        showCategory=showCategory_data,         # 显示top15条目
        # font.size=10,
        title="GO Enrichment") +
  scale_color_gradient(low="#546de5", high="#ff4757")+
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
        showCategory=1,         # 显示top15条目
        # font.size=10,
        title="GO Enrichment") +
  scale_color_gradient(low="#546de5", high="#ff4757")+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))

#predit model----
