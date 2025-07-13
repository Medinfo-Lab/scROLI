library(readxl)
library(dplyr)
library(ggplot2)
library(scMethyCA)
library(dplyr)
library(MOFA2)
library(data.table)
library(org.Mm.eg.db)
library(ggVennDiagram)
library(ggraph)
library(igraph)
library(clusterProfiler)
library(VennDiagram)



# SMARCA5

#sample----
SRP472645_CpG_meth_level_data <- read.csv("data/SRP472645/SRP472645_CpGgene_meth_level.csv",row.names = 1)
SRP472645_GpC_meth_level_data <- read.csv("data/SRP472645/SRP472645_GpCgene_meth_level.csv",row.names = 1)

SRP472645_Epis_sample <- read.csv("data/SRP472645/sample_SRP472645_nomeseq.csv")
SRP472645_Epis_sample$Run_methlevel <- paste0(SRP472645_Epis_sample$Run, ".methlevel")
SRP472645_Epis_sample$Run_site <- paste0(SRP472645_Epis_sample$Run, ".site")

SRP472645_CpG_methlevel <- Read_file_meth_colname(SRP472645_CpG_meth_level_data,"methlevel")
SRP472645_GpC_methlevel <- Read_file_meth_colname(SRP472645_GpC_meth_level_data,"methlevel")

CpG_GpC_sample <- intersect(colnames(SRP472645_CpG_methlevel),colnames(SRP472645_GpC_methlevel))
SRP472645_Epis_sample_filter <- intersect(CpG_GpC_sample,SRP472645_Epis_sample$Run_methlevel)
SRP472645_Epis_sample_filter_choose <- SRP472645_Epis_sample %>%
  filter(Run_methlevel %in% SRP472645_Epis_sample_filter)

SRP472645_CpG_methlevel_filter <- SRP472645_CpG_methlevel[,factor(SRP472645_Epis_sample_filter_choose$Run_methlevel)]
SRP472645_GpC_methlevel_filter <- SRP472645_GpC_methlevel[,factor(SRP472645_Epis_sample_filter_choose$Run_methlevel)]

colnames(SRP472645_CpG_methlevel_filter) <- SRP472645_Epis_sample_filter_choose$Title
colnames(SRP472645_GpC_methlevel_filter) <- SRP472645_Epis_sample_filter_choose$Title



#chr data----
GRCm38_data <- read.csv("data/SRP151137/GRCm38_Genes.csv")
GRCm38_data_genes <- read.csv("data/SRP151137/GRCm38_Genes_genename.csv")
GRCm38_data_paste <- Chr_region_process(GRCm38_data,"paste")
GRCm38_data_paste$gene <- GRCm38_data_genes$gene_name

# any(duplicated(GRCm38_data_paste$gene))





#NA process----
SRP472645_CpG_methlevel_na_count <- NA_count_order(SRP472645_CpG_methlevel_filter)
SRP472645_GpC_methlevel_na_count <- NA_count_order(SRP472645_GpC_methlevel_filter)

SRP472645_CpG_methlevel_na_count_nona <- SRP472645_CpG_methlevel_na_count %>%
  filter(NA_percentage<=0.5)

SRP472645_CpG_methlevel_choose <- SRP472645_CpG_methlevel_filter[rownames(SRP472645_CpG_methlevel_na_count_nona),]

SRP472645_GpC_methlevel_na_count_nona <- SRP472645_GpC_methlevel_na_count %>%
  filter(NA_percentage<=0.9)

SRP472645_GpC_methlevel_choose <- SRP472645_GpC_methlevel_filter[rownames(SRP472645_GpC_methlevel_na_count_nona),]


#CpG mofa----
SRP472645_CpG_methlevel_choose_mean <- NA_padding_mean(SRP472645_CpG_methlevel_choose)

SRP472645_CpG_methlevel_choose_mean_matrix <- as.matrix(SRP472645_CpG_methlevel_choose_mean)
mofa_list <- list()
mofa_list$file1 <- SRP472645_CpG_methlevel_choose_mean_matrix
SRP472645_CpG_mofa <- MOFA_est(mofa_list,group = SRP472645_Epis_sample_filter_choose$group1)

# SRP472645_CpG_mofa <- run_umap(SRP472645_CpG_mofa,min_dist = 1,spread = 2)
SRP472645_CpG_mofa <- run_umap(SRP472645_CpG_mofa)

plot_dimred(SRP472645_CpG_mofa, method = "UMAP",color_by = "group",dot_size = 2,label = F)+
  theme(legend.title=element_blank())

SRP472645_CpG_umap_data <- SRP472645_CpG_mofa@dim_red$UMAP
SRP472645_CpG_umap_data <- cbind(SRP472645_CpG_umap_data,SRP472645_CpG_mofa@samples_metadata$group)
colnames(SRP472645_CpG_umap_data)[4] <- "group"
rownames(SRP472645_CpG_umap_data) <- 1:nrow(SRP472645_CpG_umap_data)


ggplot(SRP472645_CpG_umap_data)+
  geom_point(SRP472645_CpG_umap_data,mapping=aes(UMAP1,UMAP2,color=group),size = 2.5,alpha = 0.7)+
  theme_classic()+
  labs(x = "UMAP1",y = "UMAP2",title = "CpG UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14))


#GpC mofa----
SRP472645_GpC_methlevel_choose_mean <- NA_padding_mean(SRP472645_GpC_methlevel_choose)

SRP472645_GpC_methlevel_choose_mean_matrix <- as.matrix(SRP472645_GpC_methlevel_choose_mean)
mofa_list <- list()



#CpG DEG----
SRP472645_CpG_DEG <- Meth_group_variance_analysis(SRP472645_CpG_methlevel_choose_mean,SRP472645_Epis_sample_filter_choose$group1,"var.test")

SRP472645_CpG_DEG_P <- SRP472645_CpG_DEG %>%
  filter(adj_P.value < 0.05)

SRP472645_CpG_DEG_choose <- SRP472645_CpG_DEG_P %>%
  group_by(group) %>%
  top_n(n = 400, wt = logFC)

SRP472645_CpG_DEG_choose_chr <- unique(SRP472645_CpG_DEG_P$chr)
SRP472645_CpG_DEG_choose_chr <- as.data.frame(SRP472645_CpG_DEG_choose_chr)
colnames(SRP472645_CpG_DEG_choose_chr) <- "chr"

SRP472645_CpG_DEG_choose_chr_gene <- inner_join(SRP472645_CpG_DEG_choose_chr, GRCm38_data_paste, by = "chr") %>%
  dplyr::select(chr, gene)  # 选择所需列


#CpG DEG enrichment----
SRP472645_CpG_DEG_gene <- bitr(SRP472645_CpG_DEG_choose_chr_gene$gene,
                               fromType = "SYMBOL",
                               toType = "ENTREZID",
                               OrgDb = org.Mm.eg.db)

ego <- enrichGO(gene = SRP472645_CpG_DEG_gene$ENTREZID,
                OrgDb = org.Mm.eg.db,
                ont = "ALL",        # "BP","MF","CC"或"ALL"
                pAdjustMethod = "BH",         # 校正方法：BH, bonferroni等
                pvalueCutoff = 0.05,         # 显著阈值
                qvalueCutoff = 0.05,
                readable = TRUE)         # 转换ID为基因名

dotplot(ego,
        showCategory=10,         # 显示top15条目
        # font.size=10,
        title="GO Enrichment") +
  scale_color_gradient(low="#546de5", high="#ff4757")+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))






#GpC DEG----
SRP472645_GpC_DEG <- Meth_group_variance_analysis(SRP472645_GpC_methlevel_choose_mean,SRP472645_Epis_sample_filter_choose$group1,"var.test")

SRP472645_GpC_DEG_P <- SRP472645_GpC_DEG %>%
  filter(adj_P.value < 0.05)

SRP472645_GpC_DEG_choose <- SRP472645_GpC_DEG_P %>%
  group_by(group) %>%
  top_n(n = 5, wt = logFC)

SRP472645_GpC_DEG_choose_chr <- unique(SRP472645_GpC_DEG_P$chr)
SRP472645_GpC_DEG_choose_chr <- as.data.frame(SRP472645_GpC_DEG_choose_chr)
colnames(SRP472645_GpC_DEG_choose_chr) <- "chr"

SRP472645_GpC_DEG_choose_chr_gene <- inner_join(SRP472645_GpC_DEG_choose_chr, GRCm38_data_paste, by = "chr") %>%
  select(chr, gene)  # 选择所需列


#GpC DEG enrichment----
# valid_keys <- keys(org.Mm.eg.db, keytype="SYMBOL")
# invalid_genes <- setdiff(SRP472645_GpC_DEG_choose_chr_gene$gene, valid_keys)
#
# SRP472645_GpC_DEG_gene <- bitr(SRP472645_GpC_DEG_choose_chr_gene$gene,
#                                fromType = "SYMBOL",
#                                toType = "ENTREZID",
#                                OrgDb = org.Mm.eg.db)
#
# ego <- enrichGO(gene = SRP472645_GpC_DEG_gene$ENTREZID,
#                 OrgDb = org.Mm.eg.db,
#                 ont = "ALL",        # "BP","MF","CC"或"ALL"
#                 pAdjustMethod = "BH",         # 校正方法：BH, bonferroni等
#                 pvalueCutoff = 0.05,         # 显著阈值
#                 qvalueCutoff = 0.05,
#                 readable = TRUE)         # 转换ID为基因名




#CpG mean meth group----
# SRP472645_CpG_methlevel_DEG <- SRP472645_CpG_methlevel_filter[SRP472645_CpG_DEG_choose_chr$chr,]
SRP472645_CpG_methlevel_SMARCA5 <- SRP472645_CpG_methlevel_choose["chr8:80737497-80741497",]

SRP472645_CpG_methlevel_choose_colmean <- colMeans(SRP472645_CpG_methlevel_SMARCA5,na.rm = T)
SRP472645_CpG_methlevel_choose_colmean <- as.data.frame(SRP472645_CpG_methlevel_choose_colmean)
SRP472645_CpG_methlevel_choose_colmean$group <- SRP472645_Epis_sample_filter_choose$group1

SRP472645_CpG_methlevel_choose_colmean$group <- factor(SRP472645_CpG_methlevel_choose_colmean$group)

SRP472645_CpG_methlevel_choose_colmean_mean_data <- SRP472645_CpG_methlevel_choose_colmean %>%
  group_by(group) %>%
  summarise(mean_value = mean(SRP472645_CpG_methlevel_choose_colmean, na.rm = TRUE))

ggplot(SRP472645_CpG_methlevel_choose_colmean,
       aes(x = group, y = SRP472645_CpG_methlevel_choose_colmean, fill = group)) +
  geom_violin(alpha = 0.8,scale = "width")+
  geom_segment(data = SRP472645_CpG_methlevel_choose_colmean_mean_data,
               aes(x = as.numeric(group) - 0.1,
                   xend = as.numeric(group) + 0.1,
                   y = mean_value, yend = mean_value),
               color = "red", size = 0.5, linetype = "solid")+
  theme_bw()+
  labs(y="Percentage (%)",
       x="",
       title = "CpG ColMean")+
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 15,hjust = 0.5))+
  # scale_fill_manual(values = group_colors)+
  guides(fill="none")



#GpC mean meth group----
# SRP472645_GpC_methlevel_DEG <- SRP472645_GpC_methlevel_filter[SRP472645_GpC_DEG_choose_chr$chr,]
SRP472645_GpC_methlevel_SMARCA5 <- SRP472645_GpC_methlevel_filter["chr8:80737497-80741497",]

SRP472645_GpC_methlevel_choose_colmean <- colMeans(SRP472645_GpC_methlevel_SMARCA5,na.rm = T)
SRP472645_GpC_methlevel_choose_colmean <- as.data.frame(SRP472645_GpC_methlevel_choose_colmean)
SRP472645_GpC_methlevel_choose_colmean$group <- SRP472645_Epis_sample_filter_choose$group1

SRP472645_GpC_methlevel_choose_colmean$group <- factor(SRP472645_GpC_methlevel_choose_colmean$group)

SRP472645_GpC_methlevel_choose_colmean_mean_data <- SRP472645_GpC_methlevel_choose_colmean %>%
  group_by(group) %>%
  summarise(mean_value = mean(SRP472645_GpC_methlevel_choose_colmean, na.rm = TRUE))

ggplot(SRP472645_GpC_methlevel_choose_colmean,
       aes(x = group, y = SRP472645_GpC_methlevel_choose_colmean, fill = group)) +
  geom_violin(alpha = 0.8,scale = "width")+
  geom_segment(data = SRP472645_GpC_methlevel_choose_colmean_mean_data,
               aes(x = as.numeric(group) - 0.1,
                   xend = as.numeric(group) + 0.1,
                   y = mean_value, yend = mean_value),
               color = "red", size = 0.5, linetype = "solid")+
  theme_bw()+
  labs(y="Percentage (%)",
       x="",
       title = "GpC ColMean")+
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        plot.title = element_text(size = 15,hjust = 0.5))+
  # scale_fill_manual(values = group_colors)+
  guides(fill="none")


#Venn genes----
genes <- list(
  RNA_Gene = SRP472646_pbmc_markers_choose$gene,
  CpG_Gene = SRP472645_CpG_DEG_choose_chr_gene$gene,
  GpC_Gene = SRP472645_GpC_DEG_choose_chr_gene$gene
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
