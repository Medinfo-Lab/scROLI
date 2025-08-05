library(ggtext)
library(ggrepel)
library(fpc)
library(cluster)
library(clusterProfiler)
library(readxl)
library(dplyr)
library(ggplot2)
library(sciTEA)
library(MOFA2)
library(data.table)
library(org.Mm.eg.db)
library(ggVennDiagram)
library(VennDiagram)
library(enrichplot)




#epis data----
SRP151137_CpG_meth_level_data <- read.csv("data/SRP151137/SRP151137_CpGgene_methlevel.csv",row.names = 1)
SRP151137_GpC_meth_level_data <- read.csv("data/SRP151137/SRP151137_GpCgene_methlevel.csv",row.names = 1)
SRP151137_CpG_meth_data <- read.csv("data/SRP151137/SRP151137_CpGgene_meth.csv",row.names = 1)
SRP151137_GpC_meth_data <- read.csv("data/SRP151137/SRP151137_GpCgene_meth.csv",row.names = 1)


SRP151137_Epis_sample <- read.csv("data/SRP151137/sample/SRP151137_CpGGpC_RNA_sample_data.csv")
# SRP151137_Epis_sample$Run_methlevel <- paste0(SRP151137_Epis_sample$Run, ".methlevel")

SRP151137_CpG_methlevel <- Read_file_meth_colname(SRP151137_CpG_meth_level_data,"methlevel")
SRP151137_GpC_methlevel <- Read_file_meth_colname(SRP151137_GpC_meth_level_data,"methlevel")

SRP151137_CpG_meth <- Read_file_meth_colname(SRP151137_CpG_meth_data,"meth")
SRP151137_CpG_UNmeth <- Read_file_meth_colname(SRP151137_CpG_meth_data,"UNmeth")
SRP151137_GpC_meth <- Read_file_meth_colname(SRP151137_GpC_meth_data,"meth")
SRP151137_GpC_UNmeth <- Read_file_meth_colname(SRP151137_GpC_meth_data,"UNmeth")

# CpG_GpC_sample <- intersect(colnames(SRP151137_CpG_methlevel),colnames(SRP151137_GpC_methlevel))

SRP151137_CpG_methlevel_filter <- SRP151137_CpG_methlevel[,factor(SRP151137_Epis_sample$Run_methlevel)]
SRP151137_GpC_methlevel_filter <- SRP151137_GpC_methlevel[,factor(SRP151137_Epis_sample$Run_methlevel)]

SRP151137_CpG_meth_filter <- SRP151137_CpG_meth[,factor(SRP151137_Epis_sample$Run_meth)]
SRP151137_CpG_UNmeth_filter <- SRP151137_CpG_UNmeth[,factor(SRP151137_Epis_sample$Run_UNmeth)]
SRP151137_GpC_meth_filter <- SRP151137_GpC_meth[,factor(SRP151137_Epis_sample$Run_meth)]
SRP151137_GpC_UNmeth_filter <- SRP151137_GpC_UNmeth[,factor(SRP151137_Epis_sample$Run_UNmeth)]

# any(colnames(SRP151137_CpG_methlevel_filter)==colnames(SRP151137_GpC_methlevel_filter))
# SRP151137_Epis_sample_filter <- SRP151137_Epis_sample %>%
#   filter(Run_methlevel %in% CpG_GpC_sample)

colnames(SRP151137_CpG_methlevel_filter) <- SRP151137_Epis_sample$Title
colnames(SRP151137_GpC_methlevel_filter) <- SRP151137_Epis_sample$Title

colnames(SRP151137_CpG_meth_filter) <- SRP151137_Epis_sample$Title
colnames(SRP151137_CpG_UNmeth_filter) <- SRP151137_Epis_sample$Title
colnames(SRP151137_GpC_meth_filter) <- SRP151137_Epis_sample$Title
colnames(SRP151137_GpC_UNmeth_filter) <- SRP151137_Epis_sample$Title


#chr data----
GRCm38_data <- read.csv("data/list/GRCm38_Genes.csv")
GRCm38_data_genes <- read.csv("data/list/GRCm38_Genes_genename.csv")
GRCm38_data_paste <- Chr_region_process(GRCm38_data,"paste")
GRCm38_data_paste$geneid <- GRCm38_data_genes$gene_id
GRCm38_data_paste$genename <- GRCm38_data_genes$gene_name

remove_after_last_dot <- function(x) {
  sub("\\.[^.]*$", "", x)  # 正则替换最后一次出现的点及后续内容
}

GRCm38_data_paste$geneid <- remove_after_last_dot(GRCm38_data_paste$geneid)


#NA process----
SRP151137_CpG_methlevel_na_count <- NA_count_order(SRP151137_CpG_methlevel_filter)
SRP151137_GpC_methlevel_na_count <- NA_count_order(SRP151137_GpC_methlevel_filter)

SRP151137_CpG_methlevel_na_count_nona <- SRP151137_CpG_methlevel_na_count %>%
  filter(NA_percentage<=0.65) #0.45

SRP151137_CpG_methlevel_choose <- SRP151137_CpG_methlevel_filter[rownames(SRP151137_CpG_methlevel_na_count_nona),]

SRP151137_GpC_methlevel_na_count_nona <- SRP151137_GpC_methlevel_na_count %>%
  filter(NA_percentage<=0.9) #0.85

SRP151137_GpC_methlevel_choose <- SRP151137_GpC_methlevel_filter[rownames(SRP151137_GpC_methlevel_na_count_nona),]

# save.image("SRP151137_plot/Gene_All_Methlevel_means.RData")
load("SRP151137_plot/Gene_All_MOFA.RData")

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

# SRP151137_CpG_mofa <- run_umap(SRP151137_CpG_mofa)
SRP151137_CpG_mofa <- run_umap(SRP151137_CpG_mofa,
                               min_dist = 2,spread=3)

plot_dimred(SRP151137_CpG_mofa, method = "UMAP",color_by = "group",dot_size = 2,label = F)+
  theme(legend.title=element_blank())

SRP151137_CpG_umap_data <- SRP151137_CpG_mofa@dim_red$UMAP
SRP151137_CpG_umap_data <- cbind(SRP151137_CpG_umap_data,SRP151137_CpG_mofa@samples_metadata$group)
colnames(SRP151137_CpG_umap_data)[4] <- "group"
rownames(SRP151137_CpG_umap_data) <- 1:nrow(SRP151137_CpG_umap_data)

# write.csv(SRP151137_CpG_umap_data,"SRP151137_plot/gene_data/SRP151137_CpGgene_umap_data.csv")

ggplot(SRP151137_CpG_umap_data)+
  geom_point(SRP151137_CpG_umap_data,mapping=aes(UMAP1,UMAP2,color=group),size = 2.5,alpha = 0.7)+
  theme_classic()+
  labs(x = "UMAP1",y = "UMAP2",title = "Gene CpG UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
        axis.title.x = element_text(face = "bold",size = 15),
        axis.title.y = element_text(face = "bold",size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))+
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
                               min_dist = 2,spread=3)

plot_dimred(SRP151137_GpC_mofa, method = "UMAP",color_by = "group",dot_size = 2,label = F)+
  theme(legend.title=element_blank())

SRP151137_GpC_umap_data <- SRP151137_GpC_mofa@dim_red$UMAP
SRP151137_GpC_umap_data <- cbind(SRP151137_GpC_umap_data,SRP151137_GpC_mofa@samples_metadata$group)
colnames(SRP151137_GpC_umap_data)[4] <- "group"
rownames(SRP151137_GpC_umap_data) <- 1:nrow(SRP151137_GpC_umap_data)

# write.csv(SRP151137_GpC_umap_data,"SRP151137_plot/gene_data/SRP151137_GpCgene_umap_data.csv")


ggplot(SRP151137_GpC_umap_data)+
  geom_point(SRP151137_GpC_umap_data,mapping=aes(UMAP1,UMAP2,color=group),size = 2.5,alpha = 0.7)+
  theme_classic()+
  labs(x = "UMAP1",y = "UMAP2",title = "Gene GpC UMAP")+
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
        axis.title.x = element_text(face = "bold",size = 15),
        axis.title.y = element_text(face = "bold",size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))+
  scale_color_manual(breaks = c("E4.5","E5.5", "E6.5"),
                     values = c("#761C78","#5B6C2F","#4682B4"))


# load("SRP151137_plot/Gene_fisher_DEG.RData")

#CpG DEG----
SRP151137_CpG_DEG <- Methlevel_group_variance_analysis(SRP151137_CpG_methlevel_choose_mean,SRP151137_Epis_sample$group1)
# write.csv(SRP151137_CpG_DEG,"SRP151137_plot/gene_data/SRP151137_CpG_DEG.csv")
SRP151137_CpG_DEG_gene <- merge(
  SRP151137_CpG_DEG,
  GRCm38_data_paste,
  by = "chr",
  all.x = TRUE
)


#CpG E4.5
SRP151137_CpG_DEG_E4.5 <- Fisher_DEG_process(
  SRP151137_CpG_meth_filter,
  SRP151137_CpG_UNmeth_filter,
  SRP151137_Epis_sample,
  "E4.5")
CpGE4.5_methlevel_diff_data <- Methlevel_diff(SRP151137_CpG_methlevel_filter,SRP151137_Epis_sample,"E4.5")
SRP151137_CpG_DEG_E4.5 <- cbind(SRP151137_CpG_DEG_E4.5,CpGE4.5_methlevel_diff_data)

SRP151137_CpG_DEG_E4.5$chr <- rownames(SRP151137_CpG_DEG_E4.5)
SRP151137_CpG_DEG_E4.5$geneid <- GRCm38_data_paste$geneid
SRP151137_CpG_DEG_E4.5$genename <- GRCm38_data_paste$genename


#CpG E5.5
SRP151137_CpG_DEG_E5.5 <- Fisher_DEG_process(
  SRP151137_CpG_meth_filter,
  SRP151137_CpG_UNmeth_filter,
  SRP151137_Epis_sample,
  "E5.5")
CpGE5.5_methlevel_diff_data <- Methlevel_diff(SRP151137_CpG_methlevel_filter,SRP151137_Epis_sample,"E5.5")
SRP151137_CpG_DEG_E5.5 <- cbind(SRP151137_CpG_DEG_E5.5,CpGE5.5_methlevel_diff_data)

SRP151137_CpG_DEG_E5.5$chr <- rownames(SRP151137_CpG_DEG_E5.5)
SRP151137_CpG_DEG_E5.5$geneid <- GRCm38_data_paste$geneid
SRP151137_CpG_DEG_E5.5$genename <- GRCm38_data_paste$genename


#CpG E6.5
SRP151137_CpG_DEG_E6.5 <- Fisher_DEG_process(
  SRP151137_CpG_meth_filter,
  SRP151137_CpG_UNmeth_filter,
  SRP151137_Epis_sample,
  "E6.5")
CpGE6.5_methlevel_diff_data <- Methlevel_diff(SRP151137_CpG_methlevel_filter,SRP151137_Epis_sample,"E6.5")
SRP151137_CpG_DEG_E6.5 <- cbind(SRP151137_CpG_DEG_E6.5,CpGE6.5_methlevel_diff_data)

SRP151137_CpG_DEG_E6.5$chr <- rownames(SRP151137_CpG_DEG_E6.5)
SRP151137_CpG_DEG_E6.5$geneid <- GRCm38_data_paste$geneid
SRP151137_CpG_DEG_E6.5$genename <- GRCm38_data_paste$genename


SRP151137_CpG_DEG_P_E4.5 <- SRP151137_CpG_DEG_E4.5 %>%
  filter(P.value < 5e-10)
SRP151137_CpG_DEG_P_E5.5 <- SRP151137_CpG_DEG_E5.5 %>%
  filter(P.value < 0.05)
SRP151137_CpG_DEG_P_E6.5 <- SRP151137_CpG_DEG_E6.5 %>%
  filter(P.value < 0.05)



SRP151137_CpG_DEG_gene_P <- SRP151137_CpG_DEG_gene %>%
  filter(P.value < 0.05)
# write.csv(SRP151137_CpG_DEG_P,"SRP151137_plot/gene_data/SRP151137_CpG_DEG_P0.05.csv")
table(SRP151137_CpG_DEG_gene_P$group)

SRP151137_CpG_DEG_gene_P %>%
  group_by(group) %>%
  top_n(n = 500, wt = logFC)


volcano_data <- SRP151137_CpG_DEG_gene_P %>%
  mutate(
    chr = SRP151137_CpG_DEG_gene_P$chr,  # 确保有基因名列
    log_pval = -log10(P.value),  # 转换校正p值
    direction = case_when(         # 标记上下调
      logFC > 0.5 & P.value < 5e-10 ~ "Up",
      logFC < -0.5 & P.value < 5e-10 ~ "Down",
      TRUE ~ "Not sig"
    )
  ) %>%
  # 过滤无限值（避免log(0)错误）
  filter(!is.infinite(log_pval))

top_genes_up <- volcano_data %>%
  filter(direction == "Up") %>%
  arrange(P.value) %>%
  head(5)

top_genes_down <- volcano_data %>%
  filter(direction == "Down") %>%
  arrange(P.value) %>%
  head(5)


ggplot(volcano_data, aes(x = logFC, y = log_pval, colour = direction)) +
  geom_point(alpha=0.7, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2", "#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(5e-10),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(Fold Change)",y="-log10 (P-value)")+
  theme_bw()+
  ggtitle("CpG Gene Methlevel DEG Volcano")+
  # 图例
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5,size = 20,face = "bold"),
        axis.title.x = element_text(face = "bold",size = 15),
        axis.title.y = element_text(face = "bold",size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  xlim(-2, 2) + ylim(0, 50)+
  geom_label_repel(data = top_genes_up, aes(label = genename),
                   size = 4,                           # 设置标签大小
                   box.padding = unit(0.8, "lines"),   # 设置标签内边距
                   point.padding = unit(0.8, "lines"), # 设置标签与点的距离
                   segment.color = "black",            # 设置标签边界线颜色
                   show.legend = FALSE,                # 不显示图例
                   max.overlaps = 10000)+             # 设置标签重叠的最大次数
  geom_label_repel(data = top_genes_down, aes(label = genename),
                   size = 4,                           # 设置标签大小
                   box.padding = unit(0.8, "lines"),   # 设置标签内边距
                   point.padding = unit(0.8, "lines"), # 设置标签与点的距离
                   segment.color = "black",            # 设置标签边界线颜色
                   show.legend = FALSE,                # 不显示图例
                   max.overlaps = 10000)               # 设置标签重叠的最大次数






#CpG DEG enrichment----
# SRP151137_CpG_DEG_gene <- bitr(SRP151137_CpG_DEG_choose_chr_gene$gene,
#                                fromType = "E",
#                                toType = "ENTREZID",
#                                OrgDb = org.Mm.eg.db)

ego <- enrichGO(gene = SRP151137_CpG_DEG_gene_P$geneid,
                keyType = "ENSEMBL",
                OrgDb = org.Mm.eg.db,
                ont = "ALL",        # "BP","MF","CC"或"ALL"
                pAdjustMethod = "BH",         # 校正方法：BH, bonferroni等
                pvalueCutoff = 0.05,         # 显著阈值
                qvalueCutoff = 0.05,
                readable = TRUE)         # 转换ID为基因名

dotplot(ego,
        showCategory=10,         # 显示top15条目
        # font.size=10,
        title="CpG Gene GO Enrichment") +
  scale_color_gradient(low="#546de5", high="#ff4757")+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))


#GpC DEG----
SRP151137_GpC_DEG <- Methlevel_group_variance_analysis(SRP151137_GpC_methlevel_choose_mean,SRP151137_Epis_sample$group1)
# write.csv(SRP151137_GpC_DEG,"SRP151137_plot/gene_data/SRP151137_GpC_DEG.csv")
SRP151137_GpC_DEG_gene <- merge(
  SRP151137_GpC_DEG,
  GRCm38_data_paste,
  by = "chr",
  all.x = TRUE
)


#GpC E4.5
SRP151137_GpC_DEG_E4.5 <- Fisher_DEG_process(
  SRP151137_GpC_meth_filter,
  SRP151137_GpC_UNmeth_filter,
  SRP151137_Epis_sample,
  "E4.5")
GpCE4.5_methlevel_diff_data <- Methlevel_diff(SRP151137_GpC_methlevel_filter,SRP151137_Epis_sample,"E4.5")
SRP151137_GpC_DEG_E4.5 <- cbind(SRP151137_GpC_DEG_E4.5,GpCE4.5_methlevel_diff_data)

SRP151137_GpC_DEG_E4.5$chr <- rownames(SRP151137_GpC_DEG_E4.5)
SRP151137_GpC_DEG_E4.5$geneid <- GRCm38_data_paste$geneid
SRP151137_GpC_DEG_E4.5$genename <- GRCm38_data_paste$genename


#GpC E5.5
SRP151137_GpC_DEG_E5.5 <- Fisher_DEG_process(
  SRP151137_GpC_meth_filter,
  SRP151137_GpC_UNmeth_filter,
  SRP151137_Epis_sample,
  "E5.5")
GpCE5.5_methlevel_diff_data <- Methlevel_diff(SRP151137_GpC_methlevel_filter,SRP151137_Epis_sample,"E5.5")
SRP151137_GpC_DEG_E5.5 <- cbind(SRP151137_GpC_DEG_E5.5,GpCE5.5_methlevel_diff_data)

SRP151137_GpC_DEG_E5.5$chr <- rownames(SRP151137_GpC_DEG_E5.5)
SRP151137_GpC_DEG_E5.5$geneid <- GRCm38_data_paste$geneid
SRP151137_GpC_DEG_E5.5$genename <- GRCm38_data_paste$genename


#GpC E6.5
SRP151137_GpC_DEG_E6.5 <- Fisher_DEG_process(
  SRP151137_GpC_meth_filter,
  SRP151137_GpC_UNmeth_filter,
  SRP151137_Epis_sample,
  "E6.5")
GpCE6.5_methlevel_diff_data <- Methlevel_diff(SRP151137_GpC_methlevel_filter,SRP151137_Epis_sample,"E6.5")
SRP151137_GpC_DEG_E6.5 <- cbind(SRP151137_GpC_DEG_E6.5,GpCE6.5_methlevel_diff_data)

SRP151137_GpC_DEG_E6.5$chr <- rownames(SRP151137_GpC_DEG_E6.5)
SRP151137_GpC_DEG_E6.5$geneid <- GRCm38_data_paste$geneid
SRP151137_GpC_DEG_E6.5$genename <- GRCm38_data_paste$genename



SRP151137_GpC_DEG_P_E4.5 <- SRP151137_GpC_DEG_E4.5 %>%
  filter(P.value < 0.5)
SRP151137_GpC_DEG_P_E5.5 <- SRP151137_GpC_DEG_E5.5 %>%
  filter(P.value < 0.5)
SRP151137_GpC_DEG_P_E6.5 <- SRP151137_GpC_DEG_E6.5 %>%
  filter(P.value < 0.5)



SRP151137_GpC_DEG_choose <- SRP151137_GpC_DEG_P_E4.5 %>%
  group_by(group) %>%
  top_n(n = 5, wt = logFC)

SRP151137_CpG_DEG_gene_P %>%
  group_by(group) %>%
  top_n(n = 500, wt = logFC)

volcano_data <- SRP151137_GpC_DEG_gene %>%
  mutate(
    chr = SRP151137_GpC_DEG_gene$chr,  # 确保有基因名列
    log_pval = -log10(P.value),  # 转换校正p值
    direction = case_when(         # 标记上下调
      logFC > 0.5 & P.value < 0.5 ~ "Up",
      logFC < -0.5 & P.value < 0.5 ~ "Down",
      TRUE ~ "Not sig"
    )
  ) %>%
  # 过滤无限值（避免log(0)错误）
  filter(!is.infinite(log_pval))

ggplot(volcano_data, aes(x = logFC, y = log_pval, colour = direction)) +
  geom_point(alpha=0.5, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2", "#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(Fold Change)",y="-log10 (P-value)")+
  theme_bw()+
  ggtitle("GpC Gene DEG Volcano")+
  # 图例
  theme(legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5,size = 18,face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.title.y = element_text(face = "bold",size = 14))+
  xlim(-2, 2) + ylim(0, 50)


#GpC DEG enrichment----
# SRP151137_GpC_DEG_gene <- bitr(SRP151137_GpC_DEG_choose_chr_gene$gene,
#                                fromType = "SYMBOL",
#                                toType = "ENTREZID",
#                                OrgDb = org.Mm.eg.db)

ego <- enrichGO(gene = SRP151137_GpC_DEG_P_E4.5$geneid,
                keyType = "ENSEMBL",
                OrgDb = org.Mm.eg.db,
                ont = "ALL",        # "BP","MF","CC"或"ALL"
                pAdjustMethod = "BH",         # 校正方法：BH, bonferroni等
                pvalueCutoff = 0.05,         # 显著阈值
                qvalueCutoff = 0.05,
                readable = TRUE)         # 转换ID为基因名

dotplot(ego,
        showCategory=10,         # 显示top15条目
        # font.size=10,
        title="CpG Gene GO Enrichment") +
  scale_color_gradient(low="#546de5", high="#ff4757")+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))



#CpG GpC DEG diff methlevel----
SRP151137_E4.5_diff_data <- DEG_diff_plot(SRP151137_CpG_DEG_E4.5,SRP151137_GpC_DEG_E4.5)
SRP151137_E4.5_diff_data_choose <- SRP151137_E4.5_diff_data %>%
  mutate(avg_p.value = (GpC_p.value + CpG_p.value)/2) %>%
  arrange(avg_p.value)
# SRP151137_E4.5_diff_data_choose <- SRP151137_E4.5_diff_data %>%
#   mutate(avg_logFC = (GpC_logFC + CpG_logFC)/2) %>%
#   arrange(avg_logFC)
SRP151137_E4.5_diff_data_choose$top10  <- FALSE
SRP151137_E4.5_diff_data_choose$top10[1:10] <- TRUE  # 前10行标记为TRUE
SRP151137_E4.5_diff_data_choose_p <- SRP151137_E4.5_diff_data_choose %>%
  # filter(avg_value <= 0.05 & logFC > 0)
  filter(avg_p.value <= 0.5)


SRP151137_E5.5_diff_data <- DEG_diff_plot(SRP151137_CpG_DEG_E5.5,SRP151137_GpC_DEG_E5.5)
SRP151137_E5.5_diff_data_choose <- SRP151137_E5.5_diff_data %>%
  mutate(avg_p.value = (GpC_p.value + CpG_p.value)/2) %>%
  arrange(avg_p.value)
# SRP151137_E5.5_diff_data_choose <- SRP151137_E5.5_diff_data %>%
#   mutate(avg_logFC = (GpC_logFC + CpG_logFC)/2) %>%
#   arrange(avg_logFC)
SRP151137_E5.5_diff_data_choose$top10  <- FALSE
SRP151137_E5.5_diff_data_choose$top10[1:10] <- TRUE  # 前10行标记为TRUE
SRP151137_E5.5_diff_data_choose_p <- SRP151137_E5.5_diff_data_choose %>%
  # filter(avg_value <= 0.05 & logFC > 0)
  filter(avg_p.value <= 0.5)


SRP151137_E6.5_diff_data <- DEG_diff_plot(SRP151137_CpG_DEG_E6.5,SRP151137_GpC_DEG_E6.5)
SRP151137_E6.5_diff_data_choose <- SRP151137_E6.5_diff_data %>%
  mutate(avg_p.value = (GpC_p.value + CpG_p.value)/2) %>%
  arrange(avg_p.value)
# SRP151137_E6.5_diff_data_choose <- SRP151137_E6.5_diff_data %>%
#   mutate(avg_logFC = (GpC_logFC + CpG_logFC)/2) %>%
#   arrange(avg_logFC)
SRP151137_E6.5_diff_data_choose$top10  <- FALSE
SRP151137_E6.5_diff_data_choose$top10[1:10] <- TRUE  # 前10行标记为TRUE
SRP151137_E6.5_diff_data_choose_p <- SRP151137_E6.5_diff_data_choose %>%
  # filter(avg_value <= 0.05 & logFC > 0)
  filter(avg_p.value <= 0.5)


# save.image("SRP151137_plot/Gene_fisher_DEG.RData")



ggplot_data <- SRP151137_E4.5_diff_data_choose_p


ggplot(data = ggplot_data,mapping = aes(x = CpG_rate_diff, y = GpC_rate_diff))+
  geom_point(size = 2,alpha = 0.9,color = "grey75")+
  geom_point(data = subset(ggplot_data, top10),
             color = "black", size = 2)+
  theme_classic()+
  geom_vline(xintercept = 0, color = "#FF7F00", linetype = "dashed") +  # 垂直参考线（x=50）
  geom_hline(yintercept = 0, color = "#FF7F00", linetype = "dashed")+   # 水平参考线（y=50）
  labs(x = "Differential methylation (%)",y = "Differential accessibility (%)",
       title = "E4.5 Gene Differential")+
  theme(legend.title=element_blank(),
        # plot.title = element_text(hjust = 0.5,size = 18,face = "bold"),
        plot.title = element_textbox(
          # family = "serif",  # 字体
          face = "bold",
          size = 20,
          linewidth = 0,    # 边框粗细
          linetype = 0,       # 实线边框
          lineheight = 1.5,
          color = "black",    # 文字颜色
          fill = "#deb4bd",   # 背景填充色(浅蓝)
          box.color = "#deb4bd",  # 边框颜色(深蓝)
          hjust = 0.5,       # 水平居中
          # vjust = 0.5,      # 垂直居中
          padding = margin(8, 15, 8, 15),  # 内边距（上、右、下、左）
          margin = margin(b = 10),   # 标题下方外边距
        ),
        axis.text.x = element_text(face = "bold",size = 13),
        axis.text.y = element_text(face = "bold",size = 13),
        axis.title.x = element_text(face = "bold",size = 15),
        axis.title.y = element_text(face = "bold",size = 15))+
  coord_cartesian(xlim = c(-100, 100), ylim = c(-100, 100))+
  geom_label_repel(
    data = subset(ggplot_data, top10),
    aes(label = genename),
    color = "#ca4b99",
    fill = alpha("white", 0.6),  # 半透明白底
    box.padding = 0.5,           # 标签与点的间距
    segment.color = "grey50",     # 连接线颜色
    max.overlaps = 20             # 允许更多重叠标签
  )



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
       title = "Gene CpG ColMean")+
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 12,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 14,face = "bold"),
        plot.title = element_text(size = 18,hjust = 0.5,face = "bold"))+
  # scale_fill_manual(values = group_colors)+
  guides(fill="none")+
  scale_fill_manual(breaks = c("E4.5","E5.5", "E6.5"),
                     values = c("#761C78","#5B6C2F","#4682B4"))



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
       title = "Gene GpC ColMean")+
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1, size = 12,face = "bold"),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title.y = element_text(size = 14,face = "bold"),
        plot.title = element_text(size = 18,hjust = 0.5,face = "bold"))+
  # scale_fill_manual(values = group_colors)+
  guides(fill="none")+
  scale_fill_manual(breaks = c("E4.5","E5.5", "E6.5"),
                    values = c("#761C78","#5B6C2F","#4682B4"))




#Venn genes----
load("SRP151137_plot/RNA_markers.RData")
# load("SRP151137_plot/Gene_All_Meth.RData")

SRP151137_pbmc_markers_E4.5 <- SRP151137_pbmc_markers_result %>%
  filter(cluster == "E4.5")
SRP151137_pbmc_markers_E5.5 <- SRP151137_pbmc_markers_result %>%
  filter(cluster == "E5.5")
SRP151137_pbmc_markers_E6.5 <- SRP151137_pbmc_markers_result %>%
  filter(cluster == "E6.5")

SRP151137_pbmc_markers_E4.5_P <- SRP151137_pbmc_markers_E4.5 %>%
  filter(p_val < 5e-4)
SRP151137_pbmc_markers_E5.5_P <- SRP151137_pbmc_markers_E5.5 %>%
  filter(p_val < 5e-5)
SRP151137_pbmc_markers_E6.5_P <- SRP151137_pbmc_markers_E6.5 %>%
  filter(p_val < 5e-6)


SRP151137_CpG_DEG_E4.5_venn <- SRP151137_CpG_DEG_E4.5 %>%
  filter(P.value <= 5e-40)
SRP151137_CpG_DEG_E5.5_venn <- SRP151137_CpG_DEG_E5.5 %>%
  filter(P.value <= 5e-4)
SRP151137_CpG_DEG_E6.5_venn <- SRP151137_CpG_DEG_E6.5 %>%
  filter(P.value <= 5e-20)


SRP151137_GpC_DEG_E4.5_venn <- SRP151137_GpC_DEG_E4.5 %>%
  filter(P.value <= 0.5)
SRP151137_GpC_DEG_E5.5_venn <- SRP151137_GpC_DEG_E5.5 %>%
  filter(P.value <= 0.5)
SRP151137_GpC_DEG_E6.5_venn <- SRP151137_GpC_DEG_E6.5 %>%
  filter(P.value <= 0.5)


#ggplot venn
genes <- list(
  RNA_Gene = SRP151137_pbmc_markers_E4.5$Geneid2,
  CpG_Gene = SRP151137_CpG_DEG_E4.5_venn$geneid,
  GpC_Gene = SRP151137_GpC_DEG_E4.5_venn$geneid
)

ggVennDiagram(genes,
              label_alpha = 0,  # 标签透明背景
              edge_size = 0.5) +
  ggplot2::scale_fill_gradient(low = "white", high = "firebrick") +
  ggplot2::labs(title = "E6.5 Gene Overlap")+
  theme(
    plot.title = element_text(color="black",hjust=0.5,vjust=0.5,size=16,face="bold")
    # axis.title.x = element_text(size = 14,face = "bold"),
    # axis.title.y = element_text(size = 14,face = "bold")
  )



inter <- get.venn.partitions(genes)
for (i in 1:nrow(inter)) {
  inter[i,'Genes'] <- paste(inter[[i,'..values..']], collapse=", ")
}

venn_genes <- inter$Genes[1]
venn_genes <- unlist(strsplit(venn_genes, ","))
# write.csv(venn_genes,"SRP151137_plot/RNA_gene/E4.5_CpG_RNA.genes.csv")


# venn_DEG_genes <- bitr(venn_genes,
#                        fromType = "SYMBOL",
#                        toType = "ENTREZID",
#                        OrgDb = org.Mm.eg.db)

ego <- enrichGO(gene = venn_genes,
                OrgDb = org.Mm.eg.db,
                keyType = "ENSEMBL",
                ont = "ALL",        # "BP","MF","CC"或"ALL"
                pAdjustMethod = "BH",         # 校正方法：BH, bonferroni等
                pvalueCutoff = 0.05,         # 显著阈值
                qvalueCutoff = 0.05,
                readable = TRUE)         # 转换ID为基因名
ego_data <- ego@result


dotplot(ego,
        showCategory=10,         # 显示top15条目
        # font.size=10,
        title="Venn Gene GO Enrichment") +
  scale_color_gradient(low="#546de5", high="#ff4757")+
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(face = "bold",size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))




#Calinski-Harabasz----
data <- SRP151137_GpC_umap_data

coords <- data[, c("UMAP1", "UMAP2")]
clustering <- as.factor(data$group)

# 计算Calinski-Harabasz指数
ch_index <- calinhara(coords, as.numeric(clustering))
ch_index


# CH_SC <- read.csv("SRP151137_plot/encode_gene/CH_SC.csv")

ggplot(CH_SC, aes(x = type, y = CH.value)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.5, fill = "#4E79A7")+
  # geom_col(position = "dodge")+
  geom_text(aes(label = round(CH.value, 2)),
            position = position_dodge(0.8),
            vjust = -0.5,
            size = 4) +
  # scale_fill_manual(values = c("CpG.ENCODE" = "#4E79A7", "GpC.ENCODE" = "#F28E2B",
  #                              "CpG.Gene" = "#5B6C2F", "GpC.Gene" = "#761C78")) +
  theme_minimal()+
  labs(title = "Calinski-Harabasz",
       x = "",
       y = "Value")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    # legend.position = "top",
    legend.title = element_blank(),
    # panel.grid.major.x = element_blank(),
    # axis.title.x = element_text(face = "bold",size = 12),
    axis.title.y = element_text(face = "bold",size = 12),
    axis.text.x = element_text(face = "bold",size = 10),
    axis.text.y = element_text(face = "bold",size = 10)
  )



#Silhouette Coefficient----
data <- SRP151137_GpC_umap_data

# 提取UMAP坐标和分组标签
umap_data <- data[, c("UMAP1", "UMAP2")]
groups <- data$group

# 计算距离矩阵（欧氏距离）
dist_matrix <- dist(umap_data, method = "euclidean")

# 将分组标签转换为数字因子（轮廓系数需要）
group_factor <- as.numeric(factor(groups))

# 计算轮廓系数
sil <- silhouette(group_factor, dist_matrix)

# 计算平均轮廓系数
avg_sil_width <- mean(sil[, "sil_width"])
cat("平均轮廓系数:", avg_sil_width, "\n")




# CH_SC <- read.csv("SRP151137_plot/CH_SC.csv")

ggplot(CH_SC, aes(x = type, y = SC.value)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.5, fill = "#4E79A7")+
  # geom_col(position = "dodge")+
  geom_text(aes(label = round(SC.value, 2)),
            position = position_dodge(0.8),
            vjust = -0.5,
            size = 4) +
  # scale_fill_manual(values = c("CpG.ENCODE" = "#4E79A7", "GpC.ENCODE" = "#F28E2B",
  #                              "CpG.Gene" = "#5B6C2F", "GpC.Gene" = "#761C78")) +
  theme_minimal()+
  labs(title = "Silhouette Coefficient",
       x = "",
       y = "Value")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    # legend.position = "top",
    legend.title = element_blank(),
    # panel.grid.major.x = element_blank(),
    # axis.title.x = element_text(face = "bold",size = 12),
    axis.title.y = element_text(face = "bold",size = 12),
    axis.text.x = element_text(face = "bold",size = 10),
    axis.text.y = element_text(face = "bold",size = 10)
  )
