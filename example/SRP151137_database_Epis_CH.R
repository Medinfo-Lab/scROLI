library(fpc)       # 提供calinhara函数计算CH指数
library(dplyr)


CpG_encode_data <- read.csv("SRP151137_plot/encode_data/SRP151137_CpGencode_umap_data.csv")
GpC_encode_data <- read.csv("SRP151137_plot/encode_data/SRP151137_GpCencode_umap_data.csv")

CpG_enhancer_data <- read.csv("SRP151137_plot/enhancer_data/SRP151137_CpGenhancer_umap_data.csv")

CpG_gene_data <- read.csv("SRP151137_plot/gene_data/SRP151137_CpGgene_umap_data.csv")
GpC_gene_data <- read.csv("SRP151137_plot/gene_data/SRP151137_GpCgene_umap_data.csv")

CpG_genetss2k_data <- read.csv("SRP151137_plot/genetss2k_data/SRP151137_CpGgenetss2k_umap_data.csv")



data <- CpG_genetss2k_data


coords <- data[, c("UMAP1", "UMAP2")]
clustering <- as.factor(data$group)
# max(as.numeric(clustering))

# 计算Calinski-Harabasz指数
ch_index <- calinhara(coords, as.numeric(clustering))
ch_index


CH_SC <- read.csv("SRP151137_plot/CH_SC2.csv")

ggplot(CH_SC, aes(x = type, y = CH.value)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, fill = "#4E79A7")+
  # geom_col(position = "dodge")+
  geom_text(aes(label = round(CH.value, 2)),
            position = position_dodge(0.8),
            vjust = -0.5,
            size = 3) +
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
