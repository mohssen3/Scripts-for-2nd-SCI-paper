library(tidyverse)
library(ComplexHeatmap)

setwd("~/anvio_enrichment_analysis_VIPs") #set to this directory
enriched_modules<-
  read.csv("./VIPs_enriched_modules.txt",sep = "\t") %>%
  filter(adjusted_q_value<=0.05 | (p_Decreasing_with_injury>=0 & p_Increasing_with_injury==0) | (p_Increasing_with_injury>=0 & p_Decreasing_with_injury==0))

increasing_or_decreasing_VIPs<-
  read.csv("./VIP_groups.txt",sep="\t") %>%
  rename(mag=sample) %>%
  mutate_all(~gsub("CAG.","CAG-",.)) %>%
  mutate_all(~gsub("X1XD8.7","1XD8-7",.))

enriched_modules_accession<- enriched_modules %>%
  select(accession)

VIP_completeness_matrix_df<- 
  read.csv("./VIPs-completeness-MATRIX.txt",
           sep="\t")

enriched_modules_matrix<-
  inner_join(VIP_completeness_matrix_df,enriched_modules_accession,by=c("module"="accession")) %>%
  pivot_longer(-module,names_to = "mag",values_to = "completeness") %>%
  mutate_all(~gsub("CAG.","CAG-",.)) %>%
  mutate_all(~gsub("X1XD8.7","1XD8-7",.)) %>%
  pivot_wider(id_cols = module,names_from = mag,values_from = completeness)


VIP_completeness_df<-
  enriched_modules_matrix %>%
  pivot_longer(-module,names_to = "mag",
               values_to = "completeness") %>%
  inner_join(increasing_or_decreasing_VIPs)
VIP_completeness_df_f<-
  inner_join(enriched_modules,VIP_completeness_df,by=c("accession"="module")) %>%
  select(mag,group,KEGG_MODULE,completeness) %>%
  pivot_wider(id_cols = c(mag,group),
              names_from = "KEGG_MODULE",
              values_from = "completeness")

VIP_completeness_df_f


VIP_presence_matrix_df<- read.csv("./VIPs-presence-MATRIX.txt",sep="\t") %>%
  pivot_longer(-module,names_to = "mag",values_to = "completeness") %>%
  mutate_all(~gsub("CAG.","CAG-",.)) %>%
  mutate_all(~gsub("X1XD8.7","1XD8-7",.)) %>%
  pivot_wider(id_cols = module,names_from = mag,values_from = completeness)
enriched_modules_presence_matrix<-
  inner_join(VIP_presence_matrix_df,enriched_modules_accession,by=c("module"="accession"))

VIP_presence_df<-
  enriched_modules_presence_matrix %>%
  pivot_longer(-module,names_to = "mag",
               values_to = "completeness") %>%
  inner_join(increasing_or_decreasing_VIPs)
VIP_presence_df_f<-
  inner_join(enriched_modules,VIP_presence_df,by=c("accession"="module")) %>%
  select(mag,group,KEGG_MODULE,completeness) %>%
  pivot_wider(id_cols = c(mag,group),
              names_from = "KEGG_MODULE",
              values_from = "completeness")
VIP_presence_df_f

rownames(VIP_presence_df_f) <- VIP_presence_df_f$mag

VIP_presence_num<-as.matrix(VIP_presence_df_f[, -1:-2])
VIP_presence_num<-matrix(as.numeric(VIP_presence_num),nrow = 19,ncol = 52)
colnames(VIP_presence_num) <-colnames(VIP_presence_df_f[, -1:-2])
rownames(VIP_presence_num) <- VIP_presence_df_f$mag

# colnames(enriched_modules_matrix) <-c("ko_id","7dpi","21dpi","35dpi","63dpi","180dpi")
# VIP_presence_matrix <- matrix(VIP_presence_num,nrow = 19,ncol=30)
# matrix(a,nrow = dim(df2)[1],ncol = 5)
#column_ha = HeatmapAnnotation(foo1 = runif(10), bar1 = anno_barplot(runif(10)))
ha = rowAnnotation(Group=VIP_presence_df_f$group,
                   col = list(Group = c("Increasing_with_injury" ="#fb8500", "Decreasing_with_injury" = "#219ebc")))#,
#                   GeneDescription=enriched_modules$KEGG_MODULE)
#                   module=VIP_presence_df$group)
library(cluster)
heatmap_obj<-Heatmap(VIP_presence_num,na_col = "black",name = "mat",cluster_rows = TRUE,cluster_columns = FALSE,
                     column_names_rot = 75,show_row_names = TRUE,right_annotation = ha,
                     clustering_method_rows = "ward.D",
                     #                     cluster_row_slices = cluster_within_group(VIP_presence_num, VIP_presence_df_f$group),
                     #                     row_dend_reorder = TRUE,
                     show_column_dend = FALSE,
                     col = c("#F4F4F4","#004346"),
                     width = unit(15, "cm"), 
                     height = unit(12, "cm"),
                     column_names_gp = gpar(fontsize = 9),
                     #                     column_names_side = "top",
                     heatmap_legend_param = list(
                       at = c(0,1),
                       labels = c("<75%", ">=75%"),
                       title = "Completeness",
                       legend_height = unit(4, "cm"),
                       title_position = "topcenter"
                     ),
                     row_names_max_width = max_text_width(rownames(VIP_presence_num))
)
heatmap_obj
pdf("enriched_pathways_VIP_anvio_presence_matrix_1.pdf", width = 13, height = 14)# Specify the dimensions
draw(heatmap_obj)
dev.off()
# 
# png(file="enriched_pathways_VIP_anvio_presence_matrix.png",
#     width=13, height=14,units = "in",res = 300)
# draw(heatmap_obj)
# dev.off()