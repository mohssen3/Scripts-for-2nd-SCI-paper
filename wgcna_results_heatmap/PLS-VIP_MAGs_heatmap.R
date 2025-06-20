library(tidyverse)
library(ComplexHeatmap)
###########
setwd("~/wgcna_results_heatmap/")
gtdb<- read_tsv("95_Clean_classification_dereplicated_taxonomy.tsv") %>%
  mutate(across(everything(), ~ str_replace_all(.x, "d__|p__|o__|c__|f__|g__|s__", ""))) %>%
  select(MAG,fam_species) %>%
  rename(bin=MAG)# %>%
#  mutate(bin=gsub("_bin.","_bin",bin))
# gtdb<- read.csv("95_Clean_classification_dereplicated_taxonomy.tsv",
#                 sep="\t") %>%
#   select(MAG,phylum,fam_species) %>%
#   rename(bin=MAG) %>%
#   mutate(bin=gsub("_bin.","_bin",bin))
VIPs_taxonomy<-read.csv("VIPs_taxonomy.txt",sep="\t",
                        header = FALSE) %>%
  rename(species=V2,bin=V1)

checkM<- read.csv("CheckM_galah95.txt",sep="\t") %>%
  rename(bin=Bin.Id) %>%
  select(bin,Completeness,Contamination)

VIPs<-
  read.csv("VIPs_group.txt",sep="\t") %>%
  rename(species=sample) %>%
  inner_join(VIPs_taxonomy) %>%
  inner_join(checkM) %>%
#  select(-bin) %>%
  inner_join(gtdb) %>%
  select(-bin,-fam_species)



T10<-read.csv("T10.csv") %>%
  pivot_wider(id_cols = species,names_from = Timepoint,values_from = ratio) %>%
  inner_join(VIPs,.)

heatmap<-as.matrix(T10[, -1:-4])
rownames(heatmap) <- T10$species

T10$phylum %>%
  unique
library(circlize)
set.seed(10395)
col_fun_compl = colorRamp2(c(80,100), c("#fdf0d5", "#780000"))
col_fun_cont=colorRamp2(c(min(T10$Contamination),max(T10$Contamination)), c("#fdf0d5", "#c1121f"))
col_fun <-colorRamp2(
  breaks = c(round(min(log2(heatmap), na.rm = TRUE)), 0, round(max(log2(heatmap), na.rm = TRUE))),
  colors = c("#219ebc", "#FFFFFF", "#fb8500")
)
#ha = HeatmapAnnotation(foo = 1:10, col = list(foo = col_fun))
ha = rowAnnotation(Direction=T10$group,
                   Completeness=T10$Completeness,
                   Contamination=T10$Contamination,
                   Phylum=T10$phylum,
                   col = list(Direction = c("Increasing with injury" ="#fb8500", "Decreasing with injury" = "#219ebc"),
                              Coompleteness=col_fun_compl,
                              Contamination=col_fun_cont,
                              Phylum=c("p__Firmicutes"="#ffffff",
                                       "p__Firmicutes_A"="#fca311",
                                       "p__Proteobacteria"="#3a86ff",
                                       "p__Bacteroidota"="#000000")))

fig4_t10_2<-Heatmap(log2(heatmap),na_col = "black",name = "T10",column_title = "T10/lam",
                     cluster_rows = TRUE,cluster_columns = FALSE,
                     right_annotation = ha,#left_annotation = ha_left,
                     column_names_rot = 90,show_row_names = TRUE,
                     show_column_dend = FALSE,
#                     col = c("#fb8500","#FFFFFF","#219ebc"),
                     col = c("#219ebc","#FFFFFF","#fb8500"),
                  
        
                     width = unit(7, "cm"), 
                     height = unit(9, "cm"),
                     row_names_side = "left", 
                     column_names_gp = gpar(fontsize = 10),#,font_face="Arial"),#row_names_gp = gpar(fontsize = 12,fontdace="Arial"),
                     column_names_side = "top",
                     heatmap_legend_param = list(
                       at = c(round(min(log2(heatmap))),log2(1),round(max(log2(heatmap)))),
                       labels=c(format(round(min(heatmap),2),nsmall=2),1,round(max(heatmap))),
#                       labels = c("Decreasing pathways with Injury vs Lam","Does not change with injury", "Increasing pathways with Injury vs Lam"),
                       title = "Ratio",
                       legend_height = unit(4, "cm"), legend_width=unit(4,"cm"),
                       title_position = "topcenter"
                     ),
                     row_names_max_width = max_text_width(rownames(heatmap)),
                     column_title_gp = gpar(fontsize = 14, fontface = "bold"))
fig4_t10_2
pdf("fig4_t10_2.pdf", width = 11, height = 6,res=600)# Specify the dimensions
draw(fig4_t10_2)
dev.off()

png(file="fig4_t10.png",
    width=11, height=6,units = "in",res = 300)
draw(fig4_t10)
dev.off()

#############
setwd("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/01_DRAM_VIPs")


T4<-read.csv("T4.csv") %>%
  pivot_wider(id_cols = species,names_from = Timepoint,values_from = ratio) %>%
  inner_join(VIPs,.)

heatmap<-as.matrix(T4[, -1:-4])
rownames(heatmap) <- T4$species

T4$phylum %>%
  unique
library(circlize)
set.seed(10395)
col_fun_compl = colorRamp2(c(80,100), c("#fdf0d5", "#780000"))
col_fun_cont=colorRamp2(c(min(T4$Contamination),max(T4$Contamination)), c("#fdf0d5", "#c1121f"))
col_fun <-colorRamp2(
  breaks = c(round(min(log2(heatmap), na.rm = TRUE)), 0, round(max(log2(heatmap), na.rm = TRUE))),
  colors = c("#219ebc", "#FFFFFF", "#fb8500")
)

#ha = HeatmapAnnotation(foo = 1:10, col = list(foo = col_fun))
ha = rowAnnotation(Direction=T4$group,
                   Completeness=T4$Completeness,
                   Contamination=T4$Contamination,
                   Phylum=T4$phylum,
                   col = list(Direction = c("Increasing with injury" ="#fb8500", "Decreasing with injury" = "#219ebc"),
                              Coompleteness=col_fun_compl,
                              Contamination=col_fun_cont,
                              Phylum=c("p__Firmicutes"="#ffffff",
                                       "p__Firmicutes_A"="#fca311",
                                       "p__Proteobacteria"="#3a86ff",
                                       "p__Bacteroidota"="#000000")))
fig4_t4_2<-Heatmap(log2(heatmap),na_col = "black",name = "T4",column_title = "T4/lam",
                  cluster_rows = TRUE,cluster_columns = FALSE,
                  right_annotation = ha,#left_annotation = ha_left,
                  column_names_rot = 90,show_row_names = TRUE,
                  show_column_dend = FALSE,
                 col = col_fun,

                  
                  width = unit(7, "cm"), 
                  height = unit(9, "cm"),
                  row_names_side = "left", 
                  column_names_gp = gpar(fontsize = 10),#,font_face="Arial"),#row_names_gp = gpar(fontsize = 12,fontdace="Arial"),
                  column_names_side = "top",
                  heatmap_legend_param = list(
                    at = c(round(min(log2(heatmap))),log2(1),round(max(log2(heatmap)))),
                    labels=c(format(round(min(heatmap),2),nsmall=2),1,round(max(heatmap))),
                    #                       labels = c("Decreasing pathways with Injury vs Lam","Does not change with injury", "Increasing pathways with Injury vs Lam"),
                    title = "Ratio",
                    legend_height = unit(4, "cm"), legend_width=unit(4,"cm"),
                    title_position = "topcenter"
                  ),
                  row_names_max_width = max_text_width(rownames(heatmap)),
                  column_title_gp = gpar(fontsize = 14, fontface = "bold")
                  )
fig4_t4_2

pdf("fig4_t4_2.pdf", width = 11, height = 6,res = 300)# Specify the dimensions
draw(fig4_t4_2)
dev.off()

png(file="fig4_t4.png",
    width=11, height=6,units = "in",res = 300)
draw(fig4_t4)
dev.off()

