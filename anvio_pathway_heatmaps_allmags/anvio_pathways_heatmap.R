library(ComplexHeatmap)
library(circlize)
#########
setwd("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/anvio_all_mags")

vip_list<-read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/anvio_all_mags/vip_list.txt",
                   sep="\t",header=FALSE) %>%
  rename(bin=V1,species=V2)
VIPs<-
  read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/01_DRAM_VIPs/VIPs_group.txt",sep="\t") %>%
  rename(species=sample) %>%
  inner_join(vip_list)
####maybe switch rows with columns
pathway_presence_matrix_df<- 
  read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/anvio_all_mags/MAGs-presence-MATRIX.txt",
           sep="\t")# %>%
#  pivot_longer(-module, names_to = "bin",values_to = "completeness")

modules_info <- read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/01_DRAM_VIPs/anvio/modules_info.txt",
                         sep="\t") %>%
  select(module,name)
heatmapdf<- pathway_presence_matrix_df %>%
  inner_join(modules_info,.) %>%
  rename(Module=module)
  
#rownames(heatmapdf) <- heatmapdf$Module
heatmap_num<-as.matrix(heatmapdf[, -1:-2])

heatmap_num<-matrix(as.numeric(heatmap_num),nrow = 145,ncol = 263)
colnames(heatmap_num) <-colnames(heatmapdf[, -1:-2])
rownames(heatmap_num) <- heatmapdf$name



highlighted_columns<-data.frame(bin=colnames(heatmapdf))
highlighted_columns$Highlight<-ifelse(highlighted_columns$bin %in% VIPs$bin,"VIP","Not a VIP")
highlighted_columns<-highlighted_columns %>%
  left_join(VIPs) %>%
  replace(is.na(.), "NA")


ha = columnAnnotation(VIP=highlighted_columns$Highlight,
                   col = list(Highlight = c("TRUE" ="#FFBA49", "FALSE" = "#2274A5")))#,

heatmap_obj<-
Heatmap(heatmap_num,na_col = "black",name = "mat",cluster_rows = FALSE,cluster_columns = TRUE,
                     column_names_rot = 90,show_row_names = TRUE,top_annotation = ha,
                     show_column_dend = FALSE,
                     col = c("#F4F4F4","#004346"),
                     width = unit(30, "cm"), 
                     height = unit(30, "cm"),
        row_names_side = "left",
                     column_names_gp = gpar(fontsize = 3), row_names_gp = gpar(fontsize = 4),
        column_names_side = "top",
                     heatmap_legend_param = list(
                       at = c(0,1),
                       labels = c("<75%", ">=75%"),
                       title = "Completeness",
                       legend_height = unit(4, "cm"),
                       title_position = "topcenter"
                     ),
                     row_names_max_width = max_text_width(rownames(heatmap))
)
pdf("all_mags_pathways_presence.pdf", width = 16, height = 16)  # Specify the dimensions
draw(heatmap_obj)
dev.off()

#How similar are these MAGs metabolically?
###########

setwd("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/anvio_all_mags")

vip_list<-read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/anvio_all_mags/vip_list.txt",
                   sep="\t",header=FALSE) %>%
  rename(bin=V1,species=V2)
VIPs<-
  read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/01_DRAM_VIPs/VIPs_group.txt",sep="\t") %>%
  rename(species=sample) %>%
  inner_join(vip_list) %>%
  # mutate_all(~gsub("CAG.","CAG-",.)) %>%
  # mutate_all(~gsub("X1XD8.7","1XD8-7",.)) %>%
  select(-species)
checkM<- read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/CheckM_galah95.txt",sep="\t") %>%
  rename(bin=Bin.Id) %>%
  select(bin,Completeness,Contamination)
gtdb<- read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/GTDBtk/Clean_classification_dereplicated_95_taxonomy.tsv",
                sep="\t") %>%
  select(MAG,phylum,fam_species) %>%
  rename(bin=MAG) %>%
  mutate(bin=gsub("_bin.","_bin",bin))


####maybe switch rows with columns
pathway_presence_matrix_df<- 
  read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/anvio_all_mags/MAGs-presence-MATRIX.txt",
           sep="\t") %>%
  pivot_longer(-module, names_to = "bin",values_to = "completeness")

heatmapdf <- read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/01_DRAM_VIPs/anvio/modules_info.txt",
                      sep="\t") %>%
  select(module,name) %>%
  inner_join(.,pathway_presence_matrix_df) %>%
  select(-module) %>%
  pivot_wider(id_cols = "bin",names_from = "name",values_from = "completeness") %>%
  inner_join(gtdb,.) %>%
  left_join(VIPs) %>%
  replace(is.na(.), "NA") %>%
  inner_join(checkM,.)
heatmapdf$VIP<-ifelse(heatmapdf$bin %in% VIPs$bin,"VIP","Not a VIP")
#rownames(heatmapdf) <- heatmapdf$bin
heatmap<- heatmapdf %>%
  select(-c(group,VIP,phylum,bin,fam_species,Completeness,Contamination)) %>%
  as.matrix()
rownames(heatmap) <- heatmapdf$fam_species


col_fun_compl = colorRamp2(c(min(heatmapdf$Completeness),100), c("#fdf0d5", "#780000"))
col_fun_cont=colorRamp2(c(min(heatmapdf$Contamination),max(heatmapdf$Contamination)), c("#FFFFFF", "#c1121f"))
#col_fun=colorRamp2(, c("#FFFFFF", "#c1121f"))
heatmapdf$phylum %>%
  unique

ha = rowAnnotation(VIP=heatmapdf$VIP,
                   Direction=heatmapdf$group,
                   MAG_Completeness=heatmapdf$Completeness,
                   MAG_Contamination=heatmapdf$Contamination,
                   Phylum=heatmapdf$phylum,
                   col = list(VIP = c("VIP" ="#780000", "Not a VIP" = "#EEEEEE"),
                              Direction = c("Increasing with injury" ="#fb8500", "Decreasing with injury" = "#219ebc",
                                        "NA"="#EEEEEE"),
                              MAG_Completeness=col_fun_compl,
                              MAG_Contamination=col_fun_cont,
                              Phylum=c("p__Firmicutes"="#ffffff",
                                       "p__Actinobacteriota"="#b08968",
                                       "p__Firmicutes_A"="#fca311",
                                       "p__Verrucomicrobiota"="#b5e48c",
                                       "p__Proteobacteria"="#3a86ff",
                                       "p__Bacteroidota"="#000000",
                                       "p__Firmicutes_B"="#0a9396")
                              ))#,
# col = list(Group = c("Increasing_with_injury" ="#F46036", "Decreasing_with_injury" = "#44CF6C")))#,

heatmap_obj<-Heatmap(heatmap,na_col = "black",name = "mat",cluster_rows = TRUE,cluster_columns = TRUE,
                     column_names_rot = 70,show_row_names = TRUE,right_annotation = ha,
                     show_column_dend = FALSE,show_row_dend = FALSE,
                     col = c("#F4F4F4","#004346"),
                     width = unit(30, "cm"), 
                     height = unit(30, "cm"),
                     row_names_side = "left",
                     column_names_gp = gpar(fontsize = 4), row_names_gp = gpar(fontsize = 4),
                     column_names_side = "top",
                     heatmap_legend_param = list(
                       at = c(0,1),
                       labels = c("<75%", ">=75%"),
                       title = "Module\n\t\tCompleteness",
                       legend_height = unit(4, "cm"),
                       title_position = "topcenter"
                     ),
                     row_names_max_width = max_text_width(rownames(heatmap))
)
pdf("all_mags_rows_pathways_presence.pdf", width = 17, height = 17)  # Specify the dimensions
draw(heatmap_obj)
dev.off()

png(file="all_mags_rows_pathways_presence.png",
    width=17, height=17,units = "in",res = 300)
draw(heatmap_obj)
dev.off()



#pathways abundance
#male abundance of pathways using mags abundance as proxy
#########################
modules_info <- read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/01_DRAM_VIPs/anvio/modules_info.txt",
                         sep="\t") %>%
  select(module,name)

pathway_presence_matrix_df<- 
  read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/anvio_all_mags/MAGs-presence-MATRIX.txt",
           sep="\t") %>%
  pivot_longer(-module, names_to = "bin",values_to = "completeness") %>%
  inner_join(modules_info) %>%
  select(-module) %>%
  rename(module=name)

m_metadata <- read_csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/Results/metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Type=="Fecal" & Sex_SurgeryDate!="M_Day1" & Sex_SurgeryDate!="F_Day4") %>%
  mutate(Group=factor(Group,levels=c("Lam", "T10","T4")),
         Timepoint=factor(Timepoint,
                          levels=c("0dpi", "7dpi","21dpi","35dpi","63dpi","6_months"),
                          labels=c("0dpi","7dpi","21dpi", "35dpi","63dpi","180dpi"))) %>%
  filter(Sex=="M")

m_tm<- read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/CoverM_genomes/Normalized_Trimmed_Mean.csv") %>%
  # pivot_longer(-Sample, names_to = "mag", values_to = "tm") %>%
  # group_by(Sample) %>%
  # mutate(rel_abund= 100* tm/sum(tm)) %>%
  # ungroup() %>%
  # select(-tm) %>%
  # pivot_wider(names_from = "mag",values_from = "rel_abund") %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  inner_join(m_metadata,.) %>%
  select(-c(3:6,8:13)) %>%
  pivot_longer(-c(Sample,Group,Timepoint),names_to = "bin",values_to = "abundance") %>%
  group_by(bin,Group,Timepoint) %>%
  summarise(median_tm=median(abundance)) %>% #median across mice
  ungroup() %>%
  inner_join(pathway_presence_matrix_df) %>%
  mutate(module_mag=completeness*median_tm) %>%
  select(-c(bin,completeness,median_tm)) %>%
  group_by(module,Group,Timepoint) %>%
  mutate(sum_median=sum(module_mag)+1) %>% #sum of pathway across all mags
  ungroup() %>%
  select(-module_mag)
#  mutate(log_sum_median=log(sum_median)) %>% #log to reduce variability
#  select(-c(sum_median,module_mag))

m_metadata0<- m_metadata %>%
  filter(Timepoint=="0dpi")
abund_long_0dpi<- m_tm %>%
  inner_join(m_metadata) %>%
  filter(Timepoint=="0dpi") %>%
  select(module,Group,Timepoint,sum_median) %>%
  rename(sum_median_0=sum_median) %>%
  rename(module_0=module) %>%
  distinct()

abund_long_folds<-left_join(m_tm,abund_long_0dpi,by=c("Group","module"="module_0")) %>%
  mutate(sum_median_fold=sum_median/sum_median_0) %>%
  select(-c(sum_median,sum_median_0,Timepoint.y)) %>%
  rename(Timepoint=Timepoint.x) %>%
  distinct
#### Trying statistical tests
abund_long_fold<-left_join(m_tm,abund_long_0dpi,by=c("Group","module"="module_0")) %>%
  mutate(sum_median_fold=sum_median/sum_median_0) %>%
  select(-c(sum_median,sum_median_0,Timepoint.y)) %>%
  rename(Timepoint=Timepoint.x) %>%
  distinct

pairwise_fold<- abund_long_fold %>%
#  filter(Group=="T10") %>%
  filter(Timepoint!="0dpi") %>%
  nest(data = -module) %>%
  mutate(pairwise_tests = map(.x=data,
                              ~pairwise.wilcox.test(x=.x$sum_median_fold,
                                                    g=.x$Group,
                                                    p.adjust.method = "bonferroni",exact=FALSE) %>%
                                tidy())) %>%
  unnest(pairwise_tests) %>%
  filter(p.value < 0.05)# %>%
#  unnest(data)
#Resulting in 12 pathways. 
####
####
lam<-abund_long_folds %>%
  filter(Group=="Lam") %>%
  rename(sum_median_fold_lam=sum_median_fold)
abund_fold <- abund_long_folds %>%
  inner_join(lam,by=c("Timepoint","module")) %>%
  mutate(change_sci_over_lam=sum_median_fold/sum_median_fold_lam) %>%
  filter(Group.x!="Lam") %>%
  select(-c(Group.y,sum_median_fold,sum_median_fold_lam)) %>%
  rename(Group=Group.x) %>%
  pivot_wider(id_cols = c(Group,module),names_from = Timepoint,values_from = change_sci_over_lam)
t10<- abund_fold %>%
  filter(Group=="T10") %>%
  select(-Group)# %>%
#   group_by(module) %>%
# #  summarise(across(everything(), sum)) %>%
#   mutate(sum_of_values = rowSums(across(everything()))) %>%
#   ungroup() %>%
#   filter(sum_of_values!=6) %>%
#   select(-sum_of_values)
#  summarize(mutate(sum=sum(c(get('7dpi'),get('21dpi'),get('35dpi'),get('63dpi'),get('180dpi')))))
  
t4<- abund_fold %>%
  filter(Group=="T4") %>%
  select(-Group)



heatmapdf<- t10 %>%
  rename(Module=module)%>%
  mutate(Total = rowSums(select(., -Module))) %>%
  filter(Total!=6) %>%
  select(-Total) %>%
 mutate(Median = apply(select(., -Module), 1, median)) %>%
 filter(Median <0.9 | Median>1.1) %>%
 select(-Median)

  
heatmap<-as.matrix(heatmapdf[, -1])
rownames(heatmap) <- heatmapdf$Module


enriched_modules<-
  read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/01_DRAM_VIPs/anvio/VIPs_enriched_modules.txt",sep = "\t") %>%
  filter(adjusted_q_value<=0.05 | (p_Decreasing_with_injury>=0 & p_Increasing_with_injury==0) | (p_Increasing_with_injury>=0 & p_Decreasing_with_injury==0)) %>%
  select(KEGG_MODULE,associated_groups) %>%
  rename(module=KEGG_MODULE)

highlighted_columns<-data.frame(module=rownames(heatmap))
#highlighted_columns$Highlight<-ifelse(highlighted_columns$module %in% enriched_modules$module,"Enriched in a VIP","Not enriched in a VIP")
highlighted_columns<-highlighted_columns %>%
  left_join(enriched_modules,by = c()) %>%
  replace(is.na(.), "NA")


ha = rowAnnotation(VIP_Enrichment=highlighted_columns$Highlight,
                   Associated_group=highlighted_columns$associated_groups,
                   col = list(#VIP_Enrichment = c("Enriched in a VIP" ="#C83E4D", 
#                                            "Not enriched in a VIP" = "#FFFFFF"),
                              Associated_group = c("Decreasing_with_injury" ="#219ebc", 
                                            "Increasing_with_injury" = "#fb8500",
                                            "NA"="#FFFFFF")))#,

heatmap_obj<-Heatmap(log10(heatmap),na_col = "black",name = "T10",column_title = "T10/lam",
                     cluster_rows = TRUE,cluster_columns = FALSE, right_annotation = ha,
                     column_names_rot = 90,show_row_names = TRUE,
                     show_column_dend = FALSE,
                     col = c("#219ebc","#FFFFFF","#fb8500"),
                     width = unit(5, "cm"), 
                     height = unit(30, "cm"),
                     row_names_side = "left",
                     column_names_gp = gpar(fontsize = 12), row_names_gp = gpar(fontsize = 4),
                     column_names_side = "top",
                     heatmap_legend_param = list(
#                       at = c(min(log10(heatmap)),log10(1),max(log10(heatmap))),
#                       labels = c("Decreasing modules with injury","Not changing", "Increasing modules with injury"),
                       at = c(log10(min(heatmap)),log10(1),log10(max(heatmap))),
                       labels = c(min(heatmap),1,max(heatmap)),
                       title = "Module's median fold change \nacross mice after injury \ncompared to Lam\n",
                       legend_height = unit(4, "cm"),
                       title_position = "topcenter"
                     ),
                     row_names_max_width = max_text_width(rownames(heatmap)),
        column_title_gp = gpar(fontsize = 14, fontface = "bold")
)
heatmap_obj
pdf("t10_all_mags_pathways_presence_wMedianFilter.pdf", width = 8, height = 14.5)  # Specify the dimensions
draw(heatmap_obj)
dev.off()
png(file="t10_all_mags_pathways_presence.png",
    width=8, height=13,units = "in",res = 300)
draw(heatmap_obj)
dev.off()


### t4

heatmapdf<- t4 %>%
  select(-c("180dpi")) %>%
  rename(Module=module)%>%
  mutate(Total = rowSums(select(., -Module))) %>%
  filter(Total!=5) %>%
  select(-Total) %>%
  mutate(Median = apply(select(., -Module), 1, median)) %>%
#  filter(Median <0.9 | Median>1.1) %>%
  select(-Median)


heatmap<-as.matrix(heatmapdf[, -1])
rownames(heatmap) <- heatmapdf$Module


enriched_modules<-
  read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/01_DRAM_VIPs/anvio/VIPs_enriched_modules.txt",sep = "\t") %>%
  filter(adjusted_q_value<=0.05 | (p_Decreasing_with_injury>=0 & p_Increasing_with_injury==0) | (p_Increasing_with_injury>=0 & p_Decreasing_with_injury==0)) %>%
  select(KEGG_MODULE,associated_groups) %>%
  rename(module=KEGG_MODULE)

highlighted_columns<-data.frame(module=rownames(heatmap))
#highlighted_columns$Highlight<-ifelse(highlighted_columns$module %in% enriched_modules$module,"Enriched in a VIP","Not enriched in a VIP")
highlighted_columns<-highlighted_columns %>%
  left_join(enriched_modules,by = c()) %>%
  replace(is.na(.), "NA")


ha = rowAnnotation(VIP_Enrichment=highlighted_columns$Highlight,
                   Associated_group=highlighted_columns$associated_groups,
                   col = list(#VIP_Enrichment = c("Enriched in a VIP" ="#C83E4D", 
                     #                                            "Not enriched in a VIP" = "#FFFFFF"),
                     Associated_group = c("Decreasing_with_injury" ="#219ebc", 
                                          "Increasing_with_injury" = "#fb8500",
                                          "NA"="#FFFFFF")))#,

heatmap_obj<-Heatmap(log10(heatmap),na_col = "black",name = "T4",column_title = "T4/lam",
                     cluster_rows = TRUE,cluster_columns = FALSE, right_annotation = ha,
                     column_names_rot = 90,show_row_names = TRUE,
                     show_column_dend = FALSE,
                     col = c("#219ebc","#FFFFFF","#fb8500"),
                     width = unit(5, "cm"), 
                     height = unit(30, "cm"),
                     row_names_side = "left",
                     column_names_gp = gpar(fontsize = 12), row_names_gp = gpar(fontsize = 4),
                     column_names_side = "top",
                     heatmap_legend_param = list(
                       #                       at = c(min(log10(heatmap)),log10(1),max(log10(heatmap))),
                       #                       labels = c("Decreasing modules with injury","Not changing", "Increasing modules with injury"),
                       at = c(log10(min(heatmap)),log10(1),log10(max(heatmap))),
                       labels = c(min(heatmap),1,max(heatmap)),
                       title = "Module's median fold change \nacross mice after injury \ncompared to Lam\n",
                       legend_height = unit(4, "cm"),
                       title_position = "topcenter"
                     ),
                     row_names_max_width = max_text_width(rownames(heatmap)),
                     column_title_gp = gpar(fontsize = 14, fontface = "bold")
)
heatmap_obj
pdf("t4_all_mags_pathways_presence_woMedianFilter.pdf", width = 8, height = 14.5)  # Specify the dimensions
draw(heatmap_obj)
dev.off()


#Female?
###############
modules_info <- read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/01_DRAM_VIPs/anvio/modules_info.txt",
                         sep="\t") %>%
  select(module,name)

pathway_presence_matrix_df<- 
  read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/anvio_all_mags/MAGs-presence-MATRIX.txt",
           sep="\t") %>%
  pivot_longer(-module, names_to = "bin",values_to = "completeness") %>%
  inner_join(modules_info) %>%
  select(-module) %>%
  rename(module=name)

f_metadata <- read_csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/Results/metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Type=="Fecal" & Sex_SurgeryDate!="M_Day1" & Sex_SurgeryDate!="F_Day4") %>%
  mutate(Group=factor(Group,levels=c("Lam", "T10","T4")),
         Timepoint=factor(Timepoint,
                          levels=c("0dpi", "7dpi","21dpi","35dpi","63dpi","6_months"),
                          labels=c("0dpi","7dpi","21dpi", "35dpi","63dpi","180dpi"))) %>%
  filter(Sex=="F")

f_tm<- read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/CoverM_genomes/Normalized_Trimmed_Mean.csv") %>%
  # pivot_longer(-Sample, names_to = "mag", values_to = "tm") %>%
  # group_by(Sample) %>%
  # mutate(rel_abund= 100* tm/sum(tm)) %>%
  # ungroup() %>%
  # select(-tm) %>%
  # pivot_wider(names_from = "mag",values_from = "rel_abund") %>%
  filter(!str_detect(Sample,"Repeat")) %>%
  inner_join(f_metadata,.) %>%
  select(-c(3:6,8:13)) %>%
  pivot_longer(-c(Sample,Group,Timepoint),names_to = "bin",values_to = "abundance") %>%
  group_by(bin,Group,Timepoint) %>%
  summarise(median_tm=median(abundance)) %>% #median across mice
  ungroup() %>%
  inner_join(pathway_presence_matrix_df) %>%
  mutate(module_mag=completeness*median_tm) %>%
  select(-c(bin,completeness,median_tm)) %>%
  group_by(module,Group,Timepoint) %>%
  mutate(sum_median=sum(module_mag)+1) %>% #sum of pathway across all mags
  ungroup() %>%
  select(-module_mag)
#  mutate(log_sum_median=log(sum_median)) %>% #log to reduce variability
#  select(-c(sum_median,module_mag))

f_metadata0<- f_metadata %>%
  filter(Timepoint=="0dpi")
abund_long_0dpi<- f_tm %>%
  inner_join(f_metadata) %>%
  filter(Timepoint=="0dpi") %>%
  select(module,Group,Timepoint,sum_median) %>%
  rename(sum_median_0=sum_median) %>%
  rename(module_0=module) %>%
  distinct()

abund_long_folds<-left_join(f_tm,abund_long_0dpi,by=c("Group","module"="module_0")) %>%
  mutate(sum_median_fold=sum_median/sum_median_0) %>%
  select(-c(sum_median,sum_median_0,Timepoint.y)) %>%
  rename(Timepoint=Timepoint.x) %>%
  distinct
#### Trying statistical tests
abund_long_fold<-left_join(f_tm,abund_long_0dpi,by=c("Group","module"="module_0")) %>%
  mutate(sum_median_fold=sum_median/sum_median_0) %>%
  select(-c(sum_median,sum_median_0,Timepoint.y)) %>%
  rename(Timepoint=Timepoint.x) %>%
  distinct

pairwise_fold<- abund_long_fold %>%
  #  filter(Group=="T10") %>%
  filter(Timepoint!="0dpi") %>%
  nest(data = -module) %>%
  mutate(pairwise_tests = map(.x=data,
                              ~pairwise.wilcox.test(x=.x$sum_median_fold,
                                                    g=.x$Group,
                                                    p.adjust.method = "bonferroni",exact=FALSE) %>%
                                tidy())) %>%
  unnest(pairwise_tests) %>%
  filter(p.value < 0.05)# %>%
#  unnest(data)
#Resulting in 12 pathways. 
####
####
lam<-abund_long_folds %>%
  filter(Group=="Lam") %>%
  rename(sum_median_fold_lam=sum_median_fold)
abund_fold <- abund_long_folds %>%
  inner_join(lam,by=c("Timepoint","module")) %>%
  mutate(change_sci_over_lam=sum_median_fold/sum_median_fold_lam) %>%
  filter(Group.x!="Lam") %>%
  select(-c(Group.y,sum_median_fold,sum_median_fold_lam)) %>%
  rename(Group=Group.x) %>%
  pivot_wider(id_cols = c(Group,module),names_from = Timepoint,values_from = change_sci_over_lam)
t10<- abund_fold %>%
  filter(Group=="T10") %>%
  select(-Group)# %>%
#   group_by(module) %>%
# #  summarise(across(everything(), sum)) %>%
#   mutate(sum_of_values = rowSums(across(everything()))) %>%
#   ungroup() %>%
#   filter(sum_of_values!=6) %>%
#   select(-sum_of_values)
#  summarize(mutate(sum=sum(c(get('7dpi'),get('21dpi'),get('35dpi'),get('63dpi'),get('180dpi')))))

t4<- abund_fold %>%
  filter(Group=="T4") %>%
  select(-Group)



heatmapdf<- t10 %>%
  rename(Module=module)%>%
  mutate(Total = rowSums(select(., -Module))) %>%
  filter(Total!=6) %>%
  select(-Total) %>%
  mutate(Median = apply(select(., -Module), 1, median)) %>%
#  filter(Median <0.9 | Median>1.1) %>%
  select(-Median)


heatmap<-as.matrix(heatmapdf[, -1])
rownames(heatmap) <- heatmapdf$Module


enriched_modules<-
  read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/01_DRAM_VIPs/anvio/VIPs_enriched_modules.txt",sep = "\t") %>%
  filter(adjusted_q_value<=0.05 | (p_Decreasing_with_injury>=0 & p_Increasing_with_injury==0) | (p_Increasing_with_injury>=0 & p_Decreasing_with_injury==0)) %>%
  select(KEGG_MODULE,associated_groups) %>%
  rename(module=KEGG_MODULE)

highlighted_columns<-data.frame(module=rownames(heatmap))
#highlighted_columns$Highlight<-ifelse(highlighted_columns$module %in% enriched_modules$module,"Enriched in a VIP","Not enriched in a VIP")
highlighted_columns<-highlighted_columns %>%
  left_join(enriched_modules,by = c()) %>%
  replace(is.na(.), "NA")


ha = rowAnnotation(VIP_Enrichment=highlighted_columns$Highlight,
                   Associated_group=highlighted_columns$associated_groups,
                   col = list(#VIP_Enrichment = c("Enriched in a VIP" ="#C83E4D", 
                     #                                            "Not enriched in a VIP" = "#FFFFFF"),
                     Associated_group = c("Decreasing_with_injury" ="#219ebc", 
                                          "Increasing_with_injury" = "#fb8500",
                                          "NA"="#FFFFFF")))#,

heatmap_obj<-Heatmap(log10(heatmap),na_col = "black",name = "T10",column_title = "T10/lam",
                     cluster_rows = TRUE,cluster_columns = FALSE, right_annotation = ha,
                     column_names_rot = 90,show_row_names = TRUE,
                     show_column_dend = FALSE,
                     col = c("#219ebc","#FFFFFF","#fb8500"),
                     width = unit(5, "cm"), 
                     height = unit(30, "cm"),
                     row_names_side = "left",
                     column_names_gp = gpar(fontsize = 12), row_names_gp = gpar(fontsize = 4),
                     column_names_side = "top",
                     heatmap_legend_param = list(
                       #                       at = c(min(log10(heatmap)),log10(1),max(log10(heatmap))),
                       #                       labels = c("Decreasing modules with injury","Not changing", "Increasing modules with injury"),
                       at = c(log10(min(heatmap)),log10(1),log10(max(heatmap))),
                       labels = c(min(heatmap),1,max(heatmap)),
                       title = "Module's median fold change \nacross mice after injury \ncompared to Lam\n",
                       legend_height = unit(4, "cm"),
                       title_position = "topcenter"
                     ),
                     row_names_max_width = max_text_width(rownames(heatmap)),
                     column_title_gp = gpar(fontsize = 14, fontface = "bold")
)
heatmap_obj
pdf("t10_all_mags_pathways_presence_female_woMedianFilter.pdf", width = 8, height = 14.5)  # Specify the dimensions
draw(heatmap_obj)
dev.off()


### t4


heatmapdf<- t4 %>%
  select(-c("180dpi")) %>%
  rename(Module=module)%>%
  mutate(Total = rowSums(select(., -Module))) %>%
  filter(Total!=5) %>%
  select(-Total) %>%
  mutate(Median = apply(select(., -Module), 1, median)) %>%
  filter(Median <0.9 | Median>1.1) %>%
  select(-Median)


heatmap<-as.matrix(heatmapdf[, -1])
rownames(heatmap) <- heatmapdf$Module


enriched_modules<-
  read.csv("/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/01_DRAM_VIPs/anvio/VIPs_enriched_modules.txt",sep = "\t") %>%
  filter(adjusted_q_value<=0.05 | (p_Decreasing_with_injury>=0 & p_Increasing_with_injury==0) | (p_Increasing_with_injury>=0 & p_Decreasing_with_injury==0)) %>%
  select(KEGG_MODULE,associated_groups) %>%
  rename(module=KEGG_MODULE)

highlighted_columns<-data.frame(module=rownames(heatmap))
#highlighted_columns$Highlight<-ifelse(highlighted_columns$module %in% enriched_modules$module,"Enriched in a VIP","Not enriched in a VIP")
highlighted_columns<-highlighted_columns %>%
  left_join(enriched_modules,by = c()) %>%
  replace(is.na(.), "NA")


ha = rowAnnotation(VIP_Enrichment=highlighted_columns$Highlight,
                   Associated_group=highlighted_columns$associated_groups,
                   col = list(#VIP_Enrichment = c("Enriched in a VIP" ="#C83E4D", 
                     #                                            "Not enriched in a VIP" = "#FFFFFF"),
                     Associated_group = c("Decreasing_with_injury" ="#219ebc", 
                                          "Increasing_with_injury" = "#fb8500",
                                          "NA"="#FFFFFF")))#,

heatmap_obj<-Heatmap(log10(heatmap),na_col = "black",name = "T4",column_title = "T4/lam",
                     cluster_rows = TRUE,cluster_columns = FALSE, right_annotation = ha,
                     column_names_rot = 90,show_row_names = TRUE,
                     show_column_dend = FALSE,
                     col = c("#219ebc","#FFFFFF","#fb8500"),
                     width = unit(5, "cm"), 
                     height = unit(30, "cm"),
                     row_names_side = "left",
                     column_names_gp = gpar(fontsize = 12), row_names_gp = gpar(fontsize = 4),
                     column_names_side = "top",
                     heatmap_legend_param = list(
                       #                       at = c(min(log10(heatmap)),log10(1),max(log10(heatmap))),
                       #                       labels = c("Decreasing modules with injury","Not changing", "Increasing modules with injury"),
                       at = c(log10(min(heatmap)),log10(1),log10(max(heatmap))),
                       labels = c(min(heatmap),1,max(heatmap)),
                       title = "Module's median fold change \nacross mice after injury \ncompared to Lam\n",
                       legend_height = unit(4, "cm"),
                       title_position = "topcenter"
                     ),
                     row_names_max_width = max_text_width(rownames(heatmap)),
                     column_title_gp = gpar(fontsize = 14, fontface = "bold")
)
heatmap_obj

pdf("t4_all_mags_pathways_presence_female_wMedianFilter.pdf", width = 8, height = 14.5)  # Specify the dimensions
draw(heatmap_obj)
dev.off()
png(file="t4_all_mags_pathways_presence_female.png",
    width=8, height=13,units = "in",res = 300)
draw(heatmap_obj)
dev.off()
