library(ComplexHeatmap)

setwd("~/anvio_pathway_heatmaps_allmags") #set to this directory


modules_info <- read.csv("modules_info.txt",
                         sep="\t") %>%
  select(module,name)

pathway_presence_matrix_df<- 
  read.csv("MAGs-presence-MATRIX.txt",
           sep="\t") %>%
  pivot_longer(-module, names_to = "bin",values_to = "completeness") %>%
  inner_join(modules_info) %>%
  select(-module) %>%
  rename(module=name)

m_metadata <- read_csv("metadata_w_surgerydate.csv") %>% 
  filter(!str_detect(Sample,"Repeat")) %>%
  filter(Type=="Fecal" & Sex_SurgeryDate!="M_Day1" & Sex_SurgeryDate!="F_Day4") %>%
  mutate(Group=factor(Group,levels=c("Lam", "T10","T4")),
         Timepoint=factor(Timepoint,
                          levels=c("0dpi", "7dpi","21dpi","35dpi","63dpi","6_months"),
                          labels=c("0dpi","7dpi","21dpi", "35dpi","63dpi","180dpi"))) %>%
  filter(Sex=="M")

m_tm<- read.csv("Normalized_Trimmed_Mean.csv") %>%
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
  select(-Group)


abr_heatmapdf<- t10 %>%
  rename(Module=module)%>%
  mutate(Total = rowSums(select(., -Module))) %>%
  filter(Total!=6) %>%
  select(-Total) %>%
  mutate(Median = apply(select(., -Module), 1, median)) %>%
  filter(Median <0.42 | Median>1.1) %>%
  select(-Median) %>%
  filter(Module!="Flavone degradation, luteolin/apigenin => DHCA/phloretate")


heatmap<-as.matrix(abr_heatmapdf[, -1])
rownames(heatmap) <- abr_heatmapdf$Module


enriched_modules<-
  read.csv("VIPs_enriched_modules.txt",sep = "\t") %>%
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

heatmap_obj<-Heatmap(heatmap^0.15,na_col = "black",name = "T10",column_title = "T10/lam",
                     cluster_rows = TRUE,cluster_columns = FALSE, right_annotation = ha,
                     column_names_rot = 90,show_row_names = TRUE,
                     show_column_dend = FALSE,
                     col = c("#219ebc","#FFFFFF","#fb8500"),
                     width = unit(5, "cm"), 
                     height = unit(10, "cm"),
                     row_names_side = "left",
                     column_names_gp = gpar(fontsize = 12), row_names_gp = gpar(fontsize = 9),
                     column_names_side = "top",
                     heatmap_legend_param = list(
                       #                       at = c(min(log10(heatmap)),log10(1),max(log10(heatmap))),
                       #                       labels = c("Decreasing modules with injury","Not changing", "Increasing modules with injury"),
                       at = c(min(heatmap)^0.15,1,max(heatmap)^0.15),
                       labels = c(round(min(heatmap),2),1,round(max(heatmap),2)),
                       title = "Module's median fold change \nacross mice after injury \ncompared to Lam\n",
                       legend_height = unit(4, "cm"),
                       title_position = "topcenter"
                     ),
                     row_names_max_width = max_text_width(rownames(heatmap)),
                     column_title_gp = gpar(fontsize = 14, fontface = "bold")
)
heatmap_obj


pdf("abriged_m_t10_all_mags_pathways_presence_wMedianFilter.pdf", width = 11, height = 6.5)  # Specify the dimensions
draw(heatmap_obj)
dev.off()