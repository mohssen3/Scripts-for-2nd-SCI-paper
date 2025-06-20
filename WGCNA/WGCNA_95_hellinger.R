# load the required libraries
library("WGCNA")
library("vegan")
library("RColorBrewer")
library("pls")
library("ggplot2")
library("plyr")
library("spls")
library("dplyr")
library(stringr)
library(tidyverse)
#enableWGCNAThreads()
disableWGCNAThreads()
#test<-read.csv("../Meta.csv")
# set the working directory

abs_path="/fs/ess/PAS2210/Mohamed2020/04_MAGs/00_Rebinning/unitem/Results/Scripts/WGCNA/"
setwd(abs_path)

Meta <- read_csv("../Meta.csv") %>%
  rename("Mouse"="Mouse_num") %>%
  filter(!str_detect(Sample,"Repeat"))
abund<- read_csv("../95_Normalized_Trimmed_Mean.csv") %>%
#  rename("PLATE"="Sample") %>%
#  select(contains(c("bin","maxbin","Sample"))) %>%
  filter(!str_detect(Sample,"Repeat"))
joined <- left_join(Meta,abund,by="Sample")
taxonomy<- read.csv("../95_Clean_classification_dereplicated_taxonomy.tsv",sep="\t") %>%
  select(c("MAG","phylum","class","order","family","genus","species"))
#sex="M";value="1";ColumnName="Group";day_excluded="1";type="Fecal"
filtering_to_picking_threshold <- function(sex,value,ColumnName,day_excluded,type="Fecal") {
  
  setwd(abs_path)
  # import the abundance table (the dataset and associated metadata were reterived from the Bowman Lab's page: www.polarmicrobes.org/weighted-gene-correlation-network-analysis-wgcna-applied-to-microbial-communities/)
  # input data is a OTU table with OTUs as rows and samples as columns, and a metadata file organized in the same order of samples as the OTU table but with samples as rows. 
  #    select(contains(c("bin","maxbin","Sample"))) %>%
  #    filter(!str_detect(Sample,"Repeat")) %>%
  #    column_to_rownames("Sample")
  
  # Sex and surgery day filtration
  data1 <- as.data.frame(joined)
  data1 <- dplyr::filter(data1, grepl(type, Type))
  data1 <- dplyr::filter(data1, !grepl('Repeat_', Sample))
  data1 <- dplyr::filter(data1, grepl(sex, Sex_SurgeryDate))
  data1 <- dplyr::filter(data1, !grepl(day_excluded, Sex_SurgeryDate))
  #data1 <- dplyr::filter(data1, grepl('F', Sex_SurgeryDate))
  #data1 <- dplyr::filter(data1, grepl('1|2|3', Sex_SurgeryDate))
  
  # time-point filtration
  data1 <- dplyr::filter(data1, grepl(value, get(ColumnName)))
  
  data1<-data1 %>% select(-c(2:14)) %>%
    column_to_rownames("Sample")

  dir1=paste0(sex,"_Day",day_excluded,"Excluded_",value)
  if (str_detect(type,"Cecum")) {
    dir1=paste0(sex,"_Day",day_excluded,"Excluded_",value,"_FC")
  }
  
  dir.create(dir1)
  setwd(dir1)
  
  write.csv(data1,paste0("Normalized_TM",dir1,".csv"))
  transformed_abund <- decostand(data1, "hellinger")

    
  gsg = goodSamplesGenes(transformed_abund, verbose = 3);
  gsg$allOK
  
  #############################################################
  ### if NOT TRUE, we remove the offending OTUs and samples from the data (by running the function below)
  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(transformed_abund)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(transformed_abund)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    transformed_abund = transformed_abund[gsg$goodSamples, gsg$goodGenes]
  }
  
  # check again that there is not too much missing data
  gsg = goodSamplesGenes(transformed_abund, verbose = 3);
  gsg$allOK
  write.csv(transformed_abund,paste0("hellinger_transformed_TM_",dir1,".csv"))
  
  # clustering SAMPLES (not OTUs) by their similarities in OTU abundance, 
  clustTree = hclust(dist(transformed_abund), method = "average");
  
  # Plot the tree: Open a graphic output window of size 12 by 9 inches
  sizeGrWindow(12,9)
  pdf(file = paste0("OTU_module_clusters_loc_",dir1,".pdf"), width = 12, height = 9);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(clustTree, main = "OTU module clusters", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  dev.off()
  
  # now we want to see how the OTUs relate to each other in terms of abundance across samples (build the network), then
  # we apply a series of power transformations on the edges of the network till we fit the scale free topology model of the network. 
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  # Call the network topology analysis function, this takes 30ish seconds
  sft = pickSoftThreshold(transformed_abund, powerVector = powers, verbose = 5, networkType = "signed")
  # Plot the results:
  sizeGrWindow(9, 5)
  pdf(paste0("Powers_",dir1,".pdf"),width=12,height=8)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit as a function of the soft-thresholding power
  # Scale free topology - imagine a group of friends (group A) go to a party, each group A person invites a number of their frends (multiple group Bs). If the room is small (little distance between people) the group B people will make multiple connections with other group B and A people, this leads to a random network. 
  # If we scale up the size of the room, there may be fewer and fewer interactions between different group B people as the distance between groups grow. At some point there will be no new connections made between different group B people because they are sticking with their group A person in various locations of the room, regardelss of the size of the room. This is scale free topology - no new topology regardelss of scale.
  # If we want to identify groups, Its helpfull to scale the data until no new network connections are made and each cluster is distinct. This is what the soft threshold power does.
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of but is somewhat arbitray, can do visually so removed
  # abline(h=0.9,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  # we select the first power at which mean connectivity remains constant
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  
  
  setwd(abs_path)
  return(transformed_abund)
}


VIP_fn <- function(module,parameter,Metadata,pops,nSamples,moduleColors,
                   transformed_abund,sex,value,ColumnName,day_excluded,type="Fecal"){
  
  weight <- as.data.frame(Metadata[,parameter])
  weight_t<- as.data.frame(weight)
  rownames(weight_t) <- rownames(transformed_abund)
#  print(weight_t)
  ###Problem is, very often we have missing data in our metadata file. The more samples you get, the more likely this will happen. R doesn't like
  #this, so we will have to clean up the data
  weight_t <-na.omit(weight_t)
#  print(weight_t)
  cleannamesweight <- rownames(weight_t)
#  print(weight_t)
  names(weight_t) = parameter
  weight <- as.data.frame(weight_t)
  MB_clean <- pops[rownames(weight_t),]
  transformed_abund_clean <- transformed_abund[rownames(weight_t),]
#  transformed_abund_clean
  #transformed_abund <- transformed_abund_clean
  
  ##check they all have the same nr of rows
  nrow(MB_clean)
  nrow(transformed_abund_clean)
  nrow(weight)
  
  ###Now we should be able to come back to our analysis as normal
  modNames = substring(names(pops), 3)
  geneModuleMembership = as.data.frame(cor(transformed_abund, pops, use = "p"));
#  print(geneModuleMembership)
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  geneTraitSignificance = as.data.frame(cor(transformed_abund_clean, weight, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
  names(GSPvalue) = paste("p.GS.", names(weight), sep="");
  
  
  # Then look at the specific module of interest
  column = match(module, modNames);
  moduleGenes = moduleColors==module
  pdf(paste(module,"Module membership vs. population significance for", parameter,".pdf"),width=12,height=8)
  #  par(mfrow=c(1,1))
  
  
  
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("OTU significance for ",parameter),
                     main = paste("Module membership vs. population significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")
  dev.off()
  test<-data.frame(row.names(geneTraitSignificance[moduleGenes,0]),abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]))
  names(test) <- c("MAG","module membership","OTU significance")
  
  write.csv(test,file=paste(module,"Module membership vs. population significance for", parameter,".csv"))
  gp1 <- test %>% ggplot(aes(x="module membership", y= "OTU significance"))+
    geom_point()
  labs(x=paste("Module Membership in", module, "module"),
       ylab = paste("OTU significance for ",parameter),
       caption=paste("Module membership vs. population significance\n"))
  ggsave(paste(module,"Module membership vs. population significance for", parameter,".png"),gp1,device="png")
  
  ##################### Then switch to PLS + VIP.
  # VIP stands for Variable Importance in the Projection.
  
  th_r2<-0.3 # We will only look at the PLS if the correlation is better than 0.3
  subnetwork<-transformed_abund_clean[,moduleGenes]
  subnetwork<- as.matrix(subnetwork)
  ncol(subnetwork)
  nrow(subnetwork)
  
  pls_result<-plsr(as.matrix(weight) ~ (subnetwork), validation="LOO",method="oscorespls")
#  print(pls_result)
  r2_vector<-R2(pls_result)
#  print(r2_vector)
  max<-0
  max_comp<--1
  for (j in 1:length(r2_vector$val)){
    if(r2_vector$val[j]>th_r2){
      if(r2_vector$val[j]>max){
        max<-r2_vector$val[j]
        max_comp<-r2_vector$comp[j]
      }
    }
  }
  
  source("http://mevik.net/work/software/VIP.R")
  
  print(paste(" the max r2 is ",max," corresponding to comp ",max_comp,sep=""))
  if(max==0){
    print ("No good correlation, we stop here")
    # Checking the normal correlation for everything
    df_normal<-as.data.frame(as.double(cor(transformed_abund_clean, Metadata$Group, use = "p")))
    colnames(df_normal)<-c("Normal_correlation")
    rownames(df_normal)<-colnames(transformed_abund_clean)
    write.csv(df_normal,file=paste("Normal_correlation_All_",parameter,".txt",sep=""))
  }else {
    print("Good correlation, we check the VIP !")
    # Checking the VIP
    output<-paste("VIP_values_",module,"_with_",parameter,".txt",sep="")
    vip_result<-VIP(pls_result)
    vip_components<-sort(vip_result[max_comp,],decreasing=TRUE)[1:100]
    for (i in 1:100){
      cat(paste("Rank ",i," we have ",names(vip_components[i])," with a VIP of ",vip_components[i],"\n",sep=""),file=output,append=TRUE)
    }
    weight_2 <- as.data.frame(weight[!is.na(weight)])
    rownames(weight_2) <-cleannamesweight
    df<-data.frame(x=weight_2,y=pls_result$validation$pred[,,max_comp])
    colnames(df)<-c("x","y")
    # Establish the correlation between predicted and measured
    cor.test(df$x,df$y)
    # Prep data for VIP - correl hive plot
    df_2<-as.data.frame(as.double(cor(subnetwork, weight, use = "p")))
    colnames(df_2)<-c("Normal_correlation")
    rownames(df_2)<-colnames(subnetwork)
    write.csv(df_2,file=paste("Normal_correlation_",module, "_",parameter,".txt",sep=""))
    df_2<-cbind(Name=row.names(df_2),df_2)
    table<-row.names(df_2)
    df_2<-cbind(df_2,Axis=as.integer(c(1)))
    df_3<-as.data.frame(as.double(vip_result[max_comp,]))
    rownames(df_3)<-colnames(subnetwork)
    colnames(df_3)<-c("VIP_correlation")
    write.csv(df_3,file=paste("VIP_correlation_",module, "_",parameter,".txt",sep=""))
    
    
    # function to extract p-value from regression
    lmp <- function (modelobject) {
      if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
      f <- summary(modelobject)$fstatistic
      p <- pf(f[1],f[2],f[3],lower.tail=F)
      attributes(p) <- NULL
      return(p)
    }
    
    regression <- lm(df$x~df$y)
    # regression <- lm(df$x[-52]~df$y[-52])
    summary(regression)$r.squared
    lmp(regression)
    
    
    g <- ggplot(data=df) + theme_bw() + geom_abline() + 
      #geom_smooth(aes(x=x,y=y),method=lm, se = FALSE) + 
      xlab("Measured") + ylab("Predicted") + ggtitle(paste("Comparison of",parameter,"measured vs predicted for module",module)) + 
      theme(plot.title = element_text(hjust = 0.5),axis.text=element_text(color="black",size=10),axis.title=element_text(color="black",size=15), axis.ticks=element_line(color="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
      geom_point(aes(x=x,y=y), size=2.5, shape=21, fill=module, color="black") +
      scale_x_continuous(breaks = pretty(df$x, n=5)) +
      scale_y_continuous(breaks = pretty(df$x, n=5)) +
      annotate("text",x=Inf,y=-Inf,hjust=1.1,vjust=-0.5,label= paste("r2 =",format(round(summary(regression)$r.squared, 2), nsmall = 2) , "\n", "p-value =", format(round(lmp(regression), 2), nsmall = 2)))
    
    # dev.print(file="test.png", device=png, width=800)
    ggsave(plot=g, filename=paste0("measured_vs_predicted_",module,"-vs-",parameter,".png"), width=6, height=6)
    # 
    # pdf(paste("measured_vs_predicted_",module,"-vs-",parameter,".pdf"))
    # 
    
    
    ###################################
    ###################################
    Normal <- read.csv(paste("Normal_correlation_",module, "_",parameter,".txt",sep=""))
    #Normal$Normal <- "Normal"
    #Normal$Normal <- Normal$Normal_correlation
    Normal$X <- NULL
    #Normal$Normal <- NULL
    VIP <- read.csv(paste("VIP_correlation_",module, "_",parameter,".txt",sep=""))
    dat <- cbind(VIP, Normal)
    colnames(dat)[2]<- "VIP"
    #let's sort our df based on the VIP correlation, so we can see which OTUs have the highest correlation
    
    dat <-arrange(dat,desc(VIP),X)
    #View(dat)
    nrow(dat)
    #then we extract the 1000 OTUs that have the highest correlation
    #dat <- dat[1:100,]
    big_contributors <- as_tibble(dat)
    #    gtdb <- read.csv("../GTDB_modified_2.csv") %>%
    #      select("user_genome","Classification","family")
    big_contributors <- left_join(big_contributors,taxonomy, by=c("X" = "MAG"))# %>%
#      rename("MAG" = "X") #modified on OSC to call dply.rename()
    
    write.csv(big_contributors, file = paste0(module, "_",parameter,"_big_contributors.csv"),quote=FALSE,row.names = FALSE)
    lengthOfDataFrame <- dim(dat)[1]
    lengthOfDataFrame
    #Plotting with ggplot2
    
    dat_10_VIP <- as.data.frame(dat$VIP)
    dat_10_VIP$X <- dat$X
    dat_10_VIP$zeroes <- 0
    dat_10_Normal <- as.data.frame(dat$Normal)
    dat_10_Normal$X <- dat$X
    dat_10_Normal$zeroes <- 0
    pdf(paste("VIPs",module,"-vs-",parameter,".pdf"))
    ggplot(dat_10_VIP, aes(x=dat_10_VIP$`dat$VIP`, y=dat_10_VIP$zeroes)) + geom_point() + geom_point(data = dat_10_Normal, mapping = aes(x=dat_10_Normal$zeroes, y=dat_10_Normal$`dat$Normal`, col=factor(X))) + geom_point(aes(col=X)) + theme_bw() + geom_curve(data= dat_10_Normal, aes(x = dat_10_Normal$zeroes, y = dat_10_Normal$`dat$Normal`, xend = dat_10_VIP$`dat$VIP`, yend = dat_10_VIP$zeroes, col=factor(X))) + theme(legend.position="none", plot.margin=unit(c(0.5,0.3,3,0.7),"cm"))
    ggplot(dat_10_VIP, aes(x=dat_10_VIP$`dat$VIP`, y=dat_10_VIP$zeroes)) + geom_point() + geom_point(data = dat_10_Normal, mapping = aes(x=dat_10_Normal$zeroes, y=dat_10_Normal$`dat$Normal`, col=factor(X))) + geom_point(aes(col=X)) + theme_bw() + geom_curve(data= dat_10_Normal, aes(x = dat_10_VIP$`dat$VIP`, y = dat_10_VIP$zeroes, xend = dat_10_Normal$zeroes, yend = dat_10_Normal$`dat$Normal`, col=factor(X))) + theme(legend.position="none", plot.margin=unit(c(0.5,0.3,3,0.7),"cm"))
    #    ggsave(plot=g1, filename=paste("VIPs",module,"-vs-",parameter,".pdf"), width=6, height=6)
    dev.off()
  }
  setwd(abs_path)  
}

building_network <- function(sft_power,data,sex,value,ColumnName,day_excluded,type="Fecal"){
  setwd(abs_path)
  
  # sft_power<- sft_power
  # sex<- sex
  # value <- value
  # ColumnName <-ColumnName
  # day_excluded<- day_excluded
  
  print(paste0(sex,"_Day",day_excluded,"Excluded_",value))
  dir1=paste0(sex,"_Day",day_excluded,"Excluded_",value)
  
  if (str_detect(type,"Cecum")) {
    dir1=paste0(sex,"_Day",day_excluded,"Excluded_",value,"_FC")
  }
  #  dir.create(dir1)
  print(dir1)
  setwd(dir1)
  
  # transformed_abund<- fread(paste0("hellinger_transformed_TM_",dir1,".csv"),check.names = FALSE, sep = ",",header = TRUE)
  # transformed_abund<- as.data.frame(transformed_abund)
  # rownames(transformed_abund) <- transformed_abund$V1
  # transformed_abund[1:2,1:2]
  # transformed_abund$V1 <- NULL
  
  print("Building Network")
  print(sft_power)
  transformed_abund <- data
  net = blockwiseModules(transformed_abund, power = sft_power,
                         corType="bicor", 
                         maxBlockSize = 1000,
                         networkType = "signed", 
                         TOMType = "signed", minModuleSize = 15,
                         reassignThreshold = 0, deepSplit = 3, 
                         mergeCutHeight = 0.1,
                         numericLabels = FALSE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE, replaceMissingAdjacencies = TRUE,
                         saveTOMFileBase = "MB", verbose = 3)
  
  print("Network done")
  # assign colors to the modules
  moduleColors = labels2colors(net$colors)
  moduleColors<-net$colors
  MAG_modules <- as.data.frame(moduleColors)
  #  rownames(MAG_modules) <- colnames(transformed_abund)
  
  #  write.csv(MAG_modules,file=paste("00_MAGS_module_membership.txt",sep=""))
  fwrite(MAG_modules,
         paste("00_MAGS_module_membership.txt",sep=""),
         row.names=T)
  # now that we have our co-occuring OTUs defined by specific modules, we want to identify the environmental features that 
  # best correlate with the variation in the module disribution
  # attach environmental factors, and match the rows of environmental factors with the selected rows of metadata
  

  Metadata <- dplyr::filter(Meta, grepl(type, Type))
  Metadata <- dplyr::filter(Metadata, !grepl('Repeat_', Sample))
  Metadata <- dplyr::filter(Metadata, grepl(sex, Sex_SurgeryDate))
  Metadata <- dplyr::filter(Metadata, !grepl(day_excluded, Sex_SurgeryDate))
  Metadata <- dplyr::filter(Metadata, grepl(value, get(ColumnName)))
  
  Metadata <- Metadata %>%
    select("Sample","Group","Reads","Surgery_date","Six_month_BMS_score","disease_status", "Timepoint_num","After0")
  Metadata<- as.data.frame(Metadata)
  # ##
  rownames(Metadata) <- Metadata$Sample
  Metadata$Sample <- NULL
  #Metadata <- Metadata$Group
  
  # calculate module eigengenes, the most representative number for each OTU abundace in a module
  MB_eig = moduleEigengenes(transformed_abund, net$colors)$eigengenes
  
  # calculate module eigengenes, the most representative number for each OTU abundace in a module
  #  MB_eig = moduleEigengenes(transformed_abund, moduleColors)$eigengenes
  # "orderMEs" reorders given eigenvectors such that similar ones (as measured by correlation) are next to each other
  pops = orderMEs(MB_eig)
  rownames(pops) <- rownames(transformed_abund)
  
  ### we correlate our modules with the metadata using the newly calculated eigengenes
  moduleTraitCor = cor(pops, Metadata, use = "p")
  nSamples = nrow(transformed_abund)
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  textMatrix = paste(signif(moduleTraitCor, 2),round(signif(moduleTraitPvalue, 1),2), sep = "\n")
  dim(textMatrix) = dim(moduleTraitCor)
  
  pdf("00_Module-trait relationships.pdf",width=20,height=14)
  par(mar=c(8,12,4,8))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(Metadata),
                 yLabels = names(pops),
                 ySymbols = substring(names(pops), 3),
                 colorLabels = FALSE,
                 colors = brewer.pal(11,"RdBu"),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 2,
                 cex.lab = 2,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  
  corrs<-as.data.frame(moduleTraitCor) %>%
    rownames_to_column("Module")%>%
    pivot_longer(cols = -Module,names_to = "Variable",values_to = "Correlation") %>%
    drop_na()
  
  pvalues<-as.data.frame(moduleTraitPvalue) %>%
    rownames_to_column("Module")%>%
    pivot_longer(cols = -Module,names_to = "Variable",values_to = "pvalue") %>%
    drop_na()
  
  corr_pvalue<- inner_join(corrs,pvalues,by=c("Module","Variable")) %>%
    mutate(Module = str_remove(Module, "ME"))
  
  corr_pvalue$File <- dir1 
  corr_pvalue <- corr_pvalue %>%
    relocate(File, .before = 1)
  
  corr_pvalue_significant<- corr_pvalue %>%
    filter(abs(Correlation)>=0.3 & pvalue<=0.05)
  
  write.csv(corr_pvalue,"00_module_trait_corr_pvalue.csv",row.names=FALSE,quote=FALSE)
  
  all_modules_vars<- ggplot(corr_pvalue, aes(Variable,Module,color=pvalue)) +
    geom_point(aes(size = Correlation)) +
    scale_color_gradient(low = "cyan", high = "maroon") +
    #    scale_size(range = c(5, 15)) +
    labs(title = "Module vs Trait", x= "Variable",y= "Module", size = "Correlation",color="pvalue")
  ggsave("00_Module_Trait_Relationship.pdf", plot = all_modules_vars, width = 16, height = 12, units = "in")
  
  sig_modules_vars<- ggplot(corr_pvalue_significant, aes(Variable,Module,color=pvalue)) +
    geom_point(aes(size = Correlation)) +
    scale_color_gradient(low = "cyan", high = "maroon") +
    #    scale_size(range = c(5, 15)) +
    labs(title = "Modules vs Variables (Significant)", x = "Variable",y="Module", size = "Correlation",color="pvalue")
  ggsave("00_Significant_Module_Trait_Relationship.pdf", plot = sig_modules_vars, width = 16, height = 12, units = "in")
  corr_pvalue_significant <- as.data.frame(corr_pvalue_significant)
  #  setwd(abs_path)
  corr_pvalue<- read.csv("00_module_trait_corr_pvalue.csv")
  corr_pvalue_significant<- corr_pvalue %>%
    filter(abs(Correlation)>=0.3 & pvalue<=0.05)
  print(unique(corr_pvalue_significant$"Module"))
  print(unique(corr_pvalue_significant$"Variable"))
  for (i in unique(corr_pvalue_significant$"Module")) {
    for (var in unique(corr_pvalue_significant$"Variable")){
      tryCatch({
        print(i)
        print(var)
        setwd(abs_path)
        setwd(dir1)
        VIP_fn(i,var,Metadata,pops,nSamples,moduleColors,transformed_abund,sex,value,ColumnName,day_excluded,type="Fecal")
        setwd(abs_path)
        setwd(dir1)
        
      })
    }
  }
  
  setwd(abs_path)
  
  #  returns <- list("Metadata" = Metadata, "pops" = pops, "nSamples" = nSamples,"moduleColors" = moduleColors)
  
  #  return(returns)
}


M0<- filtering_to_picking_threshold("M","0dpi","Timepoint","1")
M7<- filtering_to_picking_threshold("M","7dpi","Timepoint","1")
M21<- filtering_to_picking_threshold("M","21dpi","Timepoint","1")
M35<- filtering_to_picking_threshold("M","35dpi","Timepoint","1")
M63<- filtering_to_picking_threshold("M","63dpi","Timepoint","1")
M6m<- filtering_to_picking_threshold("M","6_months","Timepoint","1")
#M6m_FC<- filtering_to_picking_threshold("M","6_months","Timepoint","1","Fecal|Cecum")

# M_lam<- filtering_to_picking_threshold("M","1","Group","1")
# M_T10<- filtering_to_picking_threshold("M","2","Group","1")
# M_T4<- filtering_to_picking_threshold("M","3","Group","1")

# 
# F0<- filtering_to_picking_threshold("F","0dpi","Timepoint","4")
# F7<- filtering_to_picking_threshold("F","7dpi","Timepoint","4")
# F21<- filtering_to_picking_threshold("F","21dpi","Timepoint","4")
# F35<- filtering_to_picking_threshold("F","35dpi","Timepoint","4")
# F63<- filtering_to_picking_threshold("F","63dpi","Timepoint","4")
# F6m<- filtering_to_picking_threshold("F","6_months","Timepoint","4")
# F6m_FC<- filtering_to_picking_threshold("F","6_months","Timepoint","4","Fecal|Cecum")

# F_lam<- filtering_to_picking_threshold("F","1","Group","4")
# F_T10<- filtering_to_picking_threshold("F","2","Group","4")
# F_T4<- filtering_to_picking_threshold("F","3","Group","4")


M0_list <- building_network(10,M0,"M","0dpi","Timepoint","1")
M7_list <- building_network(12,M7,"M","7dpi","Timepoint","1")
M21_list <- building_network(12,M21,"M","21dpi","Timepoint","1")
M35_list <- building_network(14,M35,"M","35dpi","Timepoint","1")
M63_list <- building_network(16,M63,"M","63dpi","Timepoint","1")
M6m_list <- building_network(18,M6m,"M","6_months","Timepoint","1")
#M6m_FC_list <- building_network(16,M6m_FC,"M","6_months","Timepoint","1","Fecal|Cecum")

# M_lam_list<- building_network(14,"M","1","Group","1")
# M_T10_list<- building_network(14,"M","2","Group","1")
# M_T4_list<- building_network(9,"M","3","Group","1")

# F0_list <- building_network(10,F0,"F","0dpi","Timepoint","4")
# F7_list <- building_network(12,F7,"F","7dpi","Timepoint","4")
# F21_list <- building_network(14,F21,"F","21dpi","Timepoint","4")
# F35_list <- building_network(14,F35,"F","35dpi","Timepoint","4")
# F63_list <- building_network(10,F63,"F","63dpi","Timepoint","4")
# F6m_list <- building_network(14,F6m,"F","6_months","Timepoint","4")
# F6m_FC_list <- building_network(14,F6m_FC,"F","6_months","Timepoint","4","Fecal|Cecum")


# F_lam_list<- building_network(9,F_lam,"F","1","Group","4")
# F_T10_list<- building_network(14,F_T10,"F","2","Group","4")
# F_T4_list<- building_network(10,F_T4,"F","3","Group","4")
