library(tidyverse)

setwd("./M_Day1Excluded_35dpi")
###########################


table1<-read_csv("brown_Six_month_BMS_score_big_contributors.csv")
names(table1)[1] <- "MAG"

#combined_table<- rbind(table1,table2)
abund<- read_csv("hellinger_transformed_TM_M_Day1Excluded_35dpi.csv") #it is hellinger but just named RCLR
names(abund)[1] <- "Sample"
abund<- abund%>%
  column_to_rownames("Sample")


module_membership<-read_csv("00_MAGS_module_membership.txt")
adj<- adjacency(abund,power = 14,type="signed",corFnc = "bicor")

idx_brown<-match(as.array(table1$MAG),colnames(adj))
adjmat_brown<- as.data.frame(adj[idx_brown, idx_brown])
adjmat_brown[lower.tri(adjmat_brown)] <- NA
upper<-adjmat_brown%>%
  rownames_to_column("MAG") %>%
  pivot_longer(-MAG,names_to = "target",values_to="abund") %>%
  drop_na() %>%
  filter(abund<1) %>%
  inner_join(table1) %>%
  select(c(1:4))




write.csv(upper,"../Network_male_35dpi_brown_adj_upper_long.csv")

          