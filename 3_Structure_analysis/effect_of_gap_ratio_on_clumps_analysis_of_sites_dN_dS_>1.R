# These scripts are used to produce map related to evolution project
# here only the result from homolog pdb files is used, so how is about the experimental pdb files?
# 2019.2.20
# Hongzhong Lu

#load package
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(readxl)


# input the data 
file_general <- "data/CLUMPS from pdb_homo for dN_dS_max_gap_ratio_"
parameter <- c("0.2", "0.3", "0.5", "0.6", "0.7", "0.8")
parameter_as_col <- paste("gap_ratio_", parameter, sep="")
unqiue_protein_filter_0.05 <- c()
for (i in parameter){
  print(i)
  input_dir <- paste(file_general,i, "/pdb_info.txt",sep = "")
  result_general <- read.table(input_dir, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  result_general$p_value <- as.numeric(result_general$p_value)
  result_general0 <- result_general %>% filter(.,p_value <=0.0505)
  unique_protein_num <- length(unique(result_general0$locus))
  unqiue_protein_filter_0.05 <- c(unqiue_protein_filter_0.05, unique_protein_num)
}
effects_ratio <- data.frame(gap_ratio = parameter, num_protein = unqiue_protein_filter_0.05)


ggplot(data=effects_ratio, aes(x=gap_ratio, y=num_protein)) +
  geom_bar(stat="identity", color="blue",fill="blue", alpha=0.3) +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  labs(x="Gap frequency", y="Number of proteins")+
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial"),
        legend.text = element_text(size=20, family="Arial"))
ggsave(out <- paste('result/','Effect_of_gap_on_number_of_protein_with_significant_clusters','.svg', sep = ""), width=5, height=5, dpi=300)
