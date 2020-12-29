# The script is to analyze the site function type for unique mutation site related to the heat-tolerance types.
# 2020-8-30
# Hongzhong Lu

#load package
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(hongR)
library(tidyr)

site_type_unique_mutation <- read_csv("data/site_type_analysis_for_unique_mutation_analysis.csv")
site_type_unique_mutation$coordinate <- str_replace_all(site_type_unique_mutation$coordinate, "\\[", "") %>%
  str_replace_all(., "\\]", "")

site_type_unique_mutation0 <- splitAndCombine(site_type_unique_mutation$coordinate, site_type_unique_mutation$site_type, ",")
colnames(site_type_unique_mutation0 ) <- c("coordinate", "site_type")
site_type_unique_mutation1 <- site_type_unique_mutation0 %>% separate(site_type, c("gene", "none", "type", "others2"), sep = "@")
site_type_unique_mutation1$type[str_detect(site_type_unique_mutation1$type, "interface")] <- "interface"


# here based on the second independent calculation, we need to filter some genes which are not existing in the
# second calculation

# input the common gene
filter_gene <- read.table("result/commom_select_gene_for_heat_in_two_independent_calculation.txt", header=TRUE, sep = "\t", stringsAsFactors = FALSE)
site_type_unique_mutation1 <- site_type_unique_mutation1[site_type_unique_mutation1$gene %in% filter_gene$sce,]



# unify the site_type
print(unique(site_type_unique_mutation1$type))
# change  "alpha-helix" , "3-10-helix" into  "Helix"
# change  "extended strand" into "Beta strand"
site_type_unique_mutation1$type[str_detect(site_type_unique_mutation1$type, "extended strand")] <- "Beta strand"
site_type_unique_mutation1$type[str_detect(site_type_unique_mutation1$type, "alpha-helix")] <- "Helix"
site_type_unique_mutation1$type[str_detect(site_type_unique_mutation1$type, "3-10-helix")] <- "Helix"
# remove the duplicate row
site_type_unique_mutation1$combine <- paste(site_type_unique_mutation1$coordinate, site_type_unique_mutation1$gene, site_type_unique_mutation1$type, sep = "_")
site_type_unique_mutation2 <- site_type_unique_mutation1[!duplicated(site_type_unique_mutation1$combine),]

result_summary <- as.data.frame(table(site_type_unique_mutation2$type))
colnames(result_summary) <- c("Type","Count")


result_summary <- result_summary[order(result_summary$Count),]
result_summary$Type <-factor(result_summary$Type, levels=result_summary$Type)


ggplot(data=result_summary, aes(x=Type, y=Count)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=16, family="Arial"),
        axis.title=element_text(size=20, family="Arial"),
        legend.text = element_text(size=20, family="Arial"))
# output size 5 x 5 inches