# Check the conservation of the positive selected genes occured in the branch site model related to the crabtree effect and heat-tolerance

# 2020.6.4
# Hongzhong Lu

#load package
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(hongR)


sce_gene_summary <- read_excel("data/sce_gene_summary.xlsx")
sce_gene_summary$sce_gene <- str_replace_all(sce_gene_summary$sce_gene,"\\[", "") %>% str_replace_all(.,"Saccharomyces_cerevisiae@","") %>%
  str_replace_all(.,"\\]", "") %>% str_replace_all(.,"'", "")


# input more parameters related to ortholog
ortholog_occurence <- read_tsv("data/ortholog_occurence_num_all.tsv")

sce_gene_summary0 <- left_join(sce_gene_summary,ortholog_occurence,by=c("OrthologID"="ID"))

# collect OG with the positive selection
plot(density(ortholog_occurence$species_num))
positive_select_gene <- read_excel("data/positive_select_gene.xlsx")
ortholog_occurence1 <- ortholog_occurence[ortholog_occurence$ID %in% positive_select_gene$OG,] %>% select(.,ID,species_num)
ortholog_occurence1$type <- "selected"
ortholog_occurence2 <- ortholog_occurence[!(ortholog_occurence$ID %in% positive_select_gene$OG),] %>% select(., ID, species_num)
ortholog_occurence2$type <- "reference"

select_gene_occurance <- rbind.data.frame(ortholog_occurence1, ortholog_occurence2)


# for better visualization, remove og with less 7 species
select_gene_occurance <- filter(select_gene_occurance, species_num >=7)

# plot
ggplot(select_gene_occurance, aes(x=type, y=species_num, fill=type)) + geom_boxplot() +
  xlab('') + ylab('Species_num') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=16, family="Arial"),
        axis.title=element_text(size=16,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) #+







