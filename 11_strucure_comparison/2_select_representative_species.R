# Similarity analysis based on model structure, trait and evolution distance

library(ape)
library(hongR)
library(ggplot2)
library(readxl)
library(stringr)
library(readr)
library(Biostrings)
library(tidyverse)


library(ape)
library(hongR)
library(treeio)
#BiocManager::install("ggtree")
library(treeio)
library(ggtree)
library(tidytree)
library(phytools)
library(ggtreeExtra)
library(ggnewscale)





species_id <- read_excel("data_for_tree/343taxa_speicies-name_clade-name_color-code.xlsx")
quality <- read_excel("data_for_large_structure_comparison/quality_analysis_for_each_species.xlsx")
quality$clade <- getSingleReactionFormula(species_id$`Major clade`,species_id$old_speceis_names,quality$species)

species_id$ratio_hq <- getSingleReactionFormula(quality$ratio_hq,quality$species,species_id$old_speceis_names)
species_id$num_hq <- getSingleReactionFormula(quality$num_hq,quality$species,species_id$old_speceis_names)

species_id$ratio_hq <- as.numeric(species_id$ratio_hq)
species_id$num_hq <- as.numeric(species_id$num_hq)



# input the evolution distance
taxa_pairwise_dist <- read_excel("data_for_large_structure_comparison/taxa_pairwise_dist.xlsx")


taxa_pairwise_dist$ratio_hq <- getSingleReactionFormula(species_id$ratio_hq,species_id$speceis_names_fig2,taxa_pairwise_dist$species_list)
taxa_pairwise_dist$num_hq <- getSingleReactionFormula(species_id$num_hq,species_id$speceis_names_fig2,taxa_pairwise_dist$species_list)
write.table(taxa_pairwise_dist, "data_for_large_structure_comparison/taxa_pairwise_dist_update.txt")



# mapping the select species onto the tree
taxa_pairwise_dist_update <- read_excel("data_for_large_structure_comparison/taxa_pairwise_dist_update.xlsx")
taxa_pairwise_dist_update$select[which(is.na(taxa_pairwise_dist_update$select))] <- "No"



tree <- read.tree("data_for_tree/332_2408OGs_time-calibrated_phylogeny_species-names_updated.newick")

group_sign <- "select"
groupInfo <- split(taxa_pairwise_dist_update$species_list, taxa_pairwise_dist_update[,group_sign])
tree <- groupOTU(tree, groupInfo, group_sign)
ggtree(tree, layout="circular", branch.length = "none",aes(color=select)) + theme(legend.position = "top") +
  geom_tiplab2(size=1.4)

ggtree(tree, layout="circular", branch.length = "none",aes(color=select)) + theme(legend.position = "top")



# save for the following analysis
species_id$select <- getSingleReactionFormula(taxa_pairwise_dist_update$select,taxa_pairwise_dist_update$species_list,species_id$speceis_names_fig2)

write.table(species_id, "data_for_large_structure_comparison/species_id.txt", sep = "\t")







