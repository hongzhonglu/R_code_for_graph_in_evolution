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




# analyse data at the clade levels
ggplot(quality, aes(x=clade, y=ratio_hq)) + 
  geom_boxplot() +
  xlab("Clade") + 
  ylab("Ratio of high-quality structures") + 
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=12, family="Arial"),
        legend.text = element_text(size=12, family="Arial")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))




ggplot(quality, aes(x=clade, y=num_hq)) + 
  geom_boxplot() +
  xlab("Clade") + 
  ylab("Number of high-quality structures") + 
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=12, family="Arial"),
        legend.text = element_text(size=12, family="Arial")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))



# input the tree
yeast_species <- read_excel("data_for_tree/343taxa_speicies-name_clade-name_color-code.xlsx")
tree <- read.tree("data_for_tree/332_2408OGs_time-calibrated_phylogeny_species-names_updated.newick")
# manipulate the tree with data
tree_inf <- as_tibble(tree)
tree_inf$clade <- getSingleReactionFormula(yeast_species$`Major clade`,yeast_species$speceis_names_fig2, tree_inf$label)
tree_inf$Family <- getSingleReactionFormula(yeast_species$Family, yeast_species$speceis_names_fig2, tree_inf$label)
tree_inf$Genus <- getSingleReactionFormula(yeast_species$Genus, yeast_species$speceis_names_fig2, tree_inf$label)



tips <- tree_inf$label
tips <- tips[1:332]
species_id <- species_id[species_id$speceis_names_fig2 %in% tips, ]


p<-ggtree(tree,layout = "circular",size=0.1)
p

p1<-p+
  geom_tiplab(align = T,size=2)
p1






ggtree(tree, branch.length='none', layout='circular') +
  geom_text(aes(label=node), hjust=-.3, size=1.5) +
  geom_tiplab(aes(angle=angle), size=1.5) + 
  geom_hilight(node=656, fill="steelblue", alpha=.6)  


p<-ggtree(tree,layout = "circular",size=0.25)
p

p1<-p##+
  ##geom_tiplab(align = T,size=1)
p1

p3 <- p1 + 
  new_scale_fill() +  ##ggplot2只允许设置一个scale，这个函数可以添加新的scale
  geom_fruit(
    data=species_id,
    geom=geom_bar,  #这里很重要，包含树tip标签的列在映射中应为 Y
    mapping=aes(y=speceis_names_fig2, x=ratio_hq),  
    color ="white",
    fill ="red",
    pwidth=0.25,
    stat="identity",
    offset = 0.1,
    axis.params=list(
      axis="x", # add axis text of the layer.
      text.angle=-45, # the text size of axis.
      hjust=0  # adjust the horizontal position of text of axis.
    ),
    grid.params=list()
  ) 
p3
p4 <- p3 + 
  new_scale_fill() +  ##ggplot2只允许设置一个scale，这个函数可以添加新的scale
  geom_fruit(
    data=species_id,
    geom=geom_bar,  #这里很重要，包含树tip标签的列在映射中应为 Y
    mapping=aes(y=speceis_names_fig2, x=num_hq),  
    color ="white",
    fill ="blue",
    pwidth=0.25,
    stat="identity",
    offset = 0.025,
    axis.params=list(
      axis="x", # add axis text of the layer.
      text.angle=-45, # the text size of axis.
      hjust=0  # adjust the horizontal position of text of axis.
    ),
    grid.params=list()
  ) 
p4
