#---------------------------------
# unroot() in ape
# This is used to unroot tree
# and remove the bootstrap values
#---------------------------------
library(ape)
library(hongR)
library(treeio)
#BiocManager::install("ggtree")
library(treeio)
library(ggtree)
library(tidytree)
library(phytools)
library(ggplot2)
library(readxl)
library(stringr)
#devtools::install_github("tidyverse/readxl")



yeast_species <- read_excel("data_for_tree/343taxa_speicies-name_clade-name_color-code.xlsx")
# one small tasks
ura1_hgt <- read.table("data_for_tree/ura1_HGT.txt", header=FALSE, stringsAsFactors = FALSE)
ura1_hgt$V2 <- getSingleReactionFormula(yeast_species$`Major clade`, yeast_species$old_speceis_names, ura1_hgt$V1)
table(ura1_hgt$V2)


# visualize the tree based on groups
tree <- read.tree("data_for_tree/tree_XR.txt")
group_file <- read.table("data_for_tree/HGT_XR.txt",header = T,row.names = 1)
groupInfo <- split(row.names(group_file), group_file$HGTornot)
tree <- groupOTU(tree, groupInfo)
ggtree(tree, layout="circular", branch.length = "none",aes(color=group)) + geom_tiplab2(size=1) + theme(legend.position = "bottom")
ggtree(tree, layout="circular", branch.length = "none",aes(color=group)) + theme(legend.position = "top")



# tree-example URA1 HGT analysis
tree <- read.tree("data_for_tree/URA1_aa.tre")
ggtree(tree, branch.length='none', layout='circular')
ggtree(tree, branch.length='none', layout='rectangular')
annotation <- read.csv(file = 'data_for_tree/URA1_phylogeny.csv', stringsAsFactors = FALSE)
print(colnames(annotation))
# for easy observation, update the tips name
original_tip <- tree[["tip.label"]]
new_tip <- getSingleReactionFormula(annotation$strain,annotation$accession,original_tip)
tree2 <- tree
tree2[["tip.label"]] <- new_tip

group_sign <- "genus"
groupInfo <- split(annotation$strain, annotation[,group_sign])
tree <- groupOTU(tree2, groupInfo,group_sign)
ggtree(tree, layout="circular", branch.length = "none",aes(color=genus)) + theme(legend.position = "top") +
  geom_tiplab2(size=1.4)

group_sign <- "superkingdom"
groupInfo <- split(annotation$strain, annotation[,group_sign])
tree <- groupOTU(tree2, groupInfo,group_sign)
ggtree(tree, layout="circular", branch.length = "none",aes(color=superkingdom)) + theme(legend.position = "top") +
  geom_tiplab2(size=1.4)



# map the traits (crabtree and heat) onto the species tree
tree <- read.tree("data_for_tree/332_2408OGs_time-calibrated_phylogeny_species-names_updated.newick")
trait <- read_excel("data_for_tree/genome_summary_332_yeasts_heat_Ethanol_updated_02_20.xlsx", 
                    sheet = "heat")

# unify the name
tips <- tree[["tip.label"]]
trait$`Species name` <- str_replace_all(trait$`Species name` , " ", "_")

species_name_check <- data.frame(tip_from_tree=tips, stringsAsFactors = FALSE)
species_name_check$trait <- getSingleReactionFormula(trait$heat_tolerance,trait$`Species name`,species_name_check$tip_from_tree)
# two species name need to be checked
# Lachancea_fantastica_nom_nud  (from tree) --> Lachancea fantastica nom. nud. (from excel)
# Wickerhamomyces_sp._YB_2243  (from tree) --> Wickerhamomyces sp. (from excel)
trait$`Species name`[which(trait$`Species name`=="Lachancea_fantastica_nom._nud.")] <- "Lachancea_fantastica_nom_nud"
trait$`Species name`[which(trait$`Species name`=="Wickerhamomyces_sp.")] <- "Wickerhamomyces_sp._YB_2243"
# only keep three types: Yes, No, Not_sure
print(unique(trait$heat_tolerance))
trait$heat_tolerance[which(trait$heat_tolerance=="40 not available")] <- "Not_sure"
trait$heat_tolerance[which(trait$heat_tolerance=="40 variable")] <- "Not_sure"

group_sign <- "heat_tolerance"
groupInfo <- split(trait$`Species name`, trait[,group_sign])
tree <- groupOTU(tree, groupInfo, group_sign)
ggtree(tree, layout="circular", branch.length = "none",aes(color=heat_tolerance)) + theme(legend.position = "top") +
  geom_tiplab2(size=1.4)

ggtree(tree, layout="circular", branch.length = "none",aes(color=heat_tolerance)) + 
    scale_color_manual(values=c("steelblue", "grey", "firebrick")) + 
    theme(legend.position = "top")


# manipulate the tree with data
tree_inf <- as_tibble(tree)