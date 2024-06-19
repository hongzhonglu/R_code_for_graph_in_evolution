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


#####
# a method to label tree
#####
tr <- "((Pan_paniscus,Pan_troglodytes),((Homo_sapiens,Homo_erectus),Homo_abilis));"
tr <- read.tree(text = tr)
tr <- makeNodeLabel(tr, "u", nodeList = list(Pan = "Pan", Homo = "Homo"))
plot(tr, show.node.label = TRUE)
### does not erase the previous node labels:
tr <- makeNodeLabel(tr, "u", nodeList = list(Hominid = c("Pan","Homo")))
plot(tr, show.node.label = TRUE)
### the two previous commands could be combined:
L <- list(Pan = "Pan", Homo = "Homo", Hominid = c("Pan","Homo"))
tr <- makeNodeLabel(tr, "u", nodeList = L)
### combining different methods:
tr <- makeNodeLabel(tr, c("n", "u"), prefix = "#", nodeList = list(Hominid = c("Pan","Homo")))
plot(tr, show.node.label = TRUE)
#write.tree(tr, file = "/Users/luho/Documents/pan_genome/test_for_label/test2.tre")
# node lables
tr <- "((Pan_paniscus,Pan_troglodytes),((Homo_sapiens,Homo_erectus),Homo_abilis));"
tr <- read.tree(text = tr)
plot(tr, show.node.label = TRUE)
nodelabels()
all_tree <-subtrees(tr)










######## more test
#step 1 read data
tr <- "(((Hylobates_EDN , (Orang_EDN , (Gorilla_EDN , (Chimp_EDN , Human_EDN )))), (Macaq_EDN , (Cercopith_EDN , (Macaq2_EDN , Papio_EDN )))), (Orang_ECP, ((Macaq_ECP, Macaq2_ECP), (Goril_ECP, Chimp_ECP, Human_ECP))));"
interest_set0 <- c("Macaq_ECP", "Macaq2_ECP", "Goril_ECP", "Chimp_ECP", "Human_ECP")

#step2 add special label for target species
interest_set <- paste("branch@", interest_set0)

tr1 <- tr
for(i in 1:length(interest_set)){
  print(i)
  tr1 <- gsub(interest_set0[i], interest_set[i], tr1)
}
#add the node
tr <- read.tree(text = tr1)
plot(tr, show.node.label = TRUE)
tr <- makeNodeLabel(tr, "u", nodeList = list(BRANCH = "branch@"))
plot(tr, show.node.label = TRUE)
# replace the label
tip_inf <- data.frame(label = tr$tip.label, stringsAsFactors = FALSE)
tip_inf$label2 <- gsub("branch@", "", tip_inf$label)
## rename_taxa use 1st column as key and 2nd column as value by default 
tr2 <- rename_taxa(tr, tip_inf, label, label2)
plot(tr2, show.node.label = TRUE)
node_inf <- tr2[["node.label"]]
node_inf2 <- gsub("BRANCH", " #1", node_inf)
tr2[["node.label"]] <- node_inf2
plot(tr2, show.node.label = TRUE)
#save the node
#write.tree(tr2, file = "/Users/luho/Documents/pan_genome/test_for_label/test2.tre")








### tree visualization
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)

ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()
ggtree(tree, color="firebrick", size=2, linetype="dotted")
ggtree(tree, branch.length="none")

set.seed(2017-02-16)
tree <- rtree(332)
ggtree(tree)
ggtree(tree, layout="slanted") 
ggtree(tree, layout="circular")
ggtree(tree, layout="fan", open.angle=120)
ggtree(tree, layout="equal_angle")
ggtree(tree, layout="daylight")
ggtree(tree, branch.length='none')
ggtree(tree, branch.length='none', layout='circular')
ggtree(tree, layout="daylight", branch.length = 'none')



# view clade
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)
p <- ggtree(tree) + geom_tiplab()
viewClade(p, MRCA(p, "I", "L"))


p2 <- p %>% collapse(node=21) + 
  geom_point2(aes(subset=(node==21)), shape=21, size=5, fill='green')
p2 <- collapse(p2, node=23) + 
  geom_point2(aes(subset=(node==23)), shape=23, size=5, fill='red')
print(p2)
expand(p2, node=23) %>% expand(node=21)


# Viewing Selected Clade
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)
p <- ggtree(tree) + geom_tiplab()
p
viewClade(p, MRCA(p, tip=c("I", "L")))


# subtree
data(chiroptera)
tr <- drop.tip(chiroptera, 16:921, subtree = TRUE)
plot(tr, font = c(rep(3, 15), rep(2, 3)), cex = 0.8,
     no.margin = TRUE)

data(bird.families)
zoom(bird.families, 1:15, col = "grey", no.margin = TRUE,
     subtree = TRUE)

zoom(bird.families, list(1:15, 38:48), col = rep("grey", 2),
     no.margin = TRUE, font = 1, subtree = TRUE)

nodelabels() # add node numbers
tiplabels()



## ggtree subset tree by the tip label and internal node number
beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)

# by the tip label
tree2 = tree_subset(beast_tree, "A/Swine/HK/168/2012", levels_back=4)  
ggtree(beast_tree) +
  geom_tiplab() + theme_tree2() 
ggtree(tree2, aes(color=group)) +
  geom_tiplab() +  xlim(0, 4) + theme_tree2() 

# by the clade
clade <- tree_subset(beast_tree, node=121, levels_back=0)
clade2 <- tree_subset(beast_tree, node=121, levels_back=2)
ggtree(clade) + geom_tiplab() + xlim(0, 5)
ggtree(clade2, aes(color=group)) + geom_tiplab() + 
  xlim(0, 8)

data(chiroptera)
nodes <- grep("Plecotus", chiroptera$tip.label)
chiroptera <- groupOTU(chiroptera, nodes)
clade <- MRCA(chiroptera, nodes)
x <- tree_subset(chiroptera, clade, levels_back = 0)
ggtree(chiroptera, aes(colour = group)) + 
  theme(legend.position = "none")
ggtree(x) + geom_tiplab() + xlim(0, 5)




#Simulating, plotting, extracting clades, & dropping tips
## (I'm going to first set the seed for repeatability)
set.seed(1)
## simulate a birth-death tree using phytools
tree <- pbtree(b = 1, d = 0.2, n = 40)
## stopping criterion is 40 extant species, in this case
plotTree(tree, setEnv = TRUE)
nodelabels()

## ok, now extract the clade descended from node #62
tt62 <- extract.clade(tree, 62)
plotTree(tt62)

## now drop 10 tips from the tree (I'm going to pick them at random)
dtips <- sample(tree$tip.label, 10)
dt <- drop.tip(tree, dtips)
plotTree(dt)
## we could also, say, drop all tips that go extinct before the present
## this is a fun way, but not the only way to do this:
et <- fancyTree(tree, type = "droptip", tip = getExtinct(tree), cex = 0.7)
print(et)
## rotating nodes
plotTree(tree, node.numbers = T)
## first, rotate about node #50
rt.50 <- rotate(tree, 50)
plotTree(rt.50)










# use the real data
yeast_species <- read_excel("data_for_tree/343taxa_speicies-name_clade-name_color-code.xlsx")
tree <- read.tree("data_for_tree/332_2408OGs_time-calibrated_phylogeny_species-names_updated.newick")
# manipulate the tree with data
tree_inf <- as_tibble(tree)
tree_inf$clade <- getSingleReactionFormula(yeast_species$`Major clade`,yeast_species$speceis_names_fig2, tree_inf$label)
tree_inf$Family <- getSingleReactionFormula(yeast_species$Family, yeast_species$speceis_names_fig2, tree_inf$label)
tree_inf$Genus <- getSingleReactionFormula(yeast_species$Genus, yeast_species$speceis_names_fig2, tree_inf$label)


plotTree(tree, node.numbers = F)
tt62 <- extract.clade(tree, 341)

plotTree(tt62)

# get pairwise taxa-taxa distance matrix
d=cophenetic(tree)
# print distance matrix d to text file
write.table(d, file = 'taxa_pairwise_dist.txt',sep = '\t', quote = FALSE, col.names=NA)



zoom(tree, 324:332, col = "grey", no.margin = TRUE,
     subtree = TRUE)

zoom(tree, 186:192, col = "grey", no.margin = TRUE,
     subtree = TRUE)

ggtree(tree, layout="circular")



## highlight specific clade
genus0 <- "Candida"
tree_inf0 <- filter(tree_inf, Genus==genus0)
node0 <- unique(tree_inf0$parent)[1:7]

ggtree(tree) + geom_tiplab()
ggplot(tree, aes(x, y)) + geom_tree() + theme_tree() + geom_text(aes(label=node), hjust=-.3) +geom_tiplab()
ggtree(tree, branch.length='none', layout='circular') +
  geom_text(aes(label=node), hjust=-.3, size=1.5) +
  geom_tiplab(aes(angle=angle), size=1.5) + 
  geom_hilight(node=656, fill="steelblue", alpha=.6)  


# one small tasks
ura1_hgt <- read.table("data_for_tree/ura1_HGT.txt", header=FALSE, stringsAsFactors = FALSE)
ura1_hgt$V2 <- getSingleReactionFormula(yeast_species$`Major clade`, yeast_species$old_speceis_names, ura1_hgt$V1)
table(ura1_hgt$V2)


# visualize the tree based on groups
tree <- read.tree("data_for_tree/tree_XR.txt")
group_file <- read.table("data_for_tree/HGT_XR.txt",header = T,row.names = 1)
groupInfo <- split(row.names(group_file), group_file$HGTornot)
tree <- groupOTU(tree, groupInfo)
ggtree(tree, layout="rectangular", ladderize = FALSE, branch.length = "none",aes(color=group)) + geom_tiplab2(size=1)
ggtree(tree, layout="circular", branch.length = "none",aes(color=group)) + geom_tiplab2(size=1) + theme(legend.position = "bottom")
ggtree(tree, layout="circular", branch.length = "none",aes(color=group)) + theme(legend.position = "top")


# Grouping Taxa
data(iris)
rn <- paste0(iris[,5], "_", 1:150)
rownames(iris) <- rn
d_iris <- dist(iris[,-5], method="man")

tree_iris <- ape::bionj(d_iris)
grp <- list(setosa     = rn[1:50],
            versicolor = rn[51:100],
            virginica  = rn[101:150])
# procedure 1
p_iris <- ggtree(tree_iris, layout = 'circular', branch.length='none')
groupOTU(p_iris, grp, 'Species') + aes(color=Species) +
  theme(legend.position="right")
# procedure 2
tree_iris <- groupOTU(tree_iris, grp, "Species")
ggtree(tree_iris, aes(color=Species), layout = 'circular', branch.length = 'none') + 
  theme(legend.position="right")



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

ggtree(tree, layout="circular", branch.length = "none",aes(color=heat_tolerance)) + theme(legend.position = "top")






# manipulate the tree with data
tree_inf <- as_tibble(tree)


