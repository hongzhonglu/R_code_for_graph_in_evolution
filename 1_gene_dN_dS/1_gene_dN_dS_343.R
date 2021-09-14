library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(hongR)
library(readxl)
library(vioplot) # for violpot
library(caroline) # for violpot

###################################################################################################
# new result for all the filtered OGs
# gene_dn_ds_all <- read_csv("data/gene_dn_ds_macse.csv") %>% # this is early version, not used!
#  filter(., !is.na(dN_dS)) %>%
#  filter(., dN_dS < 10)
gene_dn_ds_all_new <- read_csv("data/gene_dn_ds_03_02.csv") %>%
  filter(., !is.na(dN_dS)) %>%
  filter(., dN_dS < 10)
# plot the column density plot
ggplot(gene_dn_ds_all_new, aes(dN_dS)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "blue", bins = 50) +
  xlim(0, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(colour = "black", alpha = 0.3, fill = "grey50") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  geom_vline(
    xintercept = 1, linetype = "dotted",
    color = "red", size = 1
  ) +
  geom_vline(
    xintercept = 3, linetype = "dotted",
    color = "red", size = 1
  )


###################################################################################################
# extract the OGs with sce
# Analysis one-build the conservation and dn_ds
ortholog <- read_tsv("data/ortholog_occurence_num_all.tsv")
ortholog0 <- filter(ortholog, with_sce == "with_sce")
gene_dn_ds_all_new$OG <- str_replace_all(gene_dn_ds_all_new$OG, ".out_yn00", "")
gene_dn_ds_sce <- gene_dn_ds_all_new[gene_dn_ds_all_new$OG %in% ortholog0$ID, ]

# This is early version. Here the paralog is not removed. It seems the paralog does not affect the followed analysis.
# compare the two kind of result
# gene_dn_ds_OG_with_sce <- read_csv("data/gene_dn_ds_OG_with_sce.csv")
# gene_dn_ds_OG_with_sce$OG <- str_replace_all(gene_dn_ds_OG_with_sce$OG, ".out_yn00", "")
# t.test(gene_dn_ds_sce$dN_dS, gene_dn_ds_OG_with_sce$dN_dS)
#      OGs_need_check <- setdiff(gene_dn_ds_sce$OG, gene_dn_ds_OG_with_sce$OG)
# OGs_need_check2 <- setdiff(gene_dn_ds_OG_with_sce$OG, gene_dn_ds_sce$OG)#this number is right!!!

# plot the column density plot
ggplot(gene_dn_ds_sce, aes(dN_dS)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "blue") +
  xlim(0, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(colour = "black", alpha = 0.3, fill = "grey50") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  geom_vline(
    xintercept = 1, linetype = "dotted",
    color = "red", size = 1
  ) +
  geom_vline(
    xintercept = 3, linetype = "dotted",
    color = "red", size = 1
  )


# try to compare all OGs and OGs with sce
B1 <- select(gene_dn_ds_all_new, OG, dN_dS)
B2 <- select(gene_dn_ds_sce, OG, dN_dS)
B1$group <- "All_orthologs"
B2$group <- "Orthologs_with_sce"
merge_dN_dS <- rbind.data.frame(B1, B2)
ggplot(merge_dN_dS, aes(x = dN_dS, fill = group)) + geom_histogram(alpha = 0.5, bins = 50) +
  xlim(0, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial") ) +
  theme(legend.position = c(0.75, 0.2)) +
  xlab("dN/dS") + ylab("Count") +
  theme(legend.text=element_text(size=16, family="Arial"))
ggsave(out <- paste('result/','dN_dS_for_all_OGs','.svg', sep = ""), width=8, height=6, dpi=600)

# plot the dN/dS larger than 0.5
B3 <- B1[B1$dN_dS>0.5, ]
ggplot(B3, aes(x = dN_dS, fill = group)) + geom_histogram(alpha = 0.5, bins = 50) +
  xlim(0.5, 1.5) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text=element_text(size=30, family="Arial"),
        axis.title=element_text(size=36, family="Arial") ) +
  xlab("dN/dS") + ylab("Count")+
  theme(legend.position = "none")
ggsave(out <- paste('result/','dN_dS_larger_than_one_distribution','.svg', sep = ""), width=5, height=5, dpi=600)




###################################################################################################
# dN_dS for smaller subset of three yeast species
gene_dn_ds_sce_paired_species <- read_csv("data/gene_dn_ds_all_3_species.csv") %>%
  filter(., !is.na(dN_dS)) %>%
  filter(., dN_dS < 10)
# plot the column density plot
ggplot(gene_dn_ds_sce_paired_species, aes(dN_dS)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "blue", bins = 50) +
  xlim(0, 1.5) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(colour = "black", alpha = 0.3, fill = "grey50") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  geom_vline(
    xintercept = 1, linetype = "dotted",
    color = "red", size = 1
  ) +
  geom_vline(
    xintercept = 3, linetype = "dotted",
    color = "red", size = 1
  )


###################################################################################################
# try to put the above result together
S1 <- select(gene_dn_ds_sce, OG, dN_dS)
S2 <- select(gene_dn_ds_sce_paired_species, OG, dN_dS)
S1$group <- "343_species"
S2$group <- "several_species"
merge_dN_dS <- rbind.data.frame(S1, S2)
ggplot(merge_dN_dS, aes(x = dN_dS, fill = group)) + geom_density(alpha = 0.5) +
  xlim(0, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = c(0.8, 0.2))

ggplot(merge_dN_dS, aes(x = dN_dS, fill = group)) + geom_histogram(alpha = 0.5, bins = 50) +
  xlim(0, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = c(0.8, 0.2))

ggplot(merge_dN_dS, aes(x = group, y = dN_dS, fill = group)) + geom_boxplot() +
  ylim(0, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = c(0.2, 0.8))

# here we can that the average dN_dS from 343 species is smaller than dN_dS from three species.
t.test(S1$dN_dS, S2$dN_dS)


###################################################################################################
# here we will try to classify the OG based on the species num
summary(gene_dn_ds_all_new)
gene_dn_ds_all_new$species_num <- getSingleReactionFormula(ortholog$species_num, ortholog$ID, gene_dn_ds_all_new$OG)
gene_dn_ds_all_new$protein_num <- getSingleReactionFormula(ortholog$protein_num, ortholog$ID, gene_dn_ds_all_new$OG)
gene_dn_ds_all_new$average_duplicate <- getSingleReactionFormula(ortholog$average_duplicate, ortholog$ID, gene_dn_ds_all_new$OG)

gene_dn_ds_all_new$species_num <- as.numeric(gene_dn_ds_all_new$species_num)
gene_dn_ds_all_new$protein_num <- as.numeric(gene_dn_ds_all_new$protein_num)
gene_dn_ds_all_new$average_duplicate <- as.numeric(gene_dn_ds_all_new$average_duplicate)

ggplot(gene_dn_ds_all_new, aes(species_num)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "blue") +
  xlim(5, 50)
# all the data is grouped into the followed groups
# grou1
g1 <- 7:10
g2 <- 11:50
g3 <- 51:100
g4 <- 101:150
g5 <- 151:200
g6 <- 201:250
g7 <- 251:343
gene_dn_ds_all_new$group <- NA
for (i in 1:nrow(gene_dn_ds_all_new)) {
  print(i)
  if (gene_dn_ds_all_new$species_num[i] %in% g1) {
    gene_dn_ds_all_new$group[i] <- "g1"
  } else if (gene_dn_ds_all_new$species_num[i] %in% g2) {
    gene_dn_ds_all_new$group[i] <- "g2"
  } else if (gene_dn_ds_all_new$species_num[i] %in% g3) {
    gene_dn_ds_all_new$group[i] <- "g3"
  } else if (gene_dn_ds_all_new$species_num[i] %in% g4) {
    gene_dn_ds_all_new$group[i] <- "g4"
  } else if (gene_dn_ds_all_new$species_num[i] %in% g5) {
    gene_dn_ds_all_new$group[i] <- "g5"
  } else if (gene_dn_ds_all_new$species_num[i] %in% g6) {
    gene_dn_ds_all_new$group[i] <- "g6"
  } else {
    gene_dn_ds_all_new$group[i] <- "g7"
  }
}

gene_dn_ds_all_new$group <- as.factor(gene_dn_ds_all_new$group)
ggplot(gene_dn_ds_all_new, aes(x = group, y = dN_dS, color = group)) +
  geom_boxplot() +
  ylim(0,1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = "none") +  
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial") ) +
  xlab("Orthologs with species number") + ylab("dN/dS")
ggsave(out <- paste('result/','Relation_between_species_number_and_dN_dS_distribution','.pdf', sep = ""), width=8, height=6, dpi=600)

# other versions which chould be put in the main figures
# mutiple color
ggplot(gene_dn_ds_all_new, aes(x = group, y = dN_dS, fill=group))  + 
  stat_boxplot(geom ='errorbar', width = 0.25) + 
  geom_boxplot() +
  ylim(0,1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = "none") +  
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial") ) +
  xlab("Orthologs with species number") + ylab("dN/dS")
ggsave(out <- paste('result/','Relation_between_species_number_and_dN_dS_distribution','.eps', sep = ""), width=8, height=6, dpi=600)

# single color
ggplot(gene_dn_ds_all_new, aes(x = group, y = dN_dS))  + 
  stat_boxplot(geom ='errorbar', width = 0.25) + 
  geom_boxplot(color="blue") +
  ylim(0,1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = "none") +  
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial") ) +
  xlab("Orthologs with species number") + ylab("dN/dS")
ggsave(out <- paste('result/','Relation_between_species_number_and_dN_dS_distribution','.svg', sep = ""), width=8, height=6, dpi=600)




G1 <- filter(gene_dn_ds_all_new, gene_dn_ds_all_new$group == "g1")
G2 <- filter(gene_dn_ds_all_new, gene_dn_ds_all_new$group == "g2")
G5 <- filter(gene_dn_ds_all_new, gene_dn_ds_all_new$group == "g5")
G6 <- filter(gene_dn_ds_all_new, gene_dn_ds_all_new$group == "g6")
G7 <- filter(gene_dn_ds_all_new, gene_dn_ds_all_new$group == "g7")

# t.test
t.test(G1$dN_dS, G7$dN_dS)
t.test(G2$dN_dS, G7$dN_dS)
t.test(G6$dN_dS, G7$dN_dS)
t.test(G5$dN_dS, G6$dN_dS)

# wilcon.test
wilcox.test(G6$dN_dS, G7$dN_dS, alternative = "two.sided")
wilcox.test(G5$dN_dS, G6$dN_dS, alternative = "two.sided")

# save the gene dn_ds_all_new for re-usage
write.table(gene_dn_ds_all_new, "result/gene_dn_ds_all_new.txt", row.names = FALSE, sep = "\t")


###################################################################################################
# dN_dS for 1011 sce sequence project
# The result is not good
gene_dn_ds_1011 <- read_csv("data/gene_dn_ds_1011.csv") %>% filter(., !is.na(dN_dS)) # %>%
#  filter(., dN_dS < 3 & dN_dS > 0)

# plot the column density plot
ggplot(gene_dn_ds_1011, aes(dN_dS)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "blue", bins = 50) +
  xlim(0, 3) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_density(colour = "black", alpha = 0.3, fill = "grey50") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  geom_vline(
    xintercept = 1, linetype = "dotted",
    color = "red", size = 1
  ) +
  geom_vline(
    xintercept = 3, linetype = "dotted",
    color = "red", size = 1
  )

# enrichment analysis for gene with dN/dS>1.5
# cell wall, membrane
gene_dn_ds_1011$OG <- str_replace_all(gene_dn_ds_1011$OG, ".out_yn00", "")

gene_dn_ds_filter1 <- gene_dn_ds_1011[gene_dn_ds_1011$dN_dS > 1.5, ]
paste0(gene_dn_ds_filter1$OG, collapse = ",")

# compare essential gene/non essential gene
# or core gene or non core gene
# input the gene classification of sce
sce_core_gene_list <- read_excel("data/panGene_for manual check.xlsx")
sce_essential_gene_list <- read_csv("data/sce_essential_gene_list.csv")
gene_dn_ds_1011$core_type <- getSingleReactionFormula(sce_core_gene_list$gene_type, sce_core_gene_list$gene_simple, gene_dn_ds_1011$OG)
gene_dn_ds_1011$essential_type <- getSingleReactionFormula(sce_essential_gene_list$`essentiality consensus`, sce_essential_gene_list$locus, gene_dn_ds_1011$OG)

# core gene or non core gene comparison
gene_dn_ds_1011_filter1 <- gene_dn_ds_1011[gene_dn_ds_1011$core_type != "NA", ]
gene_dn_ds_1011_filter1$core_type <- as.factor(gene_dn_ds_1011_filter1$core_type)
ggplot(gene_dn_ds_1011_filter1, aes(x = core_type, y = dN_dS, fill = core_type)) + geom_boxplot() +
  ylim(0, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = c(0.2, 0.8))

violins(
  x = gene_dn_ds_1011_filter1$dN_dS, # Variable of interest
  by = gene_dn_ds_1011_filter1$core_type, # Grouping factor
  connect = NULL, # Don't draw connecting lines (or 'strings')
  connectcol = NULL, # Nope
  las = 1, # Rotate x-axis labels
  drawRect = TRUE, # Add the boxplot
  rectCol = "black", # Makes the boxplot stand out more
  deciles = FALSE, # Don't plot deciles b/c they are too busy
  quantiles = FALSE, # Don't plot quants b/c they are too busy
  SD.or.SE = NULL, # Don't plot SD or SE b/c they are too busy
  CImed = FALSE, # Don't plot the 95% CI for the median
  col = rep("gray60", 3), # Fill color for each box
  ylab = "dN/dS", # y-axis label
  xlab = "Group of gene" # x-axis label
)


# essential gene or non essentail gene comparison
gene_dn_ds_1011_filter2 <- gene_dn_ds_1011[gene_dn_ds_1011$essential_type != "NA", ]
gene_dn_ds_1011_filter2$essential_type <- as.factor(gene_dn_ds_1011_filter2$essential_type)

ggplot(gene_dn_ds_1011_filter2, aes(x = essential_type, y = dN_dS, fill = essential_type)) + geom_boxplot() +
  ylim(0, 2) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = c(0.2, 0.8))


violins(
  x = gene_dn_ds_1011_filter2$dN_dS, # Variable of interest
  by = gene_dn_ds_1011_filter2$essential_type, # Grouping factor
  connect = NULL, # Don't draw connecting lines (or 'strings')
  connectcol = NULL, # Nope
  las = 1, # Rotate x-axis labels
  drawRect = TRUE, # Add the boxplot
  rectCol = "black", # Makes the boxplot stand out more
  deciles = FALSE, # Don't plot deciles b/c they are too busy
  quantiles = FALSE, # Don't plot quants b/c they are too busy
  SD.or.SE = NULL, # Don't plot SD or SE b/c they are too busy
  CImed = FALSE, # Don't plot the 95% CI for the median
  col = rep("gray60", 3), # Fill color for each box
  ylab = "dN/dS", # y-axis label
  xlab = "Group of gene" # x-axis label
)



##########################################################################
### analysis based on rxn info
# combine model information with dN/dS
gene_core_rxn <- read.table("data/model_info/corerxn_panid.txt", header = FALSE, stringsAsFactors = FALSE)
gene_core_rxn$V1 <- str_trim(gene_core_rxn$V1, side = "both")
gene_core_rxn$V1[!str_detect(gene_core_rxn$V1, "@")] <- paste("Saccharomyces_cerevisiae@", gene_core_rxn$V1[!str_detect(gene_core_rxn$V1, "@")], sep = "")
# find OG id
og_pan <- read_tsv("data/representatives.tsv")
gene_core_rxn$OG <- getSingleReactionFormula(og_pan$ortho_id,og_pan$representative,gene_core_rxn$V1)
print(length(intersect(gene_core_rxn$OG, gene_dn_ds_all_new$OG)))
gene_core_rxn$dN_dS <- getSingleReactionFormula(gene_dn_ds_all_new$dN_dS, gene_dn_ds_all_new$OG,gene_core_rxn$OG )
gene_core_rxn <- gene_core_rxn[!(gene_core_rxn$dN_dS=="NA"),]
gene_core_rxn$dN_dS <- as.numeric(gene_core_rxn$dN_dS)
gene_core_rxn$type <- "1.core metabolic"



gene_accessory_rxn <- read.table("data/model_info/accerxn_panid.txt", header = FALSE, stringsAsFactors = FALSE)
gene_accessory_rxn$V1 <- str_trim(gene_accessory_rxn$V1, side = "both")
gene_accessory_rxn$V1[!str_detect(gene_accessory_rxn$V1, "@")] <- paste("Saccharomyces_cerevisiae@", gene_accessory_rxn$V1[!str_detect(gene_accessory_rxn$V1, "@")], sep = "")
# find OG id
og_pan <- read_tsv("data/representatives.tsv")
gene_accessory_rxn$OG <- getSingleReactionFormula(og_pan$ortho_id,og_pan$representative,gene_accessory_rxn$V1)
gene_accessory_rxn$dN_dS <- getSingleReactionFormula(gene_dn_ds_all_new$dN_dS, gene_dn_ds_all_new$OG,gene_accessory_rxn$OG )
gene_accessory_rxn <- gene_accessory_rxn[!(gene_accessory_rxn$dN_dS=="NA"),]
gene_accessory_rxn$dN_dS <- as.numeric(gene_accessory_rxn$dN_dS)
gene_accessory_rxn$type <- "2.accessory metabolic"




gene_pan_rxn <- read.table("data/model_info/panid_panrxn.txt", header = FALSE, stringsAsFactors = FALSE)
gene_pan_rxn$V1 <- str_trim(gene_pan_rxn$V1, side = "both")
gene_pan_rxn$V1[!str_detect(gene_pan_rxn$V1, "@")] <- paste("Saccharomyces_cerevisiae@", gene_pan_rxn$V1[!str_detect(gene_pan_rxn$V1, "@")], sep = "")
# find OG id
og_pan <- read_tsv("data/representatives.tsv")
gene_pan_rxn$OG <- getSingleReactionFormula(og_pan$ortho_id,og_pan$representative,gene_pan_rxn$V1)
gene_pan_rxn$dN_dS <- getSingleReactionFormula(gene_dn_ds_all_new$dN_dS, gene_dn_ds_all_new$OG,gene_pan_rxn$OG )
gene_pan_rxn <- gene_pan_rxn[!(gene_pan_rxn$dN_dS=="NA"),]
gene_pan_rxn$dN_dS <- as.numeric(gene_pan_rxn$dN_dS)
gene_pan_rxn$type <- "3.pan metabolic"



merge_dN_dS2 <- rbind.data.frame(gene_core_rxn, gene_accessory_rxn, gene_pan_rxn)
ggplot(merge_dN_dS2, aes(x = type, y = dN_dS, fill = type)) + geom_boxplot() +
  ylim(0, 0.5) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size =14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

g1 <- gene_core_rxn$dN_dS
g2 <- gene_accessory_rxn$dN_dS
g3 <- gene_pan_rxn$dN_dS

wilcox.test(g1,g3, alternative = "two.sided")
wilcox.test(g1,g2, alternative = "two.sided")
wilcox.test(g2,g3, alternative = "two.sided")




B1 <- select(gene_dn_ds_all_new, OG, dN_dS)
M2 <- select(gene_pan_rxn, OG, dN_dS)
M3 <- select(gene_core_rxn, OG, dN_dS)
B1$group <- "All_OG"
M2$group <- "Metabolic_OG"
M3$group <- "Core_metabolic_OG"
merge_dN_dS <- rbind.data.frame(B1, M2, M3)
ggplot(merge_dN_dS, aes(x = dN_dS, fill = group)) + geom_density(alpha = 0.5) +
  xlim(0, 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial"),
        legend.text = element_text(size=16, family="Arial")) +
  theme(legend.position = c(0.75, 0.3)) +
  xlab("dN/dS") + ylab("Density") +
  theme(legend.text=element_text(size=14, family="Arial"))

write.csv(merge_dN_dS, "result/dataset_for_new_Fig.3A.csv")




# add the species number informattion
gene_core_rxn$species <- getSingleReactionFormula(ortholog$species_num,ortholog$ID,gene_core_rxn$OG)
gene_core_rxn$species <- as.numeric(gene_core_rxn$species)

gene_accessory_rxn$species <- getSingleReactionFormula(ortholog$species_num,ortholog$ID,gene_accessory_rxn$OG)
gene_accessory_rxn$species <- as.numeric(gene_accessory_rxn$species)

merge_dN_dS_species <- rbind.data.frame(gene_core_rxn, gene_accessory_rxn)


# grou1
g1 <- 0:100
g2 <- 100:250
g3 <- 250:300
g4 <- 300:343
merge_dN_dS_species$group <- NA
for (i in 1:nrow(merge_dN_dS_species)) {
  print(i)
  if (merge_dN_dS_species$species[i] %in% g1) {
    merge_dN_dS_species$group[i] <- "g1"
  } else if (merge_dN_dS_species$species[i] %in% g2) {
    merge_dN_dS_species$group[i] <- "g2"
  } else if (merge_dN_dS_species$species[i] %in% g3) {
    merge_dN_dS_species$group[i] <- "g3"
  } else{
    merge_dN_dS_species$group[i] <- "g4"
  } 
}

pd = position_dodge(width = 0.75)
ggplot(merge_dN_dS_species, aes(x=group, y=dN_dS, fill=type)) + 
  stat_boxplot(geom="errorbar", position=pd, width=0.2) +# add caps +
  geom_boxplot()+
  xlab('Species number') + ylab('dN/dS') +
  ylim(0,0.6)+
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = c(0.3,0.9)) +
  theme(axis.text=element_text(size=16, family="Arial"),
        axis.title=element_text(size=20,family="Arial"),
        legend.text = element_text(size=16, family="Arial")) +
  ggtitle('')

# statitical analysis
merge_dN_dS_species1 <- merge_dN_dS_species[merge_dN_dS_species$group=="g1", ]
wilcox.test(merge_dN_dS_species1$dN_dS[merge_dN_dS_species1$type=="1.core metabolic"], merge_dN_dS_species1$dN_dS[merge_dN_dS_species1$type=="2.accessory metabolic"], alternative = "two.sided")


merge_dN_dS_species2 <- merge_dN_dS_species[merge_dN_dS_species$group=="g2", ]
wilcox.test(merge_dN_dS_species2$dN_dS[merge_dN_dS_species2$type=="1.core metabolic"], merge_dN_dS_species2$dN_dS[merge_dN_dS_species2$type=="2.accessory metabolic"], alternative = "two.sided")


merge_dN_dS_species3 <- merge_dN_dS_species[merge_dN_dS_species$group=="g3", ]
wilcox.test(merge_dN_dS_species3$dN_dS[merge_dN_dS_species3$type=="1.core metabolic"], merge_dN_dS_species3$dN_dS[merge_dN_dS_species3$type=="2.accessory metabolic"], alternative = "two.sided")


merge_dN_dS_species4 <- merge_dN_dS_species[merge_dN_dS_species$group=="g4", ]
wilcox.test(merge_dN_dS_species4$dN_dS[merge_dN_dS_species4$type=="1.core metabolic"], merge_dN_dS_species4$dN_dS[merge_dN_dS_species4$type=="2.accessory metabolic"], alternative = "two.sided")

# compare g3 & g4
wilcox.test(merge_dN_dS_species4$dN_dS[merge_dN_dS_species4$type=="1.core metabolic"], merge_dN_dS_species3$dN_dS[merge_dN_dS_species3$type=="1.core metabolic"], alternative = "two.sided")
wilcox.test(merge_dN_dS_species4$dN_dS[merge_dN_dS_species4$type=="1.core metabolic"], merge_dN_dS_species3$dN_dS[merge_dN_dS_species3$type=="2.accessory metabolic"], alternative = "two.sided")



'
# the rxn and GPR mapping in the pan-GEM
library(readxl)
rxnMatrix_model <- read_excel("data/model_info/rxnMatrix_model.xlsx", 
                              col_names = FALSE)

colnames(rxnMatrix_model) <- c("rxn_num","rxnID","GPR")

rxnMatrix_model$GPR0 <- rxnMatrix_model$GPR %>% str_replace_all(.,"\\(", "") %>%
  str_replace_all(.,"\\)", "") %>% str_replace_all(.," and ", "&&") %>% str_replace_all(.," or ", "&&")
rxnMatrix_model01 <- splitAndCombine(rxnMatrix_model$GPR0, rxnMatrix_model$rxnID, sep0 = "&&")
colnames(rxnMatrix_model01) <- c("panID","rxnID")
rxnMatrix_model01$panID <- str_trim(rxnMatrix_model01$panID, side = "both")
rxnMatrix_model01$rxnID <- str_trim(rxnMatrix_model01$rxnID, side = "both")
rxnMatrix_model01$combine <- paste(rxnMatrix_model01$panID, rxnMatrix_model01$rxnID, sep = "&&")
rxnMatrix_model02 <- rxnMatrix_model01[!duplicated(rxnMatrix_model01$combine),]
og_rxn_num <- as.data.frame(table(rxnMatrix_model02$panID), stringsAsFactors = FALSE)'
