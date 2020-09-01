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
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = c(0.8, 0.2))

# plot the dN/dS larger than 0.5
B3 <- B1[B1$dN_dS>0.5, ]
ggplot(B3, aes(x = dN_dS, fill = group)) + geom_histogram(alpha = 0.5, bins = 50) +
  xlim(0.5, 1.5) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "none")










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
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold")) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = "none")

G1 <- filter(gene_dn_ds_all_new, gene_dn_ds_all_new$group == "g1")
G2 <- filter(gene_dn_ds_all_new, gene_dn_ds_all_new$group == "g2")
G6 <- filter(gene_dn_ds_all_new, gene_dn_ds_all_new$group == "g6")
G7 <- filter(gene_dn_ds_all_new, gene_dn_ds_all_new$group == "g7")

t.test(G1$dN_dS, G7$dN_dS)
t.test(G2$dN_dS, G7$dN_dS)
t.test(G6$dN_dS, G7$dN_dS)

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
