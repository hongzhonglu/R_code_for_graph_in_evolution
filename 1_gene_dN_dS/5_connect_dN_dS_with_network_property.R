library(readxl)
library(fitdistrplus)
library(hongR)
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)

###################################################################################################
# new result for all the filtered OGs
gene_dn_ds_all_new <- read_csv("data/gene_dn_ds_03_02.csv") %>%
  filter(., !is.na(dN_dS)) %>%
  filter(., dN_dS < 10)

gene_dn_ds_all_new$OG <- str_replace_all(gene_dn_ds_all_new$OG, ".out_yn00", "")


# input the id mapping
OG_pan <- read_tsv("data/representatives.tsv")

# input the panID and its connected rxn number
panID_linked_Unique_rxn <- read_excel("data/model_info/panID_linked_Unique_rxn.xlsx")
colnames(panID_linked_Unique_rxn) <- c("panID", "rxn_number")

for (i in 1:nrow(panID_linked_Unique_rxn)){
  if(!str_detect(panID_linked_Unique_rxn$panID[i], "@")){
    panID_linked_Unique_rxn$panID[i] <- paste("Saccharomyces_cerevisiae@", panID_linked_Unique_rxn$panID[i], sep = "")
  }
}


# get the ortholog id based on the panID
panID_linked_Unique_rxn$OG <- getSingleReactionFormula(OG_pan$ortho_id, OG_pan$representative, panID_linked_Unique_rxn$panID)

panID_linked_Unique_rxn$dN_dS <- getSingleReactionFormula(gene_dn_ds_all_new$dN_dS, gene_dn_ds_all_new$OG, panID_linked_Unique_rxn$OG)
panID_linked_Unique_rxn$dN_dS <- as.numeric(panID_linked_Unique_rxn$dN_dS)
rxn_dN_dS <- panID_linked_Unique_rxn[!is.na(panID_linked_Unique_rxn$dN_dS),]


# change the scatter into the column plot
rxn_dN_dS$rxn_type <- NA
rxn_dN_dS$rxn_type[rxn_dN_dS$rxn_number==1] <- "A.one_rxn"
rxn_dN_dS$rxn_type[rxn_dN_dS$rxn_number==2] <- "B.two_rxns"
rxn_dN_dS$rxn_type[rxn_dN_dS$rxn_number==3] <- "C.three_rxns"
rxn_dN_dS$rxn_type[rxn_dN_dS$rxn_number>=4] <- "D.over_four_rxns"
# plot
ggplot(rxn_dN_dS,aes(x=rxn_type, y=dN_dS, fill=rxn_type)) + geom_boxplot() +
  xlab('') + ylab('dN_dS') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=16, family="Arial"),
        axis.title=element_text(size=16,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) #+

g1 <- rxn_dN_dS$dN_dS[rxn_dN_dS$rxn_number >= 4]
g2 <- rxn_dN_dS$dN_dS[rxn_dN_dS$rxn_number == 1]
t.test(g1,g2)


# analyse the gene withour dN/dS calculation
rxn_no_dN_dS <- panID_linked_Unique_rxn[is.na(panID_linked_Unique_rxn$dN_dS),]
which(rxn_no_dN_dS$OG %in% gene_dn_ds_all_new$OG)
# get the species number information for these OG
OG_all <- read_tsv("data/ortholog_occurence_num_all.tsv")

OG_all_check1 <- OG_all[OG_all$ID %in% rxn_no_dN_dS$OG,]
plot(density(OG_all_check1$species_num))
OG_all_check2 <- OG_all_check1[OG_all_check1$species_num <7,]
plot(density(OG_all_check2$species_num))

write_csv(OG_all_check2, "result/OG_in_model_need_check.csv")


# calculate characteristics of network
# input OG, rxn matrix
# the connectivity is defined by the metabolites in the model according to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1779518/
uniquerxn_OG_matrix <- read_excel("data/model_info/uniquerxn_OG_matrix.xlsx")