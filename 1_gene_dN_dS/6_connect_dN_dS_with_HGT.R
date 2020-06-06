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

# input the OG with HGT event
# OG_HGT <- read.table("data/HGT_OG.txt", stringsAsFactors = FALSE)
OG_with_HGT_cell <- read_csv("data/OG_with_HGT_cell.csv")

all_OG <- unique(c(OG_with_HGT_cell$OG))

# classify all OG based on whether they have HGT
gene_dn_ds_all_new$type <- NA
gene_dn_ds_all_new$type[gene_dn_ds_all_new$OG %in% all_OG] <- "HGT"

gene_dn_ds_all_new$type[is.na(gene_dn_ds_all_new$type)] <- "non_HGT"


# plot
ggplot(gene_dn_ds_all_new,aes(x=type, y=dN_dS, fill=type)) + geom_boxplot() +
  xlab('') + ylab('dN_dS') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=16, family="Arial"),
        axis.title=element_text(size=16,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) #+

g1 <- gene_dn_ds_all_new$dN_dS[gene_dn_ds_all_new$type=="HGT"]
g2 <- gene_dn_ds_all_new$dN_dS[gene_dn_ds_all_new$type=="non_HGT"]
t.test(g1,g2)