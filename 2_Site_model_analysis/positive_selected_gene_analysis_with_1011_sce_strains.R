# These scripts are used to produce map related to evolution project
# 2019.2.20
# Hongzhong Lu

#load package
library(readr)
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(readxl)
library(hongR)

# function to standardize the snp number for each gene
calc_relative_snp <- function(geneName, snp_number, annotation = s288_sgd_annotation) {
  #gene_list <- geneGEM$geneNames
  #snp_list <- geneGEM$SNP_NUM
  gene_list <- geneName
  snp_list <- snp_number
  snp_list <- as.numeric(snp_list)
  protein_length <- getSingleReactionFormula(annotation$protein_length, annotation$systematic_name, gene_list)
  protein_length <- as.numeric(protein_length)
  snp_relative <- snp_list / protein_length
  return(snp_relative)
}


s288_sgd_annotation <- read_tsv('data/s288_genome.tsv')
fubar_02_08 <- read_csv("data/fubar_02_08.csv")
og_panid <- read_tsv("data/representatives.tsv")
fubar_detail <- left_join(fubar_02_08, og_panid, by=c("OG"="ortho_id"))
# input the dN/ds
gene_dn_ds_1011 <- read_csv("data/gene_dn_ds_1011.csv")
gene_dn_ds_1011$OG <- str_replace_all(gene_dn_ds_1011$OG,".out_yn00","")


# only choose sce
fubar_sce <- filter(fubar_detail, str_detect(fubar_detail$representative, "Saccharomyces_cerevisiae@"))
fubar_sce$representative <- str_replace_all(fubar_sce$representative , "Saccharomyces_cerevisiae@", "")

# get the snp information
all_gene_SNP <- read.table("data/all_gene_with_SNP_number.txt", header = TRUE, stringsAsFactors = FALSE)
all_gene_SNP <- filter(all_gene_SNP, nsSNP !='NA')
all_gene_SNP$nsSNP <- as.numeric(all_gene_SNP$nsSNP)
all_gene_SNP$nsSNP_relative <- calc_relative_snp(all_gene_SNP$locus_tag,all_gene_SNP$nsSNP_unique)
# connect the snp data with selected information
fubar_sce_with_snp <- left_join(fubar_sce, all_gene_SNP, by=c("representative"="locus_tag"))
fubar_sce_with_snp$dN_dS_1011 <- getSingleReactionFormula(gene_dn_ds_1011$dN_dS,gene_dn_ds_1011$OG,fubar_sce_with_snp$representative)

# classify the gene into three groups
fubar_sce_with_snp$group[fubar_sce_with_snp$positive_num==0] <- "a. No sites"
fubar_sce_with_snp$group[fubar_sce_with_snp$positive_num>=1 & fubar_sce_with_snp$positive_num <=2 ] <- "b. 1-2 sites"
fubar_sce_with_snp$group[fubar_sce_with_snp$positive_num>=3] <- "c. over 3 sites"


# plot1
fubar_sce_with_snp %>%
  ggplot(aes(x=group, y=nsSNP_relative, fill=group))+
  geom_boxplot(outlier.size=0) +
  xlab('Selected site number') + ylab( 'Unique nsSNP per codon')  +
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=12,face="bold", family="Arial"),
        axis.title=element_text(size=16,face="bold", family="Arial"),
        legend.text = element_text(size=20,face="bold", family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) #+
  #ggsave(out <- paste('result/',' choose strain based on growth of glycerol','.eps', sep = ""), width=5, height=5, dpi=300)
g1 <- fubar_sce_with_snp$nsSNP_relative[fubar_sce_with_snp$group=="a. No sites"]
g2 <- fubar_sce_with_snp$nsSNP_relative[fubar_sce_with_snp$group=="b. 1-2 sites"]
g3 <- fubar_sce_with_snp$nsSNP_relative[fubar_sce_with_snp$group=="c. over 3 sites"]

t.test(g1, g2)
t.test(g2, g3)

# plot2
fubar_sce_with_snp0 <- fubar_sce_with_snp %>% filter(., dN_dS_1011!="NA")
fubar_sce_with_snp0$dN_dS_1011 <- as.numeric(fubar_sce_with_snp0$dN_dS_1011)
fubar_sce_with_snp0 %>%
  ggplot(aes(x=group, y=dN_dS_1011, fill=group))+
  geom_boxplot(outlier.size=0) +
  xlab('Selected site number') + ylab( 'dN_dS')  +
  theme(legend.position = "none") + 
  theme(axis.text=element_text(size=12,face="bold", family="Arial"),
        axis.title=element_text(size=16,face="bold", family="Arial"),
        legend.text = element_text(size=20,face="bold", family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) #+
#ggsave(out <- paste('result/',' choose strain based on growth of glycerol','.eps', sep = ""), width=5, height=5, dpi=300)
g1 <- fubar_sce_with_snp0$nsSNP_relative[fubar_sce_with_snp0$group=="a. No sites"]
g2 <- fubar_sce_with_snp0$nsSNP_relative[fubar_sce_with_snp0$group=="b. 1-2 sites"]
g3 <- fubar_sce_with_snp0$nsSNP_relative[fubar_sce_with_snp0$group=="c. over 3 sites"]

t.test(g1, g2)
t.test(g2, g3)

