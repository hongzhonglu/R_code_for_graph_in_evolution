library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(hongR)
library(readxl)

###################################################################################################
# extract the OGs with sce
# Analysis one-build the conservation and dn_ds
ortholog <- read_tsv("data/ortholog_occurence_num_all.tsv")
ortholog0 <- filter(ortholog, with_sce == "with_sce")
gene_dn_ds_all <- read_csv("data/gene_dn_ds_03_02.csv") %>%
  filter(., !is.na(dN_dS)) %>%
  filter(., dN_dS < 10)

gene_dn_ds_all$OG <- str_replace_all(gene_dn_ds_all$OG, ".out_yn00", "")
gene_dn_ds_sce <- gene_dn_ds_all[gene_dn_ds_all$OG %in% ortholog0$ID, ]

# get the relate sce geneid based on the id mapping
sce_gene_summary <- read_tsv("data/Saccharomyces_cerevisiae.tsv")

# duplicates of gene across 1011 sce
duplicate_seq <- read_csv("data/duplicate_gene_analysis_1011_sce.csv")

# find the OG id based on the gene id
duplicate_seq$OG_ID <- getSingleReactionFormula(sce_gene_summary$ortholog_id, sce_gene_summary$gene_id, duplicate_seq$cluster)

# plot the density plot
ggplot(duplicate_seq, aes(unique_num)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "blue",  bins = 50)


#### average protein conservation score #######################################################################
conservation_score_sce_1011 <- read_csv("data/conservation_score_sce_1011.csv")
conservation_score_sce_1011$OG <- str_replace_all(conservation_score_sce_1011$OG,"_aa_aligned.fasta","")
conservation_score_sce_1011 <- conservation_score_sce_1011[!is.na(conservation_score_sce_1011$conservation_score),]
ggplot(conservation_score_sce_1011, aes(conservation_score)) +
  geom_histogram(aes(y = ..density..), alpha = 0.5, fill = "blue",  bins = 50) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20, face = "bold"))  +
  xlab("Average value of conservation score")


low_conservation <- conservation_score_sce_1011[conservation_score_sce_1011$conservation_score < 0.65, ]
low_conservation_gene <- unique(low_conservation$OG)
paste0(low_conservation_gene, collapse = ",")

high_conservation <- conservation_score_sce_1011[conservation_score_sce_1011$conservation_score > 0.844, ]
high_conservation_gene <- unique(high_conservation$OG)
paste0(high_conservation_gene, collapse = ",")


#################################################################################################################
# conservation score analysis for OGs from 343 species contains sce seq
conservation_score_sce <- read_csv("data/conservation_score_sce.csv")
conservation_score_sce <- filter(conservation_score_sce, !is.na(conservation_score))
S1 <- conservation_score_sce[,c("OG", "conservation_score")]
S2 <- conservation_score_sce_1011[,c("OG", "conservation_score")]

S1$group <- "343_species"
S2$group <- "1011_sce"
merge0 <- rbind.data.frame(S1, S2)
ggplot(merge0, aes(x = conservation_score, fill = group)) +
  geom_histogram(alpha=0.6, position="identity") +
  #geom_density(alpha = 0.5) +
  xlim(0.25, 0.9) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = c(0.2, 0.8)) +
  theme(axis.text=element_text(size=20,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial"),
        legend.text = element_text(size = 13, family = "Arial")) +  ggtitle('')
