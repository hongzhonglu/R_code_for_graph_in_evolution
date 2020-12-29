library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(hongR)
library(readxl)

sce_gene_summary <- read_excel("data/sce_gene_summary.xlsx")
sce_gene_summary$sce_gene <- str_replace_all(sce_gene_summary$sce_gene, "\\[", "") %>%
  str_replace_all(.,"\\]","") %>%
  str_replace_all(.,"\\'","") %>%
  str_replace_all(.,"Saccharomyces_cerevisiae@","")

sce_gene_OG_map <- splitAndCombine(sce_gene_summary$sce_gene,sce_gene_summary$OrthologID,",")
sce_gene_OG_map$v1 <- str_trim(sce_gene_OG_map$v1)

gene_check <- read.table("data/yeastGENEnotpresent.txt", header = FALSE, stringsAsFactors = FALSE)

sce_gene_OG_map_check <- sce_gene_OG_map[sce_gene_OG_map$v1 %in% gene_check$V1,]

# connect the OG id with the panID
OG_panid <- read_tsv("data/representatives.tsv")

sce_gene_OG_map_check$panID <- getSingleReactionFormula(OG_panid$representative,OG_panid$ortho_id,sce_gene_OG_map_check$v2)
colnames(sce_gene_OG_map_check) <- c("sce_geneid", "OG_ID", "pan_ID")

write.table(sce_gene_OG_map_check, "result/sce_gene_OG_map_check.txt", row.names = FALSE, sep = "\t")


# analyze the conservattion of sce genes in other yeast species
ortholog <- read_tsv("data/ortholog_occurence_num_all.tsv")

ggplot(sce_gene_summary, aes(other_strain_number)) + stat_ecdf() +
  xlab('Occurrence number') + ylab('Percentage') +
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('')
# output 6 x 5
