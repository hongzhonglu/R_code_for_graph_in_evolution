# note------------
# here we try to analyze the functions of genes with positive selected sites based on pathways.
# 2020-06-02
# Hongzhong Lu


library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(hongR)
library(readxl)
library (vioplot) # for violpot
library(caroline) # for violpot


# now the result only contains the FUBAR analysis result. In the future, the fel analysis results will be included.
###################################################################################################
# positive selected genes initial analysis
# read all OG information
gene_dn_ds_all_new <- read.table("result/gene_dn_ds_all_new.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# positive selected gene analysis across 343 yeast species
fubar_02_08 <- read_csv("data/fubar_02_08.csv")
select_gene_dn_ds <- merge.data.frame(gene_dn_ds_all_new, fubar_02_08, by.x = "OG", by.y="OG")
select_gene_dn_ds$select_num <- NA
select_gene_dn_ds$select_num[select_gene_dn_ds$positive_num >= 3] <- "c over_3_sites"
select_gene_dn_ds$select_num[select_gene_dn_ds$positive_num < 3 & select_gene_dn_ds$positive_num >= 1] <- "b 1-2_sites"
select_gene_dn_ds$select_num[select_gene_dn_ds$positive_num < 1] <- "a no_select_sites"
select_gene_dn_ds$select_num <- as.factor(select_gene_dn_ds$select_num)
# plot
ggplot(select_gene_dn_ds,aes(x=select_num, y=dN_dS, fill=select_num)) + geom_boxplot() +
  ylim(0,1) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = c(0.2, 0.8)) +
  theme(axis.text=element_text(size=12,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial"),
        legend.text = element_text(size = 13, family = "Arial")) +  ggtitle('') +
        theme(legend.position = "none")

SG1 <-  filter(select_gene_dn_ds, select_num == "c over_3_sites")
SG2 <-  filter(select_gene_dn_ds, select_num == "b 1-2_sites")
SG3 <-  filter(select_gene_dn_ds, select_num == "a no_select_sites")


t.test(SG1$dN_dS,SG3$dN_dS)
t.test(SG1$dN_dS,SG2$dN_dS)

# more analysis about the positive selected number
select_gene <- select_gene_dn_ds[select_gene_dn_ds$positive_num >= 1, ]
select_gene$group[select_gene$positive_num >= 3] <- "over_3_select_sites"
select_gene$group[select_gene$positive_num == 2] <- "2_select_sites"
select_gene$group[select_gene$positive_num <= 1] <- "1_select_sites"
# plot
ggplot(select_gene, aes(group)) +
  geom_bar(fill = "#FF6666") +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(legend.position = c(0.2, 0.8)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(axis.text=element_text(size=12,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial"),
        legend.text = element_text(size = 13, family = "Arial")) +  ggtitle('')


###############################################################################################################
# classification of OGs with selected sites based on KO or function annotation.
og_panID_mapping <- read_tsv("data/representatives.tsv")
select_gene <- left_join(select_gene, og_panID_mapping, by = c("OG" = "ortho_id"))

# get the KO id
panID_KO <- read.table("kegg/ko_annotation_of_all_OG.txt", header = TRUE, sep = "\t")
select_gene$KO <- getSingleReactionFormula(panID_KO$ko, panID_KO$query, select_gene$representative)

# get the pathway
KO_pathway <- read.table("kegg/ko_pathway.txt", sep = "\t", stringsAsFactors = FALSE)
pathway <- read.table("kegg/pathway_list_kegg.txt", sep = "\t", stringsAsFactors = FALSE)
KO_pathway <- filter(KO_pathway,str_detect(KO_pathway$V1, "map"))
KO_pathway$V2 <- str_replace_all(KO_pathway$V2,"ko:","")
select_gene$pathwayID <- getSingleReactionFormula(KO_pathway$V1,KO_pathway$V2, select_gene$KO)
select_gene$pathway <- getSingleReactionFormula(pathway$V2, pathway$V1, select_gene$pathway)
select_gene$pathway[select_gene$pathway=="NA"] <- NA
select_gene$pathway <- str_replace_all(select_gene$pathway, " - fly", "") %>% str_replace_all(.," - other", "") %>%
  str_replace_all(.," - yeast", "") %>% str_replace_all(.," - Caulobacter", "") %>% str_replace_all(.," - animal", "")

select_gene_filter1 <- filter(select_gene, positive_num >=2)


# here we only choose the pathway with most genes
subsystem <- as.data.frame(table(select_gene_filter1$pathway), stringsAsFactors = FALSE)
subsystem_filter <- subsystem[subsystem$Freq>=3,]
colnames(subsystem_filter) <- c("subsystem","num")

# set the factor level by the mean value of occur number
#Factor <-Result0 %>% group_by(pathway) %>% summarise(median=mean(Occur_num))
subsystem_filter <- subsystem_filter[order(subsystem_filter$num),]
subsystem_filter$subsystem <-factor(subsystem_filter$subsystem, levels=subsystem_filter$subsystem)
ggplot(subsystem_filter, aes(x=subsystem, y=num)) +
  geom_bar(stat = "identity", fill = "#FF6666") +
  xlab('') + ylab('Number of genes') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=8, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) #+



#################################################################################################################
# enrichment analysis of positive selected genes based on the reference genome of sce s288c
# new result for all the filtered OGs
gene_dn_ds_all_new <- read_csv("data/gene_dn_ds_03_02.csv") %>%
  filter(., !is.na(dN_dS)) %>%
  filter(., dN_dS < 10)

ortholog <- read_tsv("data/ortholog_occurence_num_all.tsv")
ortholog0 <- filter(ortholog, with_sce == "with_sce")
gene_dn_ds_all_new$OG <- str_replace_all(gene_dn_ds_all_new$OG, ".out_yn00", "")
gene_dn_ds_sce <- gene_dn_ds_all_new[gene_dn_ds_all_new$OG %in% ortholog0$ID, ]
# find the sce gene id for each OG
og_panID_mapping <- read_tsv("data/representatives.tsv")
gene_dn_ds_sce0 <- left_join(gene_dn_ds_sce, og_panID_mapping, by = c("OG" = "ortho_id"))
gene_dn_ds_sce0$locus <- str_replace_all(gene_dn_ds_sce0$representative, "Saccharomyces_cerevisiae@", "")
# connect the positive selected number for the OG with sce seq
fubar_02_08 <- read_csv("data/fubar_02_08.csv") %>% filter(.,positive_num >=1)
gene_dn_ds_sce0$positive_num <- getSingleReactionFormula(fubar_02_08$positive_num,fubar_02_08$OG,gene_dn_ds_sce0$OG)
gene_dn_ds_sce0$positive_num <- as.numeric(gene_dn_ds_sce0$positive_num)
# connect the sce gene name from uniprot
uniprotGeneID_mapping <- read_excel("data/uniprotGeneID_mapping.xlsx")
gene_dn_ds_sce0$uniprotID <- getSingleReactionFormula(uniprotGeneID_mapping$Entry,uniprotGeneID_mapping$GeneName,gene_dn_ds_sce0$locus)
gene_dn_ds_sce0$geneID <- getSingleReactionFormula(uniprotGeneID_mapping$`Entry name`,uniprotGeneID_mapping$GeneName,gene_dn_ds_sce0$locus)
gene_dn_ds_sce0$geneID <- str_replace_all(gene_dn_ds_sce0$geneID, "_YEAST", "")


# filter sce gene list based on positive number for the enrichment analysis
gene_inf <- filter(gene_dn_ds_sce0, positive_num >=2) # here if we choose the postive_num >=3, there are few significant enrichemnt GO/pathway terms.
gene_list <- paste0(gene_inf$locus,collapse = ",")
print(gene_list) # for DAVID
write.table(gene_dn_ds_sce0, "result/sce_gene_list.txt", row.names = FALSE, sep = "\t")
