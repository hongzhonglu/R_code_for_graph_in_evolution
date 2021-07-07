# note------------
# here we try to analyze the dN/dS distribution of genes from each sub-pathway for sce!
# 2020-06-02
# Hongzhong Lu

library(readxl)
library(hongR)
library(ggplot2)
library(dplyr)
# input the pathway definition from kegg and reactome database
# the original data is downloaded on 2019-06-03


############################################################################################################
# for sce gene only
sce_pathway_kegg <- read.table("protein_related_parameters/sce_pathway_kegg.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# merge the subsystem annotation with the manual transporter annotation
transport_gene <- read_tsv("protein_related_parameters/TransRxnGeneAnnotation.tsv")
transport_gene0 <- transport_gene[, c("gene")]
colnames(transport_gene0) <- "geneID"
transport_gene0$pathwayID <- "Transport"
transport_gene0$pathwayName <- "Transporter"

sce_pathway_kegg <- rbind.data.frame(sce_pathway_kegg, transport_gene0)


# here we just mapping dN_dS onto the existing gene-pathway relation
gene_dn_ds_sce0 <- read.table("result/sce_gene_list.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sce_pathway_kegg$dN_dS <- getSingleReactionFormula(gene_dn_ds_sce0$dN_dS, gene_dn_ds_sce0$locus, sce_pathway_kegg$geneID)
sce_pathway_kegg$dN_dS <- as.numeric(sce_pathway_kegg$dN_dS)
sce_pathway_kegg <- filter(sce_pathway_kegg, dN_dS != "NA") # remove the gene with NA SNP

# random sampling
sample_value <- sample(sce_pathway_kegg$dN_dS, size = 3545)
dN_dS_test <- data_frame(dN_dS = sample_value)
ggplot(dN_dS_test, aes(x = dN_dS, y = ..density..)) +
  geom_histogram(fill = "blue", colour = "blue", size = .2) +
  geom_density()

# test
s <- "Oxidative phosphorylation"
print(s)
gene_pathway0 <- sce_pathway_kegg %>% filter(., pathwayName == s)
ggplot(gene_pathway0, aes(x = dN_dS, y = ..density..)) +
  geom_histogram(fill = "blue", colour = "blue", size = .2) +
  geom_density() +
  ggtitle(s)


dN_dS_subsytem <- sce_pathway_kegg %>%
  group_by(pathwayName) %>%
  dplyr::summarize(Mean = mean(dN_dS, na.rm = TRUE))

# plot the result for all the subsystmes
for (i in seq_along(dN_dS_subsytem$pathwayName)) {
  print(i)
  s <- dN_dS_subsytem$pathwayName[i]
  print(s)
  gene_pathway0 <- sce_pathway_kegg %>% filter(., pathwayName == s)
  gene_pathway0 <- gene_pathway0[, c("pathwayName", "dN_dS")]
  sample_value <- sample(sce_pathway_kegg$dN_dS, size = 3545)
  dN_dS_test <- data_frame(pathwayName = "Random sampling", dN_dS = sample_value)
  # obtain the merged result
  result_merge <- rbind.data.frame(gene_pathway0, dN_dS_test)

  ggplot(result_merge, aes(dN_dS, fill = pathwayName, colour = pathwayName)) +
    geom_density(alpha = 0.2) +
    labs(x = "dN_dS in each gene") +
    theme(legend.position = c(0.85, 0.85)) +
    theme(
      axis.text = element_text(size = 20, face = "bold", family = "Arial"),
      axis.title = element_text(size = 24, face = "bold", family = "Arial"),
      legend.text = element_text(size = 13, family = "Arial")
    ) + ggtitle("") +
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 1))
  s <- str_replace_all(s, "\\ / ", "_")
  ggsave(out <- paste("result/pathway_analysis_dN_dS/", s, ".png", sep = ""), width = 8, height = 6, dpi = 300)
}



####################################################################################################################
# for all the genes
KO_function <- read.table("kegg/ko_function.txt", sep = "\t", stringsAsFactors = FALSE)
KO_pathway <- read.table("kegg/ko_pathway.txt", sep = "\t", stringsAsFactors = FALSE)
pathway <- read.table("kegg/pathway_list_kegg.txt", sep = "\t", stringsAsFactors = FALSE)
KO_pathway <- filter(KO_pathway,str_detect(KO_pathway$V1, "map"))
KO_function$V1 <- str_replace_all(KO_function$V1,"ko:","")
KO_pathway$V2 <- str_replace_all(KO_pathway$V2,"ko:","")

# classification of OGs with selected sites based on KO or function annotation.
og_panID_mapping <- read_tsv("data/representatives.tsv")
# get the KO id
panID_KO <- read.table("kegg/ko_annotation_of_all_OG.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

gene_dn_ds_all_new <- read_csv("data/gene_dn_ds_03_02.csv") %>%
  filter(., !is.na(dN_dS)) %>%
  filter(., dN_dS < 10)
gene_dn_ds_all_new$OG <- str_replace_all(gene_dn_ds_all_new$OG, ".out_yn00", "")
# get the representative id, ko id and pathway information
# it should be careful that a ko could multiple pathways
gene_dn_ds_all_new$representative <-  getSingleReactionFormula(og_panID_mapping$representative, og_panID_mapping$ortho_id, gene_dn_ds_all_new$OG)
gene_dn_ds_all_new$KO <- getSingleReactionFormula(panID_KO$ko, panID_KO$query, gene_dn_ds_all_new$representative)
gene_dn_ds_all_new <- gene_dn_ds_all_new[gene_dn_ds_all_new$KO!="NA",]
gene_dn_ds_all_new$pathwayID <- getMultipleReactionFormula(KO_pathway$V1,KO_pathway$V2, gene_dn_ds_all_new$KO)
gene_dn_ds_all_new$ko_function <- getMultipleReactionFormula(KO_function$V2,KO_function$V1, gene_dn_ds_all_new$KO)
# establish the mapping between the dN/dS and pathway (transcript factor)
dn_ds_pathway <- gene_dn_ds_all_new[!is.na(gene_dn_ds_all_new$pathwayID),]
dn_ds_pathway <- dn_ds_pathway[, c("OG","pathwayID")]
dn_ds_pathway1 <- splitAndCombine(dn_ds_pathway$pathwayID, dn_ds_pathway$OG,";")
dn_ds_pathway1$pathway <- getSingleReactionFormula(pathway$V2, pathway$V1, dn_ds_pathway1$v1)
dn_ds_pathway1$pathway[dn_ds_pathway1$pathway=="NA"] <- NA
dn_ds_pathway1$pathway <- str_replace_all(dn_ds_pathway1$pathway, " - fly", "") %>% str_replace_all(.," - other", "") %>%
  str_replace_all(.," - yeast", "") %>% str_replace_all(.," - Caulobacter", "") %>% str_replace_all(.," - animal", "")
dn_ds_pathway1 <- dn_ds_pathway1[!is.na(dn_ds_pathway1$pathway),]
# as the first step, we will focus on the genes from metabolic pathway
OG_metabolic_pathway <- dn_ds_pathway1$v2[dn_ds_pathway1$v1=="path:map01100"]
dn_ds_metabolic_pathway <- dn_ds_pathway1[dn_ds_pathway1$v2 %in% OG_metabolic_pathway, ]
# get the dn/ds
dn_ds_metabolic_pathway$dN_dS <- getSingleReactionFormula(gene_dn_ds_all_new$dN_dS, gene_dn_ds_all_new$OG, dn_ds_metabolic_pathway$v2)
dn_ds_metabolic_pathway$dN_dS <- as.numeric(dn_ds_metabolic_pathway$dN_dS)
dn_ds_metabolic_pathway$pathway <- as.factor(dn_ds_metabolic_pathway$pathway)
dn_ds_summary1 <- group_by(dn_ds_metabolic_pathway, pathway) %>% summarize(m = mean(dN_dS), sd = sd(dN_dS))
# here we only choose several interesting pathways
interest_pathway <- c("Citrate cycle (TCA cycle)","Glycolysis / Gluconeogenesis",
                      "Pentose phosphate pathway", "Oxidative phosphorylation",
                      "Biosynthesis of unsaturated fatty acids","Biosynthesis of amino acids",
                      "Biosynthesis of antibiotics",
                      "Biosynthesis of secondary metabolites",
                      "Purine metabolism", "Pyrimidine metabolism")
more_pathway <- c("Biosynthesis of secondary metabolites", "Pyrimidine metabolism", "Biosynthesis of amino acids","Purine metabolism","Carbon metabolism","Glycolysis / Gluconeogenesis","Oxidative phosphorylation")

interest_pathway2 <- unique(c(interest_pathway, more_pathway ))


dn_ds_metabolic_pathway1 <- dn_ds_metabolic_pathway[dn_ds_metabolic_pathway$pathway %in% interest_pathway2, ]

dn_ds_summary2 <- group_by(dn_ds_metabolic_pathway1, pathway) %>% summarize(m = median(dN_dS))

dn_ds_summary2 <- dn_ds_summary2[order(dn_ds_summary2$m, decreasing = FALSE),]

dn_ds_metabolic_pathway1$pathway <-factor(dn_ds_metabolic_pathway1$pathway, levels=dn_ds_summary2$pathway)



# plot
ggplot(dn_ds_metabolic_pathway1 ,aes(x=pathway, y=dN_dS)) + 
  stat_boxplot(geom ='errorbar', width = 0.25) + # add caps
  geom_boxplot(fill = "#FF6666") +
  xlab('') + ylab('dN/dS') +
  #theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=16,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", colour = "black"))  +
  coord_flip()


G1 <- filter(dn_ds_metabolic_pathway1, dn_ds_metabolic_pathway1$pathway == "Citrate cycle (TCA cycle)")
G2<- filter(dn_ds_metabolic_pathway1, dn_ds_metabolic_pathway1$pathway == "Pentose phosphate pathway")
G3 <- filter(dn_ds_metabolic_pathway1, dn_ds_metabolic_pathway1$pathway == "Biosynthesis of amino acids")
G4 <- filter(dn_ds_metabolic_pathway1, dn_ds_metabolic_pathway1$pathway == "Biosynthesis of antibiotics")
G5 <- filter(dn_ds_metabolic_pathway1, dn_ds_metabolic_pathway1$pathway == "Glycolysis / Gluconeogenesis")


# t.test
t.test(G1$dN_dS, G2$dN_dS)
t.test(G1$dN_dS, G3$dN_dS)
t.test(G1$dN_dS, G4$dN_dS)
t.test(G1$dN_dS, G5$dN_dS)

# wilcon.test
wilcox.test(G1$dN_dS, G2$dN_dS, alternative = "two.sided")
wilcox.test(G1$dN_dS, G3$dN_dS, alternative = "two.sided")
wilcox.test(G1$dN_dS, G4$dN_dS, alternative = "two.sided")
wilcox.test(G1$dN_dS, G5$dN_dS, alternative = "two.sided")
# By comparison, it could find differences between t.test and wilcon.test






