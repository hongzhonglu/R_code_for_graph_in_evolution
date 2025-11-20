library(readxl)
library(hongR)
library(stringr)
library(readr)
library(ggplot2)
library(dplyr)
library(tidyverse)

###################################################################################################
# new result for all the filtered OGs
gene_dn_ds_all_new <- read_csv("data/gene_dn_ds_03_02.csv") %>%
  filter(., !is.na(dN_dS)) %>%
  filter(., dN_dS < 10)

ortholog <- read_tsv("data/ortholog_occurence_num_all.tsv")
ortholog0 <- filter(ortholog, with_sce == "with_sce")
gene_dn_ds_all_new$OG <- str_replace_all(gene_dn_ds_all_new$OG, ".out_yn00", "")
gene_dn_ds_sce <- gene_dn_ds_all_new[gene_dn_ds_all_new$OG %in% ortholog0$ID, ]
# find the sce gene id for each OG
og_sce_mapping <- read_tsv("data/representatives.tsv")
gene_dn_ds_sce0 <- left_join(gene_dn_ds_sce, og_sce_mapping, by = c("OG" = "ortho_id"))
gene_dn_ds_sce0$locus <- str_replace_all(gene_dn_ds_sce0$representative, "Saccharomyces_cerevisiae@", "")
write.csv(gene_dn_ds_sce0, "result/dn_ds_sce.csv")



# choose the top 5% of all genes with highest dN/dS and lowest dN/dS
dn_ds_order <- gene_dn_ds_sce0[order(gene_dn_ds_sce0$dN_dS), ]
plot(density(dn_ds_order$dN_dS))
sample_num <- 222 #4459*0.05
lowest_dn_ds_group <- dn_ds_order[c(1:sample_num),]
highest_dn_ds_group <- dn_ds_order[c((4459-sample_num+1):4459),]
print(paste0(lowest_dn_ds_group$locus, collapse = ","))
print(paste0(highest_dn_ds_group$locus, collapse = ","))




################# check Saccharomyces_cerevisiae-FCC and dN/dS ####################################
kinetics_analysis <- read.table("protein_related_parameters/FCC of sce enzymes/KcatSensitivities_YEP.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(kinetics_analysis) <- c("locus", "glucose", "acetate", "ethanol", "glycerol", "sorbitol", "galactose", "ribose", "xylose")
# merge the nsSNP number with the kinetics data

geneGEM01 <- merge(gene_dn_ds_sce0, kinetics_analysis, by.x = "locus", by.y = "locus")
substrate_list <- c("glucose", "acetate", "ethanol", "glycerol", "sorbitol", "galactose", "ribose", "xylose")

#function to plot x y
plotSNP_FCC <- function(x.element, y.element='dN_dS', data=geneGEM01){

  ggplot(data, aes_string(x=x.element, y=y.element)) +
    geom_bin2d(bins = 60) +
    labs(x="",y=y.element) +
    scale_fill_gradientn(limits=c(0,40), breaks=seq(0, 40, by=10), colours=rev(rainbow(4)))+
    theme(axis.text=element_text(size=20, family="Arial"),
          axis.title=element_text(size=24, family="Arial") ) +
    ggtitle(x.element) +
    theme(plot.title=element_text(hjust = 0.8, size=20,family="Arial"))+
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 1)) +
    theme(legend.text =element_text(size=16, family="Arial"),
          legend.title =element_text(size=16, family="Arial"))+
    #theme(legend.position="none")
  ggsave(out <- paste('result/',x.element,'.tiff', sep = ""), width=8, height=6, dpi=300)

}

for (i in 1:8){
  print(substrate_list[i])
  plotSNP_FCC(x.element=substrate_list[i], y.element='dN_dS', data=geneGEM01)
}

# for the glucose as the substrate
geneGEM01$FCC_type <- NA
geneGEM01$FCC_type[geneGEM01$glucose==0] <- "A.no effect"
geneGEM01$FCC_type[geneGEM01$glucose >0 & geneGEM01$glucose <= 0.001] <- "B.0-0.001"
geneGEM01$FCC_type[geneGEM01$glucose >0.001 & geneGEM01$glucose <= 0.01] <- "C.0.001-0.01"
geneGEM01$FCC_type[geneGEM01$glucose >0.01] <- "D.0.01-0.05"
# plot
ggplot(geneGEM01,aes(x=FCC_type, y=dN_dS, fill=FCC_type)) + 
  stat_boxplot(geom ='errorbar', width = 0.25) + # add caps
  geom_boxplot() +
  xlab('Flux control coefficient') + ylab('dN/dS') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('')

ggsave(out <- paste('result/','Relation_between_FCCs_and_dN_dS_distribution','.svg', sep = ""), width=8, height=6, dpi=600)


g1 <- geneGEM01$dN_dS[geneGEM01$FCC_type == "A.no effect"]
g2 <- geneGEM01$dN_dS[geneGEM01$FCC_type == "B.0-0.001"]
g3 <- geneGEM01$dN_dS[geneGEM01$FCC_type == "C.0.001-0.01"]
g4 <- geneGEM01$dN_dS[geneGEM01$FCC_type == "D.0.01-0.05"]

# t.test
t.test(g1,g2)
t.test(g1,g3)
t.test(g1,g4)
# wilcon.test
wilcox.test(g1,g2, alternative = "two.sided")
wilcox.test(g1,g3, alternative = "two.sided")
wilcox.test(g1,g4, alternative = "two.sided")



################# check Kluyveromyces_marxianus-FCC and dN/dS ####################################
K_marxianus_ID2 <- read_tsv("data_for_K_marxianus/Kluyveromyces_marxianus.tsv")
K_marxianus_fcc <- read_excel("data_for_K_marxianus/Kluyveromyces_marxianus_limGrowth.xlsx")


K_marxianus_ID1 <- read_delim("data_for_K_marxianus/Kluyveromyces marxianus.csv", "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(K_marxianus_ID1) <- c("queryID","UniprotID","identity")
K_marxianus_ID10 <- K_marxianus_ID1 %>%
  separate(UniprotID, c("s1", "s2", "UniprotID"), "\\|")
K_marxianus_ID10$UniprotID <- str_replace_all(K_marxianus_ID10$UniprotID,"_KLUMD","")
K_marxianus_ID10$UniprotID <- paste("prot_", K_marxianus_ID10$UniprotID, sep = "")
K_marxianus_fcc$queryID <- getSingleReactionFormula(K_marxianus_ID10$queryID,K_marxianus_ID10$UniprotID,K_marxianus_fcc$Var1)
K_marxianus_fcc$OG <- getSingleReactionFormula(K_marxianus_ID2$ortholog_id,K_marxianus_ID2$protein_id,K_marxianus_fcc$queryID)

# input other metabolic genes in the model
K_marxianus_GEM_all_gene <- read_excel("data_for_K_marxianus/K_marxianus_GEM_all_gene.xlsx")
K_marxianus_GEM_all_gene$gene_GEM <- str_replace_all(K_marxianus_GEM_all_gene$gene_GEM,"_KLUMD","")
K_marxianus_GEM_all_gene$gene_GEM <- paste("prot_", K_marxianus_GEM_all_gene$gene_GEM, sep = "")

K_marxianus_GEM_all_gene1 <- K_marxianus_GEM_all_gene[str_detect(K_marxianus_GEM_all_gene$gene_GEM, "@"), ]
K_marxianus_GEM_all_gene2 <- K_marxianus_GEM_all_gene[!str_detect(K_marxianus_GEM_all_gene$gene_GEM, "@"), ]
colnames(K_marxianus_GEM_all_gene1) <- c("queryID")
K_marxianus_GEM_all_gene1$queryID <- str_replace_all(K_marxianus_GEM_all_gene1$queryID, "prot_", "")
K_marxianus_GEM_all_gene1$OG <- getSingleReactionFormula(K_marxianus_ID2$ortholog_id,K_marxianus_ID2$protein_id,K_marxianus_GEM_all_gene1$queryID)


K_marxianus_GEM_all_gene2$queryID <- getSingleReactionFormula(K_marxianus_ID10$queryID,K_marxianus_ID10$UniprotID,K_marxianus_GEM_all_gene2$gene_GEM)
K_marxianus_GEM_all_gene2$OG <- getSingleReactionFormula(K_marxianus_ID2$ortholog_id,K_marxianus_ID2$protein_id,K_marxianus_GEM_all_gene2$queryID)
K_marxianus_GEM_all_gene2 <- K_marxianus_GEM_all_gene2[,c("queryID","OG")]
K_marxianus_GEM_all_new <- rbind.data.frame(K_marxianus_GEM_all_gene1, K_marxianus_GEM_all_gene2)
K_marxianus_GEM_all_new$dN_dS <- getSingleReactionFormula(gene_dn_ds_all_new$dN_dS,gene_dn_ds_all_new$OG,K_marxianus_GEM_all_new$OG)
K_marxianus_GEM_all_new <- K_marxianus_GEM_all_new[K_marxianus_GEM_all_new$dN_dS !="NA", ]
K_marxianus_GEM_all_new$dN_dS <- as.numeric(K_marxianus_GEM_all_new$dN_dS)

K_marxianus_GEM_all_new$fcc <- getSingleReactionFormula(K_marxianus_fcc$Var5, K_marxianus_fcc$OG, K_marxianus_GEM_all_new$OG)
K_marxianus_GEM_all_new$fcc <- as.numeric(K_marxianus_GEM_all_new$fcc)
K_marxianus_GEM_all_new$fcc[is.na(K_marxianus_GEM_all_new$fcc)] <- 0


# group
K_marxianus_GEM_all_new$FCC_type <- NA
K_marxianus_GEM_all_new$FCC_type[K_marxianus_GEM_all_new$fcc==0] <- "A.no effect"
K_marxianus_GEM_all_new$FCC_type[K_marxianus_GEM_all_new$fcc >0 & K_marxianus_GEM_all_new$fcc <= 0.001] <- "B.0-0.001"
K_marxianus_GEM_all_new$FCC_type[K_marxianus_GEM_all_new$fcc >0.001 & K_marxianus_GEM_all_new$fcc <= 0.01] <- "C.0.001-0.01"
K_marxianus_GEM_all_new$FCC_type[K_marxianus_GEM_all_new$fcc >0.01] <- "D.0.01-0.05"

# plot
ggplot(K_marxianus_GEM_all_new,aes(x=FCC_type, y=dN_dS, fill=FCC_type)) + 
  stat_boxplot(geom ='errorbar', width = 0.25) + # add caps
  geom_boxplot() +
  xlab('Flux control coefficient') + ylab('dN/dS') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('')

g1 <- K_marxianus_GEM_all_new$dN_dS[K_marxianus_GEM_all_new$FCC_type == "A.no effect"]
g2 <- K_marxianus_GEM_all_new$dN_dS[K_marxianus_GEM_all_new$FCC_type == "B.0-0.001"]
g3 <- K_marxianus_GEM_all_new$dN_dS[K_marxianus_GEM_all_new$FCC_type == "C.0.001-0.01"]
g4 <- K_marxianus_GEM_all_new$dN_dS[K_marxianus_GEM_all_new$FCC_type == "D.0.01-0.05"]

# t.test
t.test(g1,g2)
t.test(g1,g3)
t.test(g1,g4)
# wilcon.test
wilcox.test(g1,g2, alternative = "two.sided")
wilcox.test(g1,g3, alternative = "two.sided")
wilcox.test(g1,g4, alternative = "two.sided")





################# check Schizosaccharomyces_pombe-FCC and dN/dS ####################################
S_pombe_ID2 <- read_tsv("data_for_S_pombe/Schizosaccharomyces_pombe.tsv")
S_pombe_fcc <- read_excel("data_for_S_pombe/S_pombe_limGrowth.xlsx")
S_pombe_ID1 <- read_delim("data_for_S_pombe/Schizosaccharomyces pombe.csv", "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(S_pombe_ID1) <- c("queryID","UniprotID","identity")
S_pombe_ID10 <- S_pombe_ID1 %>%
  separate(UniprotID, c("s1",  "UniprotID", "s2"), "\\|")
#S_pombe_ID10$UniprotID <- str_replace_all(S_pombe_ID10$UniprotID,"_KLUMD","")
S_pombe_ID10$UniprotID <- paste("prot_", S_pombe_ID10$UniprotID, sep = "")
S_pombe_fcc$queryID <- getSingleReactionFormula(S_pombe_ID10$queryID,S_pombe_ID10$UniprotID,S_pombe_fcc$Var1)
S_pombe_fcc$OG <- getSingleReactionFormula(S_pombe_ID2$ortholog_id,S_pombe_ID2$protein_id,S_pombe_fcc$queryID)

# input other metabolic genes in the model
S_pombe_GEM_all_gene <- read_excel("data_for_S_pombe/S_pombe_GEM_all_gene.xlsx")
#S_pombe_GEM_all_gene$gene_GEM <- str_replace_all(S_pombe_GEM_all_gene$gene_GEM,"_KLUMD","")
S_pombe_GEM_all_gene$gene_GEM <- paste("prot_", S_pombe_GEM_all_gene$gene_GEM, sep = "")

S_pombe_GEM_all_gene1 <- S_pombe_GEM_all_gene[str_detect(S_pombe_GEM_all_gene$gene_GEM, "@"), ]
S_pombe_GEM_all_gene2 <- S_pombe_GEM_all_gene[!str_detect(S_pombe_GEM_all_gene$gene_GEM, "@"), ]
colnames(S_pombe_GEM_all_gene1) <- c("queryID")
S_pombe_GEM_all_gene1$queryID <- str_replace_all(S_pombe_GEM_all_gene1$queryID, "prot_", "")
S_pombe_GEM_all_gene1$OG <- getSingleReactionFormula(S_pombe_ID2$ortholog_id,S_pombe_ID2$protein_id,S_pombe_GEM_all_gene1$queryID)


S_pombe_GEM_all_gene2$queryID <- getSingleReactionFormula(S_pombe_ID10$queryID,S_pombe_ID10$UniprotID,S_pombe_GEM_all_gene2$gene_GEM)
S_pombe_GEM_all_gene2$OG <- getSingleReactionFormula(S_pombe_ID2$ortholog_id,S_pombe_ID2$protein_id,S_pombe_GEM_all_gene2$queryID)
S_pombe_GEM_all_gene2 <- S_pombe_GEM_all_gene2[,c("queryID","OG")]
S_pombe_GEM_all_new <- rbind.data.frame(S_pombe_GEM_all_gene1, S_pombe_GEM_all_gene2)

S_pombe_GEM_all_new$dN_dS <- getSingleReactionFormula(gene_dn_ds_all_new$dN_dS,gene_dn_ds_all_new$OG,S_pombe_GEM_all_new$OG)
S_pombe_GEM_all_new <- S_pombe_GEM_all_new[S_pombe_GEM_all_new$dN_dS !="NA", ]
S_pombe_GEM_all_new$dN_dS <- as.numeric(S_pombe_GEM_all_new$dN_dS)

S_pombe_GEM_all_new$fcc <- getSingleReactionFormula(S_pombe_fcc$Var5, S_pombe_fcc$OG, S_pombe_GEM_all_new$OG)
S_pombe_GEM_all_new$fcc <- as.numeric(S_pombe_GEM_all_new$fcc)
S_pombe_GEM_all_new$fcc[is.na(S_pombe_GEM_all_new$fcc)] <- 0

# group
S_pombe_GEM_all_new$FCC_type <- NA
S_pombe_GEM_all_new$FCC_type[S_pombe_GEM_all_new$fcc==0] <- "A.no effect"
S_pombe_GEM_all_new$FCC_type[S_pombe_GEM_all_new$fcc >0 & S_pombe_GEM_all_new$fcc <= 0.001] <- "B.0-0.001"
S_pombe_GEM_all_new$FCC_type[S_pombe_GEM_all_new$fcc >0.001 & S_pombe_GEM_all_new$fcc <= 0.01] <- "C.0.001-0.01"
S_pombe_GEM_all_new$FCC_type[S_pombe_GEM_all_new$fcc >0.01] <- "D.0.01-0.05"

# plot
ggplot(S_pombe_GEM_all_new,aes(x=FCC_type, y=dN_dS, fill=FCC_type)) + 
  stat_boxplot(geom ='errorbar', width = 0.25) + # add caps
  geom_boxplot() +
  xlab('Flux control coefficient') + ylab('dN/dS') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('')

g1 <- S_pombe_GEM_all_new$dN_dS[S_pombe_GEM_all_new$FCC_type == "A.no effect"]
g2 <- S_pombe_GEM_all_new$dN_dS[S_pombe_GEM_all_new$FCC_type == "B.0-0.001"]
g3 <- S_pombe_GEM_all_new$dN_dS[S_pombe_GEM_all_new$FCC_type == "C.0.001-0.01"]
g4 <- S_pombe_GEM_all_new$dN_dS[S_pombe_GEM_all_new$FCC_type == "D.0.01-0.05"]

# t.test
t.test(g1,g2)
t.test(g1,g3)
t.test(g1,g4)
# wilcon.test
wilcox.test(g1,g2, alternative = "two.sided")
wilcox.test(g1,g3, alternative = "two.sided")
wilcox.test(g1,g4, alternative = "two.sided")






############ try to combine three species result#####################################
gene_sce <- geneGEM01[,c("dN_dS","FCC_type")]
gene_sce$species <- "S. cerevisiae"

gene_pombe <- S_pombe_GEM_all_new[,c("dN_dS","FCC_type")]
gene_pombe$species <- "S. pombe"

gene_marxianus <- K_marxianus_GEM_all_new[,c("dN_dS","FCC_type")]
gene_marxianus$species <- "K. marxianus"
# combine
three_species <- rbind.data.frame(gene_sce, gene_pombe, gene_marxianus)

three_species$FCC_type <- as.factor(three_species$FCC_type)

pd = position_dodge(width = 0.75)
ggplot(three_species, aes(x=FCC_type, y=dN_dS, fill=species)) + 
  stat_boxplot(geom="errorbar", position=pd, width=0.2) +# add caps +
  geom_boxplot()+
  xlab('Flux control coefficient') + ylab('dN/dS') +
  ylim(0,0.5)+
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = c(0.8,0.8)) +
  theme(axis.text=element_text(size=16, family="Arial"),
        axis.title=element_text(size=20,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('')

write.csv(three_species, "result/dataset_for_new_Fig.3C.csv")










#################### not used ################################################################
#correlation analysis between protein abundence and nsSNP number
protein_abundence <- read.table('protein_related_parameters/protein_abundance_sce/protein_abundence_paxdb.txt', header = TRUE, stringsAsFactors = FALSE)
protein_abundence$string_external_id <- str_replace_all(protein_abundence$string_external_id,'4932.','')
gene_dn_ds_sce0$abundence <- getSingleReactionFormula(protein_abundence$abundance,protein_abundence$string_external_id,gene_dn_ds_sce0$locus)
gene_dn_ds_sce0$abundence <- as.numeric(gene_dn_ds_sce0$abundence)
# plot
# plot(density(gene_dn_ds_sce0$abundence[!is.na(gene_dn_ds_sce0$abundence)]))
ggplot(gene_dn_ds_sce0,  aes(x=abundence)) +
  xlim(0,1000)+
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)

ggplot(gene_dn_ds_sce0, aes(abundence, dN_dS)) +
  geom_bin2d(bins = 80) +
  labs(x="Abundence of protein in each gene", y="dN_dS in each gene") +
  scale_fill_gradientn(limits=c(0,200), breaks=seq(0, 200, by=50), colours=rev(rainbow(4)))+
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1)) +
  theme(legend.text =element_text(size=16, family="Arial"),
        legend.title =element_text(size=16, family="Arial"))
ggsave(out <- paste('result/','dN_dS and abundence','.tiff', sep = ""), width=8, height=6, dpi=300)
# change the scatter into the column plot
gene_dn_ds_sce0$abundence_type <- NA
gene_dn_ds_sce0$abundence_type[gene_dn_ds_sce0$abundence<=50] <-"A.< 50"
gene_dn_ds_sce0$abundence_type[gene_dn_ds_sce0$abundence >50 & gene_dn_ds_sce0$abundence<=100] <-"B.50-100"
gene_dn_ds_sce0$abundence_type[gene_dn_ds_sce0$abundence >100 & gene_dn_ds_sce0$abundence<=150] <-"C.100-150"
gene_dn_ds_sce0$abundence_type[gene_dn_ds_sce0$abundence >150 & gene_dn_ds_sce0$abundence<=200] <-"D.150-200"
gene_dn_ds_sce0$abundence_type[gene_dn_ds_sce0$abundence >200 & gene_dn_ds_sce0$abundence<=250] <-"E.200-250"
gene_dn_ds_sce0$abundence_type[gene_dn_ds_sce0$abundence >250 & gene_dn_ds_sce0$abundence<=500] <-"F.250-500"
gene_dn_ds_sce0$abundence_type[gene_dn_ds_sce0$abundence >500] <-"G.>500"
gene_dn_ds_sce_pro <- gene_dn_ds_sce0[!is.na(gene_dn_ds_sce0$abundence_type),]
# plot
ggplot(gene_dn_ds_sce_pro,aes(x=abundence_type, y=dN_dS, fill=abundence_type)) + geom_boxplot() +
  xlab('') + ylab('dN_dS') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=8, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) #+
g1 <- gene_dn_ds_sce_pro$dN_dS[gene_dn_ds_sce_pro$abundence_type == "A.< 50"]
g2 <- gene_dn_ds_sce_pro$dN_dS[gene_dn_ds_sce_pro$abundence_type == "C.100-150"]
t.test(g1,g2)