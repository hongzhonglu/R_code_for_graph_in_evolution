library(readxl)
library(hongR)
library(stringr)
library(readr)
library(ggplot2)
library(dplyr)


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


# input the result for the kinetics sensitivity analysis
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
ggplot(geneGEM01,aes(x=FCC_type, y=dN_dS, fill=FCC_type)) + geom_boxplot() +
  xlab('') + ylab('dN_dS') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=8, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) #+
g1 <- geneGEM01$dN_dS[geneGEM01$FCC_type == "A.no effect"]
g2 <- geneGEM01$dN_dS[geneGEM01$FCC_type == "C.0.001-0.01"]
t.test(g1,g2)


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