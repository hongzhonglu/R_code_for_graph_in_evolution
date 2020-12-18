# note------------
# Compare the site model results from different combination of methods
# 2020-09-09
# Hongzhong Lu


library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(hongR)
library(readxl)
library(VennDiagram)

# align by macse for pruned Ortholog group
fubar_macse_prune <- read_csv("data/site_model_data/fubar_02_08.csv")
fel_macse_prune <- read_csv("data/site_model_data/fel_macse_prune_3400.csv")
fel_macse_prune$OG <- str_replace_all(fel_macse_prune$OG, ".csv", "")

# align by macse for pruned Ortholog group
fubar_guidance_prune <- read_csv("data/site_model_data/fubar_guidance_prune.csv")
fel_guidance_prune <- read_csv("data/site_model_data/fel_guidance_prune.csv")
fel_guidance_prune$OG <- str_replace_all(fel_guidance_prune$OG , ".csv", "")

# align by macse for unpruned Ortholog group
fubar_guidance_unprune <- read_csv("data/site_model_data/fubar_guidance_unprune.csv")


# filter based on the sites number
site_num <- 1
fubar_macse_prune0 <- filter(fubar_macse_prune, positive_num >= site_num)
fel_macse_prune0 <- filter (fel_macse_prune, site_pvalue_0.1 >= site_num)
#fel_macse_prune0 <- filter (fel_macse_prune, site_pvalue_0.05 >= site_num)
fubar_guidance_prune0 <- filter(fubar_guidance_prune, positive_num >= site_num)
fel_guidance_prune0 <- filter (fel_guidance_prune, site_pvalue_0.1 >= site_num)
#fel_guidance_prune0 <- filter (fel_guidance_prune, site_pvalue_0.05 >= site_num)
fubar_guidance_unprune0 <- filter(fubar_guidance_unprune, positive_num >= site_num)


# prepare the plot
type <- c("fubar_macse_prune", "fubar_guidance_prune","fubar_guidance_unprune","fel_macse_prune", 
          "fel_guidance_prune")
number <- c(length(fubar_macse_prune0$OG), length(fubar_guidance_prune0$OG), length(fubar_guidance_unprune0$OG),
            length(fel_macse_prune0$OG), length(fel_guidance_prune0$OG))
all_result <- data.frame(source=type, num=number, stringsAsFactors = FALSE)
# remove the fel_macse_prune as it is only analyze the result from 3400 OGs
all_result0 <- filter(all_result, source !="fel_macse_prune")

ggplot(data=all_result0, aes(x=source, y=num)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(legend.position = "none") #+
  #theme(panel.background = element_rect(fill = "white", color="black", size = 1)


# insect analysis
# here the cds aligned by guidance is used.
VennDiagram::venn.diagram(x= list(fubar_guidance_prune = fubar_guidance_prune0$OG,
                                  fubar_guidance_unprune = fubar_guidance_unprune0$OG, fel_guidance_prune =fel_guidance_prune0$OG), 
                          filename = "result/site_model_results_using_different_method_combination3.png", height = 1000, width = 1000,resolution =300, imagetype="png", col="transparent",
                          fill=c("blue","red", "green"),alpha = 0.50, cex=0.45, cat.cex=0.45)

# there are 862 positive selected genes with selected amino acid sites
Combine0  <- intersect(fubar_guidance_prune0$OG, fubar_guidance_unprune0$OG)
Combine1 <- intersect(Combine0, fel_guidance_prune0$OG)


#####################################################################################
# further analysis these OGs with positive selected sites
gene_dn_ds_all_new <- read.table("result/gene_dn_ds_all_new.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gene_dn_ds_all_new$Type <- NA
gene_dn_ds_all_new$Type[gene_dn_ds_all_new$OG %in% Combine1] <- "Positive_selection"
gene_dn_ds_all_new$Type[!(gene_dn_ds_all_new$OG %in% Combine1)] <- "Negative_selection"
# plot
ggplot(gene_dn_ds_all_new,aes(x=Type, y=dN_dS, fill=Type)) + geom_boxplot() +
  ylim(0,1) +
  #theme(panel.background = element_rect(fill = "white", colour = "black")) + # generate the whole border
  theme_bw() +                                       # control background and border
  theme(axis.line = element_line(colour = "black"),  # control background and border
        panel.grid.major = element_blank(),          # control background and border   
        panel.grid.minor = element_blank(),          # control background and border
        panel.border = element_blank(),              # control background and border
        panel.background = element_blank()) +        # control background and border
  
  theme(legend.position = c(0.2, 0.8)) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial"),
        legend.text = element_text(size = 13, family = "Arial")) +  ggtitle('') +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

t.test(gene_dn_ds_all_new$dN_dS[gene_dn_ds_all_new$Type=="Positive_selection"], gene_dn_ds_all_new$dN_dS[gene_dn_ds_all_new$Type=="Negative_selection"])










# classification of OGs with selected sites based on KO or function annotation.
og_panID_mapping <- read_tsv("data/representatives.tsv")

select_gene <- og_panID_mapping[og_panID_mapping$ortho_id %in% Combine1,]

# get the KO id
panID_KO <- read.table("kegg/ko_annotation_of_all_OG.txt", header = TRUE, sep = "\t")
select_gene$KO <- getSingleReactionFormula(panID_KO$ko, panID_KO$query, select_gene$representative)

# get the pathway
KO_pathway <- read.table("kegg/ko_pathway.txt", sep = "\t", stringsAsFactors = FALSE)
pathway <- read.table("kegg/pathway_list_kegg.txt", sep = "\t", stringsAsFactors = FALSE)
KO_pathway <- filter(KO_pathway,str_detect(KO_pathway$V1, "map"))
KO_pathway$V2 <- str_replace_all(KO_pathway$V2,"ko:","")

select_gene$pathwayID <- getMultipleReactionFormula(KO_pathway$V1,KO_pathway$V2, select_gene$KO)
select_gene_pathway0 <- select_gene[!is.na(select_gene$pathwayID),]
select_gene_pathway1 <- splitAndCombine(select_gene_pathway0$pathwayID,select_gene_pathway0$ortho_id,";")
select_gene_pathway1$pathway <- getSingleReactionFormula(pathway$V2, pathway$V1, select_gene_pathway1$v1)



select_gene_pathway1$pathway[select_gene_pathway1$pathway=="NA"] <- NA
select_gene_pathway1$pathway <- str_replace_all(select_gene_pathway1$pathway, " - fly", "") %>% str_replace_all(.," - other", "") %>%
  str_replace_all(.," - yeast", "") %>% str_replace_all(.," - Caulobacter", "") %>% str_replace_all(.," - animal", "")


# here we only choose the pathway with most genes
subsystem <- as.data.frame(table(select_gene_pathway1$pathway), stringsAsFactors = FALSE)
subsystem_filter <- subsystem[subsystem$Freq>=4,]
colnames(subsystem_filter) <- c("subsystem","num")

# set the factor level by the mean value of occur number
#Factor <-Result0 %>% group_by(pathway) %>% summarise(median=mean(Occur_num))
subsystem_filter <- subsystem_filter[order(subsystem_filter$num),]
subsystem_filter$subsystem <-factor(subsystem_filter$subsystem, levels=subsystem_filter$subsystem)
# remove two wrong
subsystem_filter0 <- filter(subsystem_filter, !str_detect(subsystem_filter$subsystem, "disease"))
# remove the metabolic pathways
subsystem_filter0 <- filter(subsystem_filter0, !str_detect(subsystem_filter0$subsystem, "Metabolic pathways"))
ggplot(subsystem_filter0, aes(x=subsystem, y=num)) +
  geom_bar(stat = "identity", fill = "#FF6666") +
  xlab('') + ylab('Number of genes') +
  theme_bw() +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=8, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  coord_flip()



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
og_panID_mapping0 <- og_panID_mapping[og_panID_mapping$ortho_id %in% Combine1, ]
# choose the OGs with sce genes
og_panID_sce <- og_panID_mapping0[str_detect(og_panID_mapping0$representative, "Saccharomyces_cerevisiae"),]
og_panID_sce$representative <- str_replace_all(og_panID_sce$representative,"Saccharomyces_cerevisiae@", "")
print(paste0(og_panID_sce$representative, collapse = ",")) # for DAVID