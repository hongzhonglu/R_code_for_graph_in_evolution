# These scripts are used to produce supplementary figure 9 for the evolution project
# 2020.4.15
# Hongzhong Lu

#load package
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(hongR)



sce_gene_summary <- read_excel("data/sce_gene_summary.xlsx")
sce_gene_summary$sce_gene <- str_replace_all(sce_gene_summary$sce_gene,"\\[", "") %>% str_replace_all(.,"Saccharomyces_cerevisiae@","") %>%
  str_replace_all(.,"\\]", "") %>% str_replace_all(.,"'", "")

# input more parameters related to ortholog
ortholog_occurence <- read_tsv("data/ortholog_occurence_num_all.tsv")
sce_gene_summary0 <- left_join(sce_gene_summary,ortholog_occurence,by=c("OrthologID"="ID") )



################################################################################
# WGD gene analysis
# add the whole genome duplication information for sce gene
gene_WGD <- read_excel("data/yeast_gene_annotation_SGD.xlsx", sheet = "gene_WGD")
sce_gene_summary0$WGD_gene_num <- NA


for (i in 1: nrow(sce_gene_summary0)){
  print(i)
  gene0 <- sce_gene_summary0$sce_gene[i]
  gene1 <- unlist(str_split(gene0,",")) %>% str_trim(., side = "both")
  exist0 <- gene1 %in% gene_WGD$locus
  num_WGD <- sum(as.numeric(exist0))
  sce_gene_summary0$WGD_gene_num[i] <- num_WGD
  
}

sce_gene_summary0$WGD[sce_gene_summary0$WGD_gene_num > 0] <- "WGD"
sce_gene_summary0$WGD[sce_gene_summary0$WGD_gene_num == 0] <- "non_WGD"


# plot
sce_gene_summary0 %>%
  ggplot(aes(x=WGD,y=species_num, fill=WGD, alpha=0.5)) +
  geom_boxplot(alpha = 0.1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1))

ggplot(sce_gene_summary0, aes(species_num, fill = WGD, colour = WGD)) +
  geom_density(alpha = 0.1) +
  theme(legend.position = c(0.2, 0.8)) +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1))
# statistical analysis
# as the whole, the non_WGD gene is more conserved across different species
g1 <- sce_gene_summary0$species_num[sce_gene_summary0$WGD=="non_WGD"]
g2 <- sce_gene_summary0$species_num[sce_gene_summary0$WGD=="WGD"]

t.test(g1,g2)




####################################################################################
# gene function annotate based on yeast GEM
core_metabolic_gene <- read_excel("data/yeastGEM_october.xls",  sheet = "core_rxn")
core_metabolic_gene <- core_metabolic_gene[!is.na(core_metabolic_gene$GPR), ]


###### function to splite all the related genes
getSpliteGene <- function(genelist, reactionName){
  ##split of or relation
  #genelist <- yeast_7_7$Gene.reaction.association
  #reactionName <-yeast_7_7$Rxn.name
  
  rxnNum <- length(reactionName)
  
  
  GR <- list() # gene relation
  ss <- list()
  
  
  Add_and <- function(x1){
    tt0 <- vector()
    tt <- x1
    tt_length <- length(tt)
    for (i in 1:tt_length){
      tt0[i] <- paste(i, tt[i], sep = ";")
    }
    return(tt0)
  }
  
  
  Add_or <- function(x1){ #x1 <- GR[[1]]
    tt0 <- vector()
    tt <- x1
    tt_length <- length(tt)
    for (i in 1:tt_length){
      tt0[i] <- paste(i, tt[i], sep = ";")
    }
    return(tt0)
  }
  
  
  ## split the first order of gene relation
  
  splitOrRelation <- function(genelist, reactionName){
    #genelist <- "((YAL023C and YDL095W) or YDL093W or YJR143C or YOR321W)"
    #reactionName <- "r_0016"
    GR <- list() # gene relation
    ss <- list()
    if (!(str_detect(genelist, "\\) or \\(") | str_detect(genelist, "\\) or") | str_detect(genelist, "or \\(") | str_detect(genelist, "\\) and \\(") | str_detect(genelist, "\\) and") | str_detect(genelist, "and \\("))){
      
      GR[1] <- genelist
      ss[[1]] <- paste(reactionName,GR, sep = "@none;" )
      
    } else {
      
      if (str_detect(genelist, "\\) or \\(")){
        GR[1] <- str_split(genelist,"\\) or \\(" )
        
      } else{
        GR[1] <- genelist
      }
      
      if (str_detect(paste0(unlist(GR[1]),collapse = "@@"), "\\) or")) {
        GR[1] <- str_split(paste0(unlist(GR[1]),collapse = "@@"),"\\) or")
      } else{
        GR[1] <- GR[1]
      }
      
      if (str_detect(paste0(unlist(GR[1]),collapse = "@@"), "or \\(")){
        GR[1] <- str_split(paste0(unlist(GR[1]),collapse = "@@"),"or \\(" )
        
      } else{
        GR[1] <- GR[1]
      }
      
      
      
      if (str_detect(genelist, "\\) and \\(")){
        GR[1] <- str_split(genelist,"\\) and \\(" )
        
      } else{
        GR[1] <- GR[1]
      }
      
      if (str_detect(paste0(unlist(GR[1]),collapse = "@@"), "\\) and")) {
        GR[1] <- str_split(paste0(unlist(GR[1]),collapse = "@@"),"\\) and")
      } else{
        GR[1] <- GR[1]
      }
      
      if (str_detect(paste0(unlist(GR[1]),collapse = "@@"), "and \\(")){
        GR[1] <- str_split(paste0(unlist(GR[1]),collapse = "@@"),"and \\(" )
        
      } else{
        GR[1] <- GR[1]
      }
      
      
      
      ##for genelist <- "((YAL023C and YDL095W) or YDL093W or YJR143C or YOR321W)" 
      ##reactionName <- "r_0016"
      ##tt <- GR[1]
      
      sl <-length(GR[[1]])
      for (i in 1:sl){
        if( (str_detect(GR[[1]][i],"\\(\\(")) | (str_detect(GR[[1]][i],"\\)\\)"))){
          GR[[1]][i] <- GR[[1]][i]
        } else if (str_detect(GR[[1]][i], "or")){
          GR[[1]][i]  <- str_split(GR[[1]][i], "or")
        } 
      }
      ##### above is for special condition
      
      
      
      GR[1] <- paste0(unlist(GR[1]),collapse = "@@")
      GR[1] <- str_split(GR[1],"@@")
      
      if (str_detect(genelist, "\\) or \\(") | str_detect(genelist, "\\) or") | str_detect(genelist, "or \\(")){
        ss[[1]] <- paste(reactionName,Add_or(GR[[1]]), sep = "@or" )
      } else {
        ss[[1]] <- paste(reactionName,Add_and(GR[[1]]), sep = "@and" )
      }
    }
    
    return(ss[[1]])
    
  }
  
  for (i in 1:rxnNum){
    
    ss[[i]] <- splitOrRelation(genelist[i],reactionName[i])
    
  }
  
  
  GR1 <- unlist(ss)
  
  GR2 <-  str_replace_all(GR1,"\\(", "")
  GR2 <-  str_replace_all(GR2,"\\)", "")
  GR3 <- str_split(GR2,";")
  tt0 <- length(GR3)
  GR4 <- data.frame( ID= character(tt0),gene=character(tt0), stringsAsFactors = FALSE )
  for (i in 1:tt0){
    GR4$ID[i] <- GR3[[i]][1]
    GR4$gene[i] <- GR3[[i]][2]
    
  }
  
  ## split the first order of gene relation
  rxnNum2 <- length(GR4$ID)
  GR_2 <- list()
  ss2 <- list()
  for (i in 1:rxnNum2){
    if (str_detect(GR4$gene[i], " or ")){
      GR_2[i] <- str_split(GR4$gene[i]," or " )
      ss2[[i]] <- paste(GR4$ID[i],Add_or(GR_2[[i]]), sep = "@or" )
      
    } else if(str_detect(GR4$gene[i], " and ")){
      GR_2[i] <- str_split(GR4$gene[i]," and " )
      ss2[[i]] <- paste(GR4$ID[i],Add_and(GR_2[[i]]), sep = "@and" )
      
    } else{
      GR_2[i] <- GR4$gene[i]
      ss2[[i]] <- paste(GR4$ID[i],GR_2[[i]], sep = "@none;" )
      
    }
  }
  
  
  ## obtain the final formula
  GR_3 <- unlist(ss2)
  GR_4 <- str_split(GR_3,";")
  tt <- length(GR_4)
  GR_5 <- data.frame( ID= character(tt),gene=character(tt), stringsAsFactors = FALSE )
  
  for (i in 1:tt){
    GR_5$ID[i] <- GR_4[[i]][1]
    GR_5$gene[i] <- GR_4[[i]][2]
  }
  
  GR_or_and <- str_split(GR_5$ID,"@")
  
  for (i in 1:tt){
    GR_5$IDnew[i] <- GR_or_and[[i]][1]
    GR_5$R1[i] <- GR_or_and[[i]][2]
    GR_5$R2[i] <- GR_or_and[[i]][3]
  }
  
  GPs_redesign <- select(GR_5,IDnew, R1, R2, gene) ##obtain the new format of GPRs, the can be base to add the new GPRs or correct GPRs
  
  return(GPs_redesign)
  
}

GPs_redesign_yeast <- getSpliteGene(core_metabolic_gene$GPR, core_metabolic_gene$Abbreviation)
GPs_redesign_yeast$Description <- getSingleReactionFormula(core_metabolic_gene$Description, core_metabolic_gene$Abbreviation,GPs_redesign_yeast$IDnew)
GPs_redesign_yeast$Subsystem_new<- getSingleReactionFormula(core_metabolic_gene$Subsystem_new, core_metabolic_gene$Abbreviation,GPs_redesign_yeast$IDnew)
GPs_redesign_yeast$gene <- str_trim(GPs_redesign_yeast$gene, side = "both")
# save the core gene list
core_gene_df <- data.frame(core_gene=unique(GPs_redesign_yeast$gene), stringsAsFactors = FALSE)
write.table(core_gene_df, "result/core_metabolic_gene.txt", row.names = FALSE, sep = "\t")


# then based on ortholog result we will calculate the gene copy of obove genes from core metabolic pathways
all_file <- list.files("ortholog_find_result/")
for (i in seq_along(all_file)) {
  print(i)
  in_dir <- paste("ortholog_find_result/", all_file[i], sep = "")
  input <- read_delim(in_dir, "\t", escape_double = FALSE, trim_ws = TRUE)
  input <- input[!is.na(input$Saccharomyces_cerevisiae), ]
  col0 <- colnames(input)
  col0 <- col0[col0 != "Saccharomyces_cerevisiae"]
  col0 <- col0[col0 != "...1"]
  col0 <- col0[col0 != "X1"]
  input[["Saccharomyces_cerevisiae"]] <- str_replace_all(input[["Saccharomyces_cerevisiae"]], "\\,", "@@")
  # then we build a mapping between sce gene and other genes
  input_new <- splitAndCombine(input[["Saccharomyces_cerevisiae"]], input[[col0[1]]], sep0 = "@@")
  GPs_redesign_yeast[[col0[1]]] <- getSingleReactionFormula(input_new$v2, input_new$v1, GPs_redesign_yeast$gene)
}

# Combine the isoenzyme as one group
# Here with the glucose as an example

#complex
#test <- GPs_redesign_yeast[GPs_redesign_yeast$IDnew=="r_0886", 1:ncol(GPs_redesign_yeast)]
#test <- test[1,] #r_0886 YMR205C
#test <- test[2,] #r_0886 YGR240C

#isoenzyme
test <- GPs_redesign_yeast[GPs_redesign_yeast$IDnew=="r_1166", 1:ncol(GPs_redesign_yeast)]
#r_1166 glucose transporter reaction
#r_0450 FBA1
#r_0486 TDH
#r_0892 PGK1
#r_0356 GPM1
#r_0366	enolase ENO1
#r_0962	pyruvate kinase PYK1 
#r_0958	pyruvate carboxylase PYC1
#r_2115 alcohol dehydrogenase, (acetaldehyde to ethanol)
#r_0163 alcohol dehydrogenase, (ethanol to acetaldehyde)

#TCA: r_0300	citrate synthase

#PPP: r_0091
#PPP: r_0888
#PPP: r_0889
#PPP: r_0907
#PPP: r_1048
#PPP: r_1049
#PPP: r_1050

# if we assume the gene annotation in yeast GEM is right, then the gene from other genomes should also
# have the similar function
species <- colnames(test)[7:ncol(test)]
number_paralog <- c()
for (i in seq_along(species)){
s0 <- species[i]
combine0 <- test[[s0]]
combine0 <- paste0(combine0[combine0!="NA"], collapse = ",")
combine1 <- unlist(str_split(combine0, ",")) %>% unique(str_trim(.,side = "both"))
num0 <- length(combine1)
number_paralog <- c(number_paralog, num0)
}

df0 <- data.frame(species=species, stringsAsFactors = FALSE)
df0[, "protein_homolog_number"] <- number_paralog
plot(density(df0[, "protein_homolog_number"]))

df0[nrow(df0) + 1,] = c("Saccharomyces_cerevisiae", length(unique(test$gene))-1)
df0$protein_homolog_number <- as.numeric(df0$protein_homolog_number)

# connect the glucose transportor with crabtree effect
yeast_species_classification <- read_csv("data/yeast_species_classification.csv")
yeast_species_classification0 <- yeast_species_classification
yeast_species_classification0$protein_homolog_number <- getSingleReactionFormula(df0$protein_homolog_number,df0$species,yeast_species_classification0$old_species_id)
yeast_species_classification0$protein_homolog_number <- as.numeric(yeast_species_classification0$protein_homolog_number)


# explore the relation between the duplication of glucose transportor and Crabtree effect
yeast_species_classification1 <- yeast_species_classification0[!is.na(yeast_species_classification0$crabtree_effect), ]
yeast_species_classification1$crabtree_effect <- as.factor(yeast_species_classification1$crabtree_effect)
# plot
yeast_species_classification1 %>%
  ggplot(aes(x=crabtree_effect,y=protein_homolog_number, fill=crabtree_effect)) +
  stat_boxplot(geom ='errorbar', width = 0.25) + # add caps
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1))

g1 <- yeast_species_classification1$protein_homolog_number[yeast_species_classification1$crabtree_effect=="Yes"]
g2 <- yeast_species_classification1$protein_homolog_number[yeast_species_classification1$crabtree_effect=="No"]
t.test(g1,g2)

# wilcon.test
wilcox.test(g1,g2, alternative = "two.sided")



# plot2
# explore the relation between the duplication of glucose transportor and whole genome duplication
yeast_species_classification0 %>%
  ggplot(aes(x=WGD,y=protein_homolog_number, fill=WGD, alpha=0.5)) +
  geom_boxplot(alpha = 0.1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text=element_text(size=12, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') +
  theme(legend.position = "none") +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1))

w1 <- yeast_species_classification0$protein_homolog_number[yeast_species_classification0$WGD=="Yes"]
w2 <- yeast_species_classification0$protein_homolog_number[yeast_species_classification0$WGD=="No"]
t.test(w1,w2)


