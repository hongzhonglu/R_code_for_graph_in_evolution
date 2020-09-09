# These scripts are used to compare the site wise evolution rate calculated using different methods.
# 2019.2.20
# Hongzhong Lu

#load package
#library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(ggExtra)


# function
plotXYdotGraph2 <- function (data_frame, paraX, paraY, xlab_name, ylab_name, title=''){
  # This function is used to plot the dot plot between two parameters
  # A bar column on the x and y axis will be ueed to do the statistical analysis
  # Input
  # data_frame: A dataframe
  # paraX: the column name of the parameter
  # paraY: the column name of another parameter
  
  p <- ggplot(data_frame, aes_string(x=paraX, y=paraY)) +
    geom_point(color="#69b3a2", alpha=0.8) +
    xlab(xlab_name) + 
    ylab(ylab_name) +
    theme(axis.text=element_text(size=20, family="Arial"),
          axis.title=element_text(size=24, family="Arial"),
          legend.text = element_text(size=20, family="Arial")) +
    ggtitle(title) +
    theme(panel.background = element_rect(fill = "white", color="black", size = 1),
          plot.margin = margin(1, 1, 1, 1, "cm")) 
  ggExtra::ggMarginal(p, type = "histogram", color="white")
  #+ggsave(out <- paste('result/','Metabolic gene number distribution for strains specific model from biocyc','.eps', sep = ""), width=5, height=5, dpi=300)
}
plotXYdotGraph <- function (data_frame, paraX, paraY, xlab_name, ylab_name, title=''){
  # This function is used to plot the dot plot between two parameters
  # Input
  # data_frame: A dataframe
  # paraX: the column name of the parameter
  # paraY: the column name of another parameter
  
  ggplot(data_frame, aes_string(x=paraX, y=paraY)) +
    geom_point() +
   # stat_smooth(method = "loess",
  #              col = "#C42126",
   #             se = FALSE,
    #            size = 1) +
    xlab(xlab_name) + 
    ylab(ylab_name) +
    theme(axis.text=element_text(size=20, family="Arial"),
          axis.title=element_text(size=24, family="Arial"),
          legend.text = element_text(size=20, family="Arial")) +
    ggtitle(title) +
    theme(panel.background = element_rect(fill = "white", color="black", size = 1),
          plot.margin = margin(1, 1, 1, 1, "cm")) 
  
  #+ggsave(out <- paste('result/','Metabolic gene number distribution for strains specific model from biocyc','.eps', sep = ""), width=5, height=5, dpi=300)
}

## part 1
## for the protein consercation score calculated using consurf
msa_aa_variety_percentage <- read_excel("data/msa_aa_variety_percentage_consurf.xlsx", sheet = "Sheet1") # This file has grades score of each site
score <- msa_aa_variety_percentage[,c('pos','ConSurf Grade')]
colnames(score) <- c('pos','consurf_grade')
# bar plot
hist(score$`consurf_grade`)
ggplot(score, aes(x=`consurf_grade`)) + geom_histogram(binwidth=0.5)
#  column annotation
#  Amino Acid Conservation Scores
#- POS: The position of the AA in the SEQRES derived sequence.
#- SEQ: The SEQRES derived sequence in one letter code.
#- SCORE: The normalized conservation scores.
#- COLOR: The color scale representing the conservation scores (9 - conserved, 1 - variable).
#- CONFIDENCE INTERVAL: When using the bayesian method for calculating rates, a confidence interval is assigned to each of the inferred evolutionary conservation scores.
#- CONFIDENCE INTERVAL COLORS: When using the bayesian method for calculating rates. The color scale representing the lower and upper bounds of the confidence interval.
#- B/E: Burried (b) or Exposed (e) residue.
#- FUNCTION: functional (f) or structural (s) residue (f - highly conserved and exposed, s - highly conserved and burried).
#- MSA DATA: The number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position.
#- RESIDUE VARIETY: The residues variety at each position of the multiple sequence alignment.
# *Below the confidence cut-off - The calculations for this site were performed on less than 6 non-gaped homologue sequences,					
# or the confidence interval for the estimated score is equal to- or larger than- 4 color grades.					
conservation_score <- read_excel("data/OG5327_conservation_score_consurf.xlsx") # This file has detailed score of each site
hist(conservation_score$SCORE)
ggplot(conservation_score, aes(x=SCORE)) + geom_histogram(binwidth=0.1) +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  xlab('Amino acid conservation score')



## part 2
## for the site-wise relative rates with HyPhy
relative_rate <- read_excel("data/OG5327_relative_rate_YPR001W.xlsx")
# bar plot
hist(relative_rate$Rate)
ggplot(relative_rate, aes(x=Rate)) + geom_histogram(binwidth=0.5)
# density plot
j <- 'Rate'
ggplot(relative_rate, aes_string(j)) +
  geom_density(fill="lightblue", alpha=1) +
  xlab("Site wise relative rate") + 
  ylab("Density") +
  theme(axis.text=element_text(size=20,face="bold", family="Arial"),
        axis.title=element_text(size=24,face="bold", family="Arial"),
        legend.text = element_text(size=20,face="bold", family="Arial")) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1))




## part 3
## for the site specific dN/dS
dNdS_YPR001W <- read_excel("data/datamonkey_OG5327_dNdS_YPR001W.xlsx")
#dNdS_YPR001W <- read_excel("data/extracted_OG5327_dNdS_YPR001W.xlsx")

j <- '`p-value`'

ggplot(dNdS_YPR001W, aes_string(j)) +
  geom_density(fill="lightblue", alpha=1) +
  geom_vline(xintercept=0.1, linetype="dotted", color = 'red', size=1) + 
  geom_text(x=0.20, y=5, label="P value = 0.1") + 
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial"),
        legend.text = element_text(size=20, family="Arial")) +
  xlab('P value')

# for omega
j <- 'omega'
dNdS_YPR001W$omega<- as.numeric(dNdS_YPR001W$omega)
dNdS_YPR001W0 <- filter(dNdS_YPR001W, omega < 10)
ggplot(dNdS_YPR001W0, aes_string(j)) +
  geom_density(fill="lightblue", alpha=1) +
  geom_vline(xintercept=1, linetype="dotted", color = 'red', size=1) + 
  geom_text(x=1, y=2, label="ω=1: neutral evolution") +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=24, family="Arial"),
        legend.text = element_text(size=20, family="Arial")) +
  xlab('dN/dS')

# for correlation between dS and dN
dNdS_YPR001W0 <- filter(dNdS_YPR001W, alpha < 10) %>%
  filter(.,beta < 10)
plotXYdotGraph2(data_frame=dNdS_YPR001W0, paraX='alpha', paraY='beta', xlab_name='dS', ylab_name='dN', title='')


## part 4
## for the conservation score calcualted by JDS
## part 3
## for the site specific dN/dS
conservation_jds <- read_excel("data/OG5327_conservation_score_YPR001W_jsd.xlsx")





## part 5
## how is going if we put different paramter together
score_summary <- conservation_score[, c('POS','SEQ','SCORE')]
score_summary$relative_rate <- relative_rate$Rate
score_summary$gap_ratio <- relative_rate$gap_ratio
score_summary$dnds <- dNdS_YPR001W$omega
score_summary$SCORE_jds <- as.numeric(conservation_jds$score)


# remove the unnormal dnds scores
score_summary_filter0 <- filter(score_summary, dnds < 100)

# filter based on the gap ratio
# the gap_ration cut-off is set at 0.8
# when further decrease the gap_ration. there is no obvious improvement in the correlation coefficient
# betweet two parameters
score_summary_filter <- filter(score_summary_filter0, gap_ratio <=0.8)



plotXYdotGraph(data_frame=score_summary_filter, paraX='SCORE', paraY='relative_rate', xlab_name='Conservation score', ylab_name='Site wise relative rate', title='')
# correlation test
res <- cor.test(score_summary_filter$SCORE, score_summary_filter$relative_rate, method = "pearson")
unlist(res['p.value'])
unlist(res$estimate)

plotXYdotGraph(data_frame=score_summary_filter, paraX='dnds', paraY='SCORE', xlab_name='dN/dS', ylab_name='Conservation score', title='')
# correlation test
res <- cor.test(score_summary_filter$SCORE, score_summary_filter$dnds, method = "pearson")
unlist(res['p.value'])
unlist(res$estimate)

plotXYdotGraph(data_frame=score_summary_filter, paraX='dnds', paraY='relative_rate', xlab_name='dN/dS', ylab_name='Site wise relative rate', title='')
# correlation test
res <- cor.test(score_summary_filter$relative_rate, score_summary_filter$dnds, method = "pearson")
unlist(res['p.value'])
unlist(res$estimate)


# to compare the protein residue conservation score using the JDS methods
score_summary_filter1 <- filter(score_summary_filter0, gap_ratio <=0.3)
score_summary_filter1$SCORE_jds <- - score_summary_filter1$SCORE_jds
plotXYdotGraph(data_frame=score_summary_filter1, paraX='SCORE_jds', paraY='SCORE', xlab_name='Conservation score JSD', ylab_name='Conservation score_consurf', title='')
res <- cor.test(score_summary_filter1$SCORE, score_summary_filter1$SCORE_jds, method = "pearson")
unlist(res['p.value'])
unlist(res$estimate)

plotXYdotGraph(data_frame=score_summary_filter1, paraX='SCORE_jds', paraY='relative_rate', xlab_name='Conservation score JSD', ylab_name='Site wise relative rate', title='')
res <- cor.test(score_summary_filter1$relative_rate, score_summary_filter1$SCORE_jds, method = "pearson")
unlist(res['p.value'])
unlist(res$estimate)


plotXYdotGraph(data_frame=score_summary_filter1, paraX='SCORE_jds', paraY='dnds', xlab_name='Conservation score JSD', ylab_name='dN/dS', title='')
res <- cor.test(score_summary_filter1$dnds, score_summary_filter1$SCORE_jds, method = "pearson")
unlist(res['p.value'])
unlist(res$estimate)











