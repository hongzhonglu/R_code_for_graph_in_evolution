# This script is used to explore the distribution of pan-gene occurance
library(ggplot2)
library(readr)
library(dplyr)
library(stringr)
library(hongR)


ortholog <- read_tsv("data/ortholog_occurence_num_all.tsv")
ortholog$max_seq_ratio <- ortholog$max_seq_nums_from_single_tax / ortholog$protein_num

# further to explore the relation between the species number and the coveraged total ortholog number
total_ortho <- vector()
species_num <- c(1, 2, 4, 6, 8, 10, 20, 40, 100, 200, 300, 342)
for (i in species_num) {
  print(i)
  s1 <- filter(ortholog, species_num >= i)
  num <- length(s1$X1)
  total_ortho <- c(total_ortho, num)
}
plot(species_num, total_ortho)


total_ortho <- vector()
species_num <- 1:30
for (i in species_num) {
  print(i)
  s1 <- filter(ortholog, species_num >= i)
  num <- length(s1$X1)
  total_ortho <- c(total_ortho, num)
}
plot(species_num, total_ortho)





# filter1
ortholog1 <- filter(ortholog, species_num > 6)
# plot the column density plot
ggplot(ortholog1, aes(species_num)) +
  geom_histogram(aes(y = ..density..), alpha = 0.7, fill = "#333333") +
  theme(panel.background = element_rect(fill = "#ffffff"))


ggplot(ortholog1, aes(protein_num)) +
  geom_histogram(aes(y = ..density..), alpha = 0.7, fill = "#333333") +
  theme(panel.background = element_rect(fill = "#ffffff"))


ggplot(ortholog1, aes(average_duplicate)) +
  geom_histogram(aes(y = ..density..), alpha = 0.7, fill = "#333333") +
  theme(panel.background = element_rect(fill = "#ffffff"))

# plot the density plot
ggplot(ortholog1, aes(x = species_num)) +
  geom_density()


# plot the accumulative plot
ggplot(ortholog1, aes(species_num)) +
  stat_ecdf() #+
# geom_vline(xintercept = 4)



# filter2
ortholog2 <- filter(ortholog1, average_duplicate <= 1.5)
ggplot(ortholog2, aes(max_seq_nums_from_single_tax)) +
  geom_histogram(aes(y = ..density..), alpha = 0.7, fill = "#333333") +
  theme(panel.background = element_rect(fill = "#ffffff")) +
  xlim(0, 10)


# filter3
ortholog3 <- filter(ortholog2, max_seq_nums_from_single_tax <= 3)

ggplot(ortholog3, aes(species_num)) +
  geom_histogram(aes(y = ..density..), alpha = 0.7, fill = "#333333") +
  theme(panel.background = element_rect(fill = "#ffffff")) +
  xlim(5, 20)




####### OG duplication analysis ##############################
# extract the OGs with sce
# Analysis one-build the conservation and dn_ds
ortholog <- read_tsv("data/ortholog_occurence_num_all.tsv")
ortholog$seq_species_ratio <- ortholog$protein_num / ortholog$species_num
# filter one
ortholog_special <- filter(ortholog, max_seq_nums_from_single_tax >= 4 & max_seq_nums_from_single_tax <= 5)

# filter two
ortholog_special1 <- filter(ortholog_special, species_num >= 7 & seq_species_ratio <= 1.5)

# here if we increase the max duplication number from 3 to 5 while other parameters are the same, we can increase
# 284 proteins to analysis
# however, if we further increase the max duplication number, the new protein clusters number will be not increased
# obviously accordingly.

# save the result for further analysis
write.table(ortholog_special1, "result/ortholog_special_second_batch.txt", row.names = FALSE, sep = "\t")

