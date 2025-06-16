# Heatmap of E. coli and Associated Protists in Conventional Soil (-CHX)

library(pheatmap)
library(phyloseq)
library(tidyverse)

rm(list = ls())

# Set working directory and import 18S metadata and OTU table
setwd('C:/Users/User/OneDrive/Desktop/qian/R/CYC-with E.coil+unknow/18dna-cs')
metadata18 <- read.table("group(-cyc)1.txt", sep = "\t", header = TRUE)
otu18 <- read.table("level-8.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample")

# Merge and clean 18S table
mergedtable18 <- merge(otu18, metadata18, by = "sample") %>%
  dplyr::select(-c(sample, group)) %>%
  relocate(Replicate, .before = everything())

# Recode time labels
metadata18 <- metadata18 %>%
  mutate(time = case_when(
    grepl("CS-Day0\\(CHX-free\\)", group) ~ "Day0",
    grepl("CS-Day14\\(CHX-free\\)", group) ~ "Day14",
    grepl("CS-Day28\\(CHX-free\\)", group) ~ "Day28",
    TRUE ~ NA_character_
  )) %>%
  dplyr::select(-group)

# Import 16S data (only E. coli kept)
setwd('C:/Users/User/OneDrive/Desktop/qian/R/CYC-with E.coil+unknow/16dna-cs')
metadata16 <- read.table("group(-cyc)1.txt", sep = "\t", header = TRUE)
otu16 <- read.table("genus.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample")

# Merge 16S table (only E. coli)
mergedtable16 <- merge(otu16, metadata16, by = "sample") %>%
  dplyr::select(-c(sample, group)) %>%
  dplyr::select(Replicate, `g__Escherichia-Shigella`)

# Merge 16S and 18S tables
mergedtableall <- merge(mergedtable18, mergedtable16, by = "Replicate")
mergedtableall <- mergedtableall[order(mergedtableall$Replicate), ]
rownames(mergedtableall) <- mergedtableall$Replicate
mergedtableall <- mergedtableall[, -1]

# Z-score scaling and cleanup
mergedtableall <- as.data.frame(t(scale(mergedtableall)))
mergedtableall[] <- lapply(mergedtableall, as.numeric)
mergedtableall[is.na(mergedtableall)] <- 0

# Clean taxon names
mergedtableall <- rownames_to_column(mergedtableall, var = "Row_Name")
mergedtableall$Row_Name <- gsub("d8__", "", mergedtableall$Row_Name)
mergedtableall$Row_Name <- gsub("st__un_d[2-7]__", "un_", mergedtableall$Row_Name)
mergedtableall$Row_Name <- gsub("g__Escherichia-Shigella", "E.coli", mergedtableall$Row_Name)

# Import network file
setwd('C:/Users/User/OneDrive/Desktop/qian/R/network/18cs-cyc-ecoli')
network <- read.csv("ecoli_related_taxa_trophic.csv", sep = ",", header = TRUE, row.names = 1)

# Keep only Consumers and E.coli
network_filtered <- network[network$Trophic %in% c("Consumer", "E.coli"), ]
colnames(network_filtered)[colnames(network_filtered) == "Rank8"] <- "Row_Name"

# Subset merged data
mergedtableall <- mergedtableall %>%
  inner_join(network_filtered %>% dplyr::select(Row_Name), by = "Row_Name") %>%
  distinct()

mergedtableall <- column_to_rownames(mergedtableall, var = "Row_Name")

# Format metadata
metadata18$group <- factor(metadata18$time, levels = c("Day0", "Day14", "Day28"))
rownames(metadata18) <- metadata18$Replicate
metadata18 <- metadata18[, -3] %>% dplyr::select(group)
colnames(metadata18) <- c("Time")

# Separate E.coli and protists
group_ecoli <- mergedtableall[rownames(mergedtableall) == "E.coli", , drop = FALSE]
group_protists <- mergedtableall[rownames(mergedtableall) != "E.coli", , drop = FALSE]

# Cluster protists
clustered_protists <- group_protists[hclust(dist(group_protists))$order, ]
combined_data <- rbind(group_ecoli, clustered_protists)

# Annotate row groups
annotation_row <- data.frame(Taxon = rep(c("E.coli", "Bacterivorous protists"), c(1, nrow(clustered_protists))))
rownames(annotation_row) <- rownames(combined_data)

# Plot heatmap
heatmap_plot <- pheatmap(combined_data,
                         annotation_col = metadata18,
                         annotation_row = annotation_row,
                         cutree_cols = 3,
                         cellwidth = 20,
                         cellheight = 20,
                         cluster_rows = FALSE,
                         cluster_cols = FALSE,
                         gaps_col = c(5, 10),
                         gaps_row = 1)

# Save the heatmap
ggsave("18cs-cyc-heatmap-network-L8.jpg", plot = heatmap_plot, device = "jpeg", dpi = 600, width = 23, height = 14, units = "cm")
ggsave("18cs-cyc-heatmap-network-L8.tiff", plot = heatmap_plot, device = "tiff", dpi = 600, width = 23, height = 14, units = "cm")
