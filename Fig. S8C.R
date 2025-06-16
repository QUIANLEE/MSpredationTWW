library(pheatmap)
library(phyloseq)
library(tidyverse)

rm(list = ls())

# Set working directory for 18S data with CHX
setwd('C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os')

# Load metadata and OTU table
metadata18 <- read.table("group(+cyc)1.txt", sep = "\t", header = TRUE)
otu18 <- read.table("level-8.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>%
  t() %>% as.data.frame() %>% rownames_to_column("sample")

# Merge metadata and OTU table
mergedtable18 <- merge(otu18, metadata18, by = "sample") %>%
  dplyr::select(-c(sample, group)) %>%
  relocate(Replicate, .before = everything())

# Add time labels for CHX-treated groups
metadata18 <- metadata18 %>%
  mutate(time = case_when(
    grepl("OS-Day0\\(CHX-free\\)", group) ~ "Day0",
    grepl("OS-Day14\\(with CHX\\)", group) ~ "Day14(CHX)",
    grepl("OS-Day28\\(with CHX\\)", group) ~ "Day28(CHX)",
    TRUE ~ NA_character_
  )) %>%
  dplyr::select(-group)

# Set working directory for 16S data
setwd('C:/Users/User/OneDrive/Desktop/qian/R/CYC-with E.coil+unknow/16dna-os')

# Load 16S metadata and OTU table
metadata16 <- read.table("group(+cyc)1.txt", sep = "\t", header = TRUE)
otu16 <- read.table("genus.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>%
  t() %>% as.data.frame() %>% rownames_to_column("sample")

# Extract Escherichia-Shigella data and merge
mergedtable16 <- merge(otu16, metadata16, by = "sample") %>%
  dplyr::select(-c(sample, group)) %>%
  dplyr::select(Replicate, `g__Escherichia-Shigella`)

# Merge 16S and 18S data
mergedtableall <- merge(mergedtable18, mergedtable16, by = "Replicate")
mergedtableall <- mergedtableall[order(mergedtableall$Replicate), ]
rownames(mergedtableall) <- mergedtableall$Replicate
mergedtableall <- mergedtableall[, -1]

# Normalize, scale, and clean the data
mergedtableall <- as.data.frame(t(scale(mergedtableall)))
mergedtableall[] <- lapply(mergedtableall, as.numeric)
mergedtableall[is.na(mergedtableall)] <- 0

# Clean row names
mergedtableall <- rownames_to_column(mergedtableall, var = "Row_Name")
mergedtableall$Row_Name <- gsub("d8__", "", mergedtableall$Row_Name)
mergedtableall$Row_Name <- gsub("st__un_d[2-7]__", "un_", mergedtableall$Row_Name)
mergedtableall$Row_Name <- gsub("g__Escherichia-Shigella", "E.coli", mergedtableall$Row_Name)

# Set working directory for network data
setwd('C:/Users/User/OneDrive/Desktop/qian/R/network/18os+cyc-ecoli')
network <- read.csv("ecoli_related_taxa_trophic.csv", sep = ",", header = TRUE, row.names = 1)

# Filter network data to include only E.coli and Consumers
network_filtered <- network[network$Trophic %in% c("Consumer", "E.coli"), ]
colnames(network_filtered)[colnames(network_filtered) == "Rank8"] <- "Row_Name"

# Keep only taxa in the filtered network
mergedtableall <- mergedtableall %>%
  inner_join(network_filtered %>% dplyr::select(Row_Name), by = "Row_Name") %>%
  distinct()
mergedtableall <- column_to_rownames(mergedtableall, var = "Row_Name")

# Format metadata
metadata18$group <- factor(metadata18$time, levels = c("Day0", "Day14(CHX)", "Day28(CHX)"))
rownames(metadata18) <- metadata18$Replicate
metadata18 <- metadata18[, -3] %>% dplyr::select(group)
colnames(metadata18) <- c("Time")

# Separate E.coli and protist groups
group_ecoli <- mergedtableall[rownames(mergedtableall) == "E.coli", , drop = FALSE]
group_protists <- mergedtableall[rownames(mergedtableall) != "E.coli", , drop = FALSE]

# Cluster protist group
clustered_protists <- group_protists[hclust(dist(group_protists))$order, ]
combined_data <- rbind(group_ecoli, clustered_protists)

# Annotate heatmap rows
annotation_row <- data.frame(Taxon = rep(c("E.coli", "Bacterivorous protists"), c(1, nrow(clustered_protists))))
rownames(annotation_row) <- rownames(combined_data)

# Plot heatmap
a <- pheatmap(combined_data,
              annotation_col = metadata18,
              annotation_row = annotation_row,
              cutree_cols = 3,
              cellwidth = 20,
              cellheight = 20,
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              gaps_col = c(5, 10),
              gaps_row = 1)

# Save plots
ggsave("18os+cyc-heatmap-network-L8.jpg", plot = a, device = "jpeg", dpi = 600, width = 21, height = 14, units = "cm")
ggsave("18os+cyc-heatmap-network-L8.tiff", plot = a, device = "tiff", dpi = 600, width = 21, height = 14, units = "cm")
