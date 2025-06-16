library(pheatmap)
library(phyloseq)
library(tidyverse)

rm(list = ls())

# Set working directory
setwd("C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os")

# Load metadata for 18S data (CHX-free groups)
metadata18 <- read.table("group(-cyc)1.txt", sep = "\t", header = TRUE)

# Load OTU table for 18S (level-8) and transpose
otu18 <- read.table("level-8.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("sample")

# Merge OTU with metadata
mergedtable18 <- merge(otu18, metadata18, by = "sample") %>%
  dplyr::select(-sample, -group) %>%
  relocate(Replicate, .before = everything())

# Reformat metadata18 with time point information
metadata18 <- metadata18 %>%
  mutate(time = case_when(
    grepl("OS-Day0\\(CHX-free\\)", group) ~ "Day0",
    grepl("OS-Day14\\(CHX-free\\)", group) ~ "Day14",
    grepl("OS-Day28\\(CHX-free\\)", group) ~ "Day28",
    TRUE ~ NA_character_
  )) %>%
  dplyr::select(-group)

# Load 16S metadata and OTU table (focusing on Escherichia-Shigella)
setwd("C:/Users/User/OneDrive/Desktop/qian/R/CYC-with E.coil+unknow/16dna-os")
metadata16 <- read.table("group(-cyc)1.txt", sep = "\t", header = TRUE)
otu16 <- read.table("genus.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("sample")
mergedtable16 <- merge(otu16, metadata16, by = "sample") %>%
  dplyr::select(-sample, -group) %>%
  dplyr::select(Replicate, `g__Escherichia-Shigella`)

# Merge 16S and 18S tables
mergedtableall <- merge(mergedtable18, mergedtable16, by = "Replicate")
mergedtableall <- mergedtableall[order(mergedtableall$Replicate), ]
rownames(mergedtableall) <- mergedtableall$Replicate
mergedtableall <- mergedtableall[, -1] %>%
  t() %>% scale() %>% t() %>% as.data.frame()

mergedtableall[] <- lapply(mergedtableall, as.numeric)
mergedtableall[is.na(mergedtableall)] <- 0

# Format taxon names
mergedtableall <- rownames_to_column(mergedtableall, var = "Row_Name") %>%
  mutate(Row_Name = str_replace_all(Row_Name, c(
    "d8__" = "",
    "st__un_d[2-7]__" = "un_",
    "g__Escherichia-Shigella" = "E.coli"
  )))

# Load network-filtered taxa
setwd("C:/Users/User/OneDrive/Desktop/qian/R/network/18os-cyc-ecoli")
network <- read.csv("ecoli_related_taxa_trophic.csv", sep = ",", header = TRUE, row.names = 1)
network_filtered <- network[network$Trophic %in% c("Consumer", "E.coli"), ]
colnames(network_filtered)[colnames(network_filtered) == "Rank8"] <- "Row_Name"

# Join network-filtered taxa with scaled abundance data
mergedtableall <- mergedtableall %>%
  inner_join(network_filtered %>% dplyr::select(Row_Name), by = "Row_Name") %>%
  distinct() %>%
  column_to_rownames(var = "Row_Name")

# Format metadata
metadata18 <- metadata18 %>%
  column_to_rownames(var = "Replicate") %>%
  dplyr::select(group = time) %>%
  mutate(group = factor(group, levels = c("Day0", "Day14", "Day28")))
colnames(metadata18) <- "Time"

# Cluster protist rows but keep E.coli on top
group_ecoli <- mergedtableall[rownames(mergedtableall) == "E.coli", , drop = FALSE]
group_protists <- mergedtableall[rownames(mergedtableall) != "E.coli", , drop = FALSE]
clustered_protists <- group_protists[hclust(dist(group_protists))$order, ]
combined_data <- rbind(group_ecoli, clustered_protists)

# Annotate rows
annotation_row <- data.frame(Taxon = rep(c("E.coli", "Bacterivorous protists"), c(1, nrow(clustered_protists))))
rownames(annotation_row) <- rownames(combined_data)

# Plot heatmap
heatmap_plot <- pheatmap(
  combined_data,
  annotation_col = metadata18,
  annotation_row = annotation_row,
  cutree_cols = 3,
  cellwidth = 20,
  cellheight = 20,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = c(5, 10),
  gaps_row = 1
)

# Save heatmap in high resolution
output_paths <- ".../18os-cyc-heatmap-network-L8"


for (path in output_paths) {
  ggsave(paste0(path, ".jpg"), plot = heatmap_plot, device = "jpeg", dpi = 600, width = 23, height = 10, units = "cm")
  ggsave(paste0(path, ".tiff"), plot = heatmap_plot, device = "tiff", dpi = 600, width = 23, height = 10, units = "cm")
}
