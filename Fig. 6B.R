# Load necessary libraries
library(dplyr)
library(tibble)
library(MASS)
library(reshape2)
library(Hmisc)
library(igraph)
library(ggplot2)
library(ggpubr)
library(ggsignif)

# ---------------------------
# Data Preparation: 18S (protists)
# ---------------------------
rm(list = ls())
setwd("C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os")

metadata18 <- read.table("group(+cyc).txt", sep = "\t", header = TRUE)
otu18 <- read.table("level-8.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>%
  t() %>% as.data.frame() %>% rownames_to_column("sample")

merged18 <- merge(otu18, metadata18, by = "sample") %>%
  select(-sample, -group) %>%
  relocate(Replicate, .before = everything()) %>%
  t() %>% as.data.frame()

colnames(merged18) <- as.character(merged18[1, ])
merged18 <- merged18[-1, ] %>% mutate_all(as.numeric)

merged18_pseudo <- merged18 + 1
clr18 <- microbiome::transform(merged18_pseudo, "clr") %>% as.data.frame()

# ---------------------------
# Data Preparation: 16S (E. coli only)
# ---------------------------
setwd("C:/Users/User/OneDrive/Desktop/qian/R/CYC-with E.coil+unknow/16dna-os")

metadata16 <- read.table("group(+cyc).txt", sep = "\t", header = TRUE)
otu16 <- read.table("genus.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>%
  t() %>% as.data.frame() %>% rownames_to_column("sample")

merged16 <- merge(otu16, metadata16, by = "sample") %>%
  select(-sample) %>%
  relocate(Replicate, .before = everything()) %>%
  t() %>% as.data.frame()

colnames(merged16) <- as.character(merged16[1, ])
merged16 <- merged16[-1, ] %>% mutate_all(as.numeric)

merged16_pseudo <- merged16 + 1
clr16 <- microbiome::transform(merged16_pseudo, "clr") %>% as.data.frame()
clr16 <- t(clr16) %>% as.data.frame() %>% select("g__Escherichia-Shigella") %>% t() %>% as.data.frame()

# ---------------------------
# Standardize Taxa Names
# ---------------------------
clr18 <- rownames_to_column(clr18, var = "Row_Name")
clr16 <- rownames_to_column(clr16, var = "Row_Name")
clr18$Row_Name <- gsub("d8__", "", clr18$Row_Name)
clr18$Row_Name <- gsub("st__un_d[2-7]__", "un_", clr18$Row_Name)
clr16$Row_Name <- gsub("g__Escherichia-Shigella", "E.coli", clr16$Row_Name)

clr18 <- column_to_rownames(clr18, var = "Row_Name")
clr16 <- column_to_rownames(clr16, var = "Row_Name")

clr18 <- clr18 %>% t() %>% as.data.frame() %>% rownames_to_column("sample")
clr16 <- clr16 %>% t() %>% as.data.frame() %>% rownames_to_column("sample")

# ---------------------------
# Merge CLR-transformed 18S and 16S Data
# ---------------------------
merged_all <- merge(clr18, clr16, by = "sample")
merged_all <- t(merged_all) %>% as.data.frame()
colnames(merged_all) <- as.character(merged_all[1, ])
merged_all <- merged_all[-1, ] %>% mutate_all(as.numeric)

# ---------------------------
# Filter Taxa by Prevalence (>=10%)
# ---------------------------
dt <- merged_all
dt[dt > 1] <- 1
valid_taxa <- which(rowSums(dt) / ncol(dt) > 0.1)
merged_all <- merged_all[valid_taxa, ] %>% as.matrix()

# ---------------------------
# Spearman Correlation
# ---------------------------
cor_result <- rcorr(t(merged_all), type = "spearman")
r_mat <- cor_result$r
p_mat <- cor_result$P

write.csv(r_mat, "C:/Users/User/OneDrive/Desktop/qian/R/network/18os+cyc-ecoli/rvaluecor_level8.csv")
write.csv(p_mat, "C:/Users/User/OneDrive/Desktop/qian/R/network/18os+cyc-ecoli/pvaluecor_level8.csv")

r_mat[p_mat >= 0.05 | abs(r_mat) < 0.6] <- 0
p_mat[p_mat >= 0.05 | abs(r_mat) < 0.6] <- 0
p_mat[is.na(p_mat)] <- 1

write.csv(r_mat, "C:/Users/User/OneDrive/Desktop/qian/R/network/18os+cyc-ecoli/Rvalue_filter_level8.csv")
write.csv(p_mat, "C:/Users/User/OneDrive/Desktop/qian/R/network/18os+cyc-ecoli/Pvalue_filter_level8.csv")

# ---------------------------
# Extract Significant E. coli Correlations
# ---------------------------
r_ecoli <- as.data.frame(r_mat) %>% filter(abs(E.coli) > 0) %>% t() %>% as.data.frame() %>% filter(abs(E.coli) > 0)
p_ecoli <- as.data.frame(p_mat) %>% filter(abs(E.coli) != 0) %>% t() %>% as.data.frame() %>% filter(abs(E.coli) != 0)

diag(r_ecoli) <- 0
diag(p_ecoli) <- NA

write.csv(r_ecoli, "C:/Users/User/OneDrive/Desktop/qian/R/network/18os+cyc-ecoli/R.ecoli_level8_sig.csv")
write.csv(p_ecoli, "C:/Users/User/OneDrive/Desktop/qian/R/network/18os+cyc-ecoli/P.ecoli_level8_sig.csv")

# ---------------------------
# Filter Abundance Matrix by Significant Correlations
# ---------------------------
r_ecoli <- read.csv("C:/Users/User/OneDrive/Desktop/qian/R/network/18os+cyc-ecoli/R.ecoli_level8_sig.csv", row.names = 1)
merged_all <- merged_all %>% as.data.frame() %>% mutate(X = row.names(.))
merged_all <- merged_all[merged_all$X %in% rownames(r_ecoli), ] %>% select(-X) %>% mutate_all(as.numeric)
merged_all <- as.matrix(merged_all)

r_ecoli <- as.matrix(r_ecoli)

# ---------------------------
# Network Construction
# ---------------------------
graph_net <- graph_from_adjacency_matrix(r_ecoli, weight = TRUE, mode = "undirected")
graph_net <- simplify(graph_net)

# Assign edge widths
E(graph_net)$width <- abs(E(graph_net)$weight) * 5

# Node abundance and labeling
node_abundance <- rowSums(merged_all)
node_log10 <- log10(node_abundance)

V(graph_net)$Abundance <- node_log10[V(graph_net)$name]

# Assign taxonomy groups as node color
tax_annot <- read.table("C:/Users/User/OneDrive/Desktop/qian/R/network/18os+cyc-ecoli/level2-8.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)
tax_annot <- rownames_to_column(tax_annot, var = "Row_Name")
tax_annot <- semi_join(tax_annot, data.frame(Row_Name = rownames(as.data.frame(node_abundance))), by = "Row_Name")
tax_annot <- column_to_rownames(tax_annot, var = "Row_Name")

V(graph_net)$G <- as.character(tax_annot[V(graph_net)$name, ])
V(graph_net)$label <- V(graph_net)$name
V(graph_net)$degree <- degree(graph_net)

plot(graph_net, vertex.size = 20, vertex.label.cex = 0.8)

write_graph(graph_net, "C:/Users/User/OneDrive/Desktop/qian/R/network/18os+cyc-ecoli/os+cyc-L8-ecoli.graphml", format = "graphml")
write_graph(graph_net, "C:/Users/User/OneDrive/Desktop/qian/R/network/18os+cyc-ecoli/os+cyc-L8-ecoli.gml", format = "gml")