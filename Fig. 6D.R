library(dplyr)
library(tibble)
library(MASS)
library(reshape2)
library(Hmisc)
library(igraph)
library(ggplot2)
library(ggpubr)
library(ggsignif)

############################# Data Processing: CS + CHX (Day 0 to Day 28) ##########################
# Step 1: Load 18S data (OTU table)
rm(list = ls())
setwd('C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-cs')
metadata18 <- read.table("group(+cyc).txt", sep = "\t", header = TRUE)
otu18 <- read.table("level-8.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("sample")

mergedtable18 <- merge(otu18, metadata18, by = "sample") %>% 
  dplyr::select(-sample) %>% 
  relocate(Replicate, .before = everything()) %>% 
  t() %>% as.data.frame()

column_names <- as.character(mergedtable18[1, ])
mergedtable18 <- mergedtable18[-1, ]
colnames(mergedtable18) <- column_names
mergedtable18 <- mergedtable18 %>% mutate_all(as.numeric)

# Add pseudocount and apply CLR transformation
mergedtable18_clr <- microbiome::transform(mergedtable18 + 1, "clr") %>% as.data.frame()

# Step 2: Load 16S data and retain only Escherichia-Shigella (E. coli)
setwd('C:/Users/User/OneDrive/Desktop/qian/R/CYC-with E.coil+unknow/16dna-cs')
metadata16 <- read.table("group(+cyc).txt", sep = "\t", header = TRUE)
otu16 <- read.table("genus.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("sample")

mergedtable16 <- merge(otu16, metadata16, by = "sample") %>% 
  dplyr::select(-sample) %>% 
  relocate(Replicate, .before = everything()) %>% 
  t() %>% as.data.frame()

column_names <- as.character(mergedtable16[1, ])
mergedtable16 <- mergedtable16[-1, ]
colnames(mergedtable16) <- column_names
mergedtable16 <- mergedtable16 %>% mutate_all(as.numeric)

mergedtable16_clr <- microbiome::transform(mergedtable16 + 1, "clr") %>% as.data.frame() %>% 
  t() %>% as.data.frame() %>% dplyr::select(`g__Escherichia-Shigella`) %>% 
  t() %>% as.data.frame()

# Step 3: Standardize taxonomic labels and merge tables
mergedtable18_clr <- rownames_to_column(mergedtable18_clr, var = "Row_Name")
mergedtable16_clr <- rownames_to_column(mergedtable16_clr, var = "Row_Name")

mergedtable18_clr$Row_Name <- gsub("d8__|st__un_d[2-7]__", "un_", mergedtable18_clr$Row_Name)
mergedtable16_clr$Row_Name <- gsub("g__Escherichia-Shigella", "E.coli", mergedtable16_clr$Row_Name)

mergedtable18_clr <- column_to_rownames(mergedtable18_clr, var = "Row_Name")
mergedtable16_clr <- column_to_rownames(mergedtable16_clr, var = "Row_Name")

mergedtable18_clr <- mergedtable18_clr %>% t() %>% as.data.frame() %>% rownames_to_column("sample")
mergedtable16_clr <- mergedtable16_clr %>% t() %>% as.data.frame() %>% rownames_to_column("sample")

mergedtableall <- merge(mergedtable18_clr, mergedtable16_clr, by = "sample") %>% 
  t() %>% as.data.frame()

column_names <- as.character(mergedtableall[1, ])
mergedtableall <- mergedtableall[-1, ]
colnames(mergedtableall) <- column_names
mergedtableall <- mergedtableall %>% mutate_all(as.numeric)

# Step 4: Filter taxa with >10% prevalence
dt <- mergedtableall
dt[dt > 1] <- 1
prevalent <- which(rowSums(dt) / ncol(dt) > 0.1)
mergedtableall <- mergedtableall[prevalent, ] %>% as.matrix()

# Step 5: Compute Spearman correlation
sp.cor <- rcorr(t(mergedtableall), type = "spearman")
r.cor <- sp.cor$r
p.cor <- sp.cor$P
write.csv(r.cor, ".../rvaluecor_level8.csv")
write.csv(p.cor, ".../Pvaluecor_level8.csv")

# Step 6: Multiple testing correction and threshold filtering
p.adj <- matrix(p.adjust(p.cor, method = "BH"), nrow = nrow(p.cor), dimnames = dimnames(p.cor))
r.cor[p.cor >= 0.05 | abs(r.cor) < 0.6] <- 0
p.cor[p.cor >= 0.05 | abs(r.cor) < 0.6] <- 0
p.cor[is.na(p.cor)] <- 1

write.csv(r.cor, ".../Rvalue_filter_level8.csv")
write.csv(p.cor, ".../Pvalue_filter_level8.csv")

# Step 7: Extract significant correlations involving E. coli
r.ecoli <- r.cor %>% as.data.frame() %>% filter(abs(E.coli) > 0) %>% t() %>% 
  as.data.frame() %>% filter(abs(E.coli) > 0)
diag(r.ecoli) <- 0
write.csv(r.ecoli, ".../R.ecoli_level8_sig.csv")

# Step 8: Create undirected network
r.ecoli <- read.csv(".../R.ecoli_level8_sig.csv", row.names = 1)
colnames(r.ecoli)[colnames(r.ecoli) == "un_Rhogostoma.lineage"] <- "un_Rhogostoma-lineage"
r.ecoli <- as.matrix(r.ecoli)

A <- graph_from_adjacency_matrix(r.ecoli, weight = TRUE, mode = "undirected")
E(A)$width <- abs(E(A)$weight) * 5

# Step 9: Annotate network with abundance and taxonomic data
data <- rowSums(mergedtableall)
data1 <- as.data.frame(data)
vertex_names <- V(A)$name
df_igraph_size <- data1[vertex_names, , drop = FALSE]
V(A)$Abundance <- log10(df_igraph_size$data)

# Load taxonomy info
data2 <- read.table(".../level2-8.txt", sep = "\t", header = TRUE, check.names = FALSE, row.names = 1) %>% 
  rownames_to_column("Row_Name") %>% semi_join(data1 %>% rownames_to_column("Row_Name"), by = "Row_Name") %>% 
  column_to_rownames("Row_Name")

V(A)$G <- as.character(data2[V(A)$name, ])

# Add graph attributes
V(A)$label <- V(A)$name
V(A)$degree <- degree(A)

# Export network
write_graph(A, ".../cs+cyc-L8-ecoli.graphml", format = "graphml")
write_graph(A, ".../cs+cyc-L8-ecoli.gml", format = "gml")
