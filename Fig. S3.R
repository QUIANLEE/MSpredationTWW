library(phyloseq)
library(vegan)
rm(list = ls())
############Conventional soil#######
setwd(dir = 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-cs')

otu <- "featuretable.qza"
rep <- "rep-seqs.qza"
tree <- "rooted-tree.qza"
tax <- "taxonomy.qza"
sample <- "sample-metadata-dnacs.tsv"
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,
                          mapfilename=sample,treeio=tree)

ps <-ps_dada2 

# Extract OTU table as matrix
otu_mat <- as(otu_table(ps), "matrix")

# If taxa are rows instead of columns, transpose the matrix
if (taxa_are_rows(ps)) {
  otu_mat <- t(otu_mat)
}

# Plot rarefaction curves
tiff("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/Conventional_rarefaction_curves.tiff", 
     width = 13, height = 12, units = "cm", res = 600)
rarecurve(otu_mat, 
          step = 100, 
          sample = min(rowSums(otu_mat)), 
          col = "blue", 
          cex = 0.6, 
          label = FALSE)
title(main = "Conventional soil")
dev.off()  # ✔️ Close the TIFF device

# ✔️ Save the same plot as JPEG
jpeg("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/Conventional_rarefaction_curves.jpg", 
     width = 13, height = 12, units = "cm", res = 600)
rarecurve(otu_mat, 
          step = 100, 
          sample = min(rowSums(otu_mat)), 
          col = "blue", 
          cex = 0.6, 
          label = FALSE)
title(main = "Conventional soil")
dev.off()  # ✔️ Close the JPEG device

############Organic soil#######
setwd(dir = 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os')
otu <- "featuretable.qza"
rep <- "rep-seqs.qza"
tree <- "rooted-tree.qza"
tax <- "taxonomy.qza"
sample <- "sample-metadata-dnaos.tsv"
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,
                          mapfilename=sample,treeio=tree)
ps <- ps_dada2

# Extract OTU table as matrix
otu_mat <- as(otu_table(ps), "matrix")

# If taxa are rows instead of columns, transpose the matrix
if (taxa_are_rows(ps)) {
  otu_mat <- t(otu_mat)
}

# Plot rarefaction curves
tiff("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/Organic_rarefaction_curves.tiff", 
     width = 13, height = 12, units = "cm", res = 600)
rarecurve(otu_mat, 
          step = 100, 
          sample = min(rowSums(otu_mat)), 
          col = "blue", 
          cex = 0.6, 
          label = FALSE)
title(main = "Organic soil")
dev.off()  # ✔️ Close the TIFF device

# ✔️ Save the same plot as JPEG
jpeg("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/Organic_rarefaction_curves.jpg", 
     width = 13, height = 12, units = "cm", res = 600)
rarecurve(otu_mat, 
          step = 100, 
          sample = min(rowSums(otu_mat)), 
          col = "blue", 
          cex = 0.6, 
          label = FALSE)
title(main = "Organic soil")
dev.off()  # ✔️ Close the JPEG device

############Code to Merge Two Boxplots into One Figure#####################
merged_plot <- ggarrange(b+
                           theme(
                             #axis.title.x = element_blank(),
                             axis.title.y = element_text(size = 15, family = "Arial", color = "black",margin = margin(r = 15)),
                             plot.margin = margin(5, -5, 5, 5)  # Remove all margins
                           ),
                         a+
                           theme(
                             axis.title.y = element_blank(),
                             #axis.text.y = element_blank(),
                             #axis.ticks.y = element_blank(),
                             #axis.title.x = element_blank(),
                             plot.margin = margin(5, 5, 5, -5)  # Remove all margins
                           ),
                         labels = c("A", "B"),  # Add labels for identification
                         font.label = list(size = 20),
                         ncol = 2, nrow = 1,  # Arrange in one row, two columns
                         align = "hv")  # Align horizontally and vertically


print(merged_plot)

# Save the merged plot as an image
ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/merged_rarefaction_curves.tiff", plot = merged_plot, device = "tiff",
       dpi = 600, width = 25, height = 12, units = "cm")

ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/merged_rarefaction_curves.jpg", plot = merged_plot, device = "jpeg",
       dpi = 600, width = 25, height = 12, units = "cm")