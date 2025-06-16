library(MicrobiotaProcess)
library(phyloseq)
library(tidyverse)
library(vegan)
library(multcompView)
library(ggpubr)
library(RColorBrewer)
library(pacman)
library(magrittr)
library(reshape2)
library(ggsci)
library(FSA)  # For dunnTest
library(scales) 
library(ggpubr)
rm(list = ls())
############################Shannon index####################
##########Conventional soil######
setwd(dir = 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-cs')

otu <- "featuretable.qza"
rep <- "rep-seqs.qza"
tree <- "rooted-tree.qza"
tax <- "taxonomy.qza"
sample <- "sample-metadata-dnacs.tsv"
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,
                          mapfilename=sample,treeio=tree)

ps <-ps_dada2 

# Rarefy the phyloseq object 
set.seed(123)  # for reproducibility
ps_rarefied <- rarefy_even_depth(ps_dada2, 
                                 sample.size = 1274, 
                                 rngseed = 123, 
                                 replace = FALSE, 
                                 verbose = TRUE)

alphaobj <- get_alphaindex(ps_rarefied)
head(as.data.frame(alphaobj))

# Step 1: Extract Shannon index and sample data from alphaobj
alpha_shannon <- data.frame(
  Shannon = alphaobj@alpha$Shannon,
  group = alphaobj@sampleda$group,
  row.names = rownames(alphaobj@sampleda)  # Set row names to match the original alphaobj
)

shannon <- alpha_shannon %>%
  mutate(
    Treatment = case_when(
      grepl("c-14", group) ~ "Day14(+CHX)",   # If Sample starts with "CSc-"
      grepl("n-14", group) ~ "Day14",    # If Sample starts with "CSn-"
      grepl("c-28", group) ~ "Day28(+CHX)",   # If Sample starts with "CSc-"
      grepl("n-28", group) ~ "Day28",    # If Sample starts with "CSn-"
      grepl("n-0", group) ~ "Day0"
    )
  )%>%
  dplyr::select(-group)

# Step 2: Check the structure of the new data frame
head(alpha_shannon)
alpha_shannon$group <- gsub("CSn-14", "14", alpha_shannon$group)
alpha_shannon$group <- gsub("CSn-28", "28", alpha_shannon$group)
alpha_shannon$group <- gsub("CSn-0", "0", alpha_shannon$group)
alpha_shannon$group <- gsub("CSc-14", "14\n(CHX)", alpha_shannon$group)
alpha_shannon$group <- gsub("CSc-28", "28\n(CHX)", alpha_shannon$group)
alpha_shannon$group <- factor(alpha_shannon$group, levels = c("0", "14","28","14\n(CHX)","28\n(CHX)"))

#####################Kruskal–Wallis and Dunn's Test
# Step 1: Perform Kruskal–Wallis Test
kruskal_result <- kruskal.test(Shannon ~ group, data = alpha_shannon)
print(kruskal_result)

# Step 2: Perform Dunn's Test for Post Hoc Comparisons
dunn_result <- dunnTest(Shannon ~ group, data = alpha_shannon, method = "bonferroni")

# Extract Dunn's results
dunn_df <- as.data.frame(dunn_result$res)

# Rename columns for clarity
colnames(dunn_df) <- c("Comparison", "Z", "p.unadj", "p.adj")

# Step 3: Assign Significance Letters Correctly
# Build a pairwise matrix for multcompLetters
all_groups <- levels(alpha_shannon$group)
comparison_split <- strsplit(dunn_df$Comparison, " - ")
group1 <- sapply(comparison_split, `[`, 1)
group2 <- sapply(comparison_split, `[`, 2)

p_mat <- matrix(1, nrow = length(all_groups), ncol = length(all_groups),
                dimnames = list(all_groups, all_groups))

for (i in seq_along(group1)) {
  g1 <- group1[i]
  g2 <- group2[i]
  p_val <- dunn_df$p.adj[i]
  p_mat[g1, g2] <- p_val
  p_mat[g2, g1] <- p_val
}

letters_obj <- multcompLetters(p_mat, threshold = 0.05)

# Create a data frame with group labels
group_labelsTH <- data.frame(
  group = names(letters_obj$Letters),
  label = letters_obj$Letters
)

# Merge with all groups and handle missing labels
group_labelsTH <- merge(data.frame(group = all_groups), group_labelsTH, by = "group", all.x = TRUE)

if (any(is.na(group_labelsTH$label))) {
  most_common_label <- names(sort(table(group_labelsTH$label), decreasing = TRUE))[1]
  group_labelsTH$label[is.na(group_labelsTH$label)] <- most_common_label
}

# Step 4: Merge with Maximum Shannon Index for Plotting
group_max_values <- alpha_shannon %>%
  group_by(group) %>%
  summarise(max_value = max(Shannon), .groups = "drop")

group_labelsTH <- group_labelsTH %>%
  left_join(group_max_values, by = "group") %>%
  mutate(y_position = max_value + 0.1)

# Print final result
print(group_labelsTH)

##########################plotting
# Generate the box plot for Shannon index using ggplot2

p_alpha <- ggplot(alpha_shannon, aes(x = group, y = Shannon, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FEB4A9","#9AC9DB", "#40588F", "#FFF795","#F9C00F")) +
  ylim(3, 5) +  # Set the y-axis range
  theme_bw() +  # Set a base theme and dont change the postion to the end, otherwise it overrides all previous settings
  theme(
    legend.position = "none",  # Remove the legend
    axis.text.x = element_text(size = 15, color = "black", family = "Arial"),
    axis.text.y = element_text(size = 15, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 15, color = "black", family = "Arial"),
    axis.title.x = element_text(size = 15, color = "black", family = "Arial"),
    plot.title = element_text(size = 20, hjust = 0.5, family = "Arial"),
    panel.grid = element_blank()  # Remove all grid lines
  ) +
  labs(title = "Conventional soil", x = "Time(Days)", y = "Shannon Index")
print(p_alpha)


# Add significance letters (Dunn's test) to the plot
alpha_plot_cs <- p_alpha +
  geom_text(data = group_labelsTH, aes(x = group, y = y_position, label = label),
            color = "black", fontface = "bold", size = 5)

# Always display "p < 0.05" annotation from Kruskal-Wallis
p_alpha_TH_cs <- alpha_plot_cs +
  annotate("text", x = 1, y = 4.95,
           label = "padj< 0.05",
           size = 5)

# Print the plot
print(p_alpha_TH_cs)


############Organic soil################
setwd(dir = 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os')
otu <- "featuretable.qza"
rep <- "rep-seqs.qza"
tree <- "rooted-tree.qza"
tax <- "taxonomy.qza"
sample <- "sample-metadata-dnaos.tsv"
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax,refseqqza=rep,
                          mapfilename=sample,treeio=tree)
ps <- ps_dada2

# Rarefy the phyloseq 
set.seed(123)  # for reproducibility
ps_rarefied <- rarefy_even_depth(ps_dada2, 
                                 sample.size = 1435, 
                                 rngseed = 123, 
                                 replace = FALSE, 
                                 verbose = TRUE)


alphaobj <- get_alphaindex(ps_rarefied)
head(as.data.frame(alphaobj))

# Step 1: Extract Shannon index and sample data from alphaobj
alpha_shannon <- data.frame(
  Shannon = alphaobj@alpha$Shannon,
  group = alphaobj@sampleda$group,
  row.names = rownames(alphaobj@sampleda)  # Set row names to match the original alphaobj
)

shannon <- alpha_shannon %>%
  mutate(
    Treatment = case_when(
      grepl("c-14", group) ~ "Day14(+CHX)",   # If Sample starts with "CSc-"
      grepl("n-14", group) ~ "Day14",    # If Sample starts with "CSn-"
      grepl("c-28", group) ~ "Day28(+CHX)",   # If Sample starts with "CSc-"
      grepl("n-28", group) ~ "Day28",    # If Sample starts with "CSn-"
      grepl("n-0", group) ~ "Day0"
    )
  )%>%
  dplyr::select(-group)

write.table(shannon, 'C:/Users/User/OneDrive/Desktop/qian/1. PhD project/Project1-Valli data/3.Paper 1 results/shannon_index_OS.txt', row.names = FALSE, sep = '\t', quote = FALSE)

# Step 2: Check the structure of the new data frame
head(alpha_shannon)

alpha_shannon$group <- gsub("OSn-14", "14", alpha_shannon$group)
alpha_shannon$group <- gsub("OSn-28", "28", alpha_shannon$group)
alpha_shannon$group <- gsub("OSn-0", "0", alpha_shannon$group)
alpha_shannon$group <- gsub("OSc-14", "14\n(CHX)", alpha_shannon$group)
alpha_shannon$group <- gsub("OSc-28", "28\n(CHX)", alpha_shannon$group)
alpha_shannon$group <- factor(alpha_shannon$group, levels = c("0", "14", "28","14\n(CHX)","28\n(CHX)"))


#####################Kruskal–Wallis and Dunn's Test
# Step 1: Perform Kruskal–Wallis Test
kruskal_result <- kruskal.test(Shannon ~ group, data = alpha_shannon)
print(kruskal_result)

# Step 2: Perform Dunn's Test for Post Hoc Comparisons
dunn_result <- dunnTest(Shannon ~ group, data = alpha_shannon, method = "bonferroni")

# Extract Dunn's results
dunn_df <- as.data.frame(dunn_result$res)

# Rename columns for clarity
colnames(dunn_df) <- c("Comparison", "Z", "p.unadj", "p.adj")

# Export as a .txt file
write.table(dunn_df, file = "C:/Users/User/OneDrive/Desktop/qian/1. PhD project/Project1-Valli data/3.Paper 1 results/shannon_DunnTest_OS.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Step 3: Assign Significance Letters Correctly
all_groups <- levels(alpha_shannon$group)
comparison_split <- strsplit(dunn_df$Comparison, " - ")
group1 <- sapply(comparison_split, `[`, 1)
group2 <- sapply(comparison_split, `[`, 2)

p_mat <- matrix(1, nrow = length(all_groups), ncol = length(all_groups),
                dimnames = list(all_groups, all_groups))

for (i in seq_along(group1)) {
  g1 <- group1[i]
  g2 <- group2[i]
  p_val <- dunn_df$p.adj[i]
  p_mat[g1, g2] <- p_val
  p_mat[g2, g1] <- p_val
}

letters_obj <- multcompLetters(p_mat, threshold = 0.05)

# Create a data frame with group labels
group_labelsTH <- data.frame(
  group = names(letters_obj$Letters),
  label = letters_obj$Letters
)

# Merge with all groups and handle missing labels
group_labelsTH <- merge(data.frame(group = all_groups), group_labelsTH, by = "group", all.x = TRUE)

if (any(is.na(group_labelsTH$label))) {
  most_common_label <- names(sort(table(group_labelsTH$label), decreasing = TRUE))[1]
  group_labelsTH$label[is.na(group_labelsTH$label)] <- most_common_label
}

# Step 4: Merge with Maximum Shannon Index for Plotting
group_max_values <- alpha_shannon %>%
  group_by(group) %>%
  summarise(max_value = max(Shannon), .groups = "drop")

group_labelsTH <- group_labelsTH %>%
  left_join(group_max_values, by = "group") %>%
  mutate(y_position = max_value + 0.1)

# Print final group labels
print(group_labelsTH)

##########################plotting
# Generate the box plot for Shannon index using ggplot2

p_alpha_os <- ggplot(alpha_shannon, aes(x = group, y = Shannon, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FEB4A9","#9AC9DB", "#40588F", "#FFF795","#F9C00F")) +
  ylim(3, 5) +  # Set the y-axis range
  theme_bw() +  # Set a base theme and dont change the postion to the end, otherwise it overrides all previous settings
  theme(
    legend.position = "none",  # Remove the legend
    axis.text.x = element_text(size = 15, color = "black", family = "Arial"),
    axis.text.y = element_text(size = 15, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 15, color = "black", family = "Arial"),
    axis.title.x = element_text(size = 15, color = "black", family = "Arial"),
    plot.title = element_text(size = 20, hjust = 0.5, family = "Arial"),
    panel.grid = element_blank()  # Remove all grid lines
  ) +
  labs(title = "Organic soil", x = "Time(Days)", y = "Shannon Index")
print(p_alpha_os)

# Add significance letters (Dunn's test) to the plot
alpha_plot_os <- p_alpha_os +
  geom_text(data = group_labelsTH, aes(x = group, y = y_position, label = label),
            color = "black", fontface = "bold", size = 5)

# Add Kruskal-Wallis p-value annotation (fixed as "p < 0.05")
p_alpha_TH_os <- alpha_plot_os +
  annotate("text", x = 1, y = 4.95,
           label = "padj< 0.05",
           size = 5)

# Print the plot
print(p_alpha_TH_os)



#####################PCoA#####################
############Conventional soil######
# Set working directory
setwd(dir = 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-cs')

# Import Qiime2 data
otu <- "featuretable.qza"
rep <- "rep-seqs.qza"
tree <- "rooted-tree.qza"
tax <- "taxonomy.qza"
sample <- "sample-metadata-dnacs.tsv"
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax, refseqqza=rep,
                          mapfilename=sample, treeio=tree)
ps <- ps_dada2

# Rarefy the phyloseq object 
set.seed(123)  # for reproducibility
ps_rarefied <- rarefy_even_depth(ps_dada2, 
                                 sample.size = 1274, 
                                 rngseed = 123, 
                                 replace = FALSE, 
                                 verbose = TRUE)

# Calculate Bray-Curtis distance matrix
dist_matrix <- phyloseq::distance(ps_rarefied, method = "bray")

# Extract sample data
sample_data_df <- data.frame(sample_data(ps_rarefied))

# Set grouping variable with correct factor levels
sample_data_df$group <- gsub("CSn-14", "14", sample_data_df$group)
sample_data_df$group <- gsub("CSn-28", "28", sample_data_df$group)
sample_data_df$group <- gsub("CSn-0", "0", sample_data_df$group)
sample_data_df$group <- gsub("CSc-14", "14(CHX)", sample_data_df$group)
sample_data_df$group <- gsub("CSc-28", "28(CHX)", sample_data_df$group)
sample_data_df$group <- factor(sample_data_df$group, 
                               levels = c("0", "14", "28","14(CHX)","28(CHX)"))

# Update phyloseq object with updated sample data
sample_data(ps_rarefied) <- sample_data(sample_data_df)

# Define the grouping variable
group <- sample_data_df$group

# Perform PERMANOVA test
library(vegan)
permanova_result <- adonis2(dist_matrix ~ group, data = sample_data_df, permutations = 999)
print(permanova_result)

# Extract R² and p-value
R2_value <- round(permanova_result$R2[1], 3)
p_value <- permanova_result$`Pr(>F)`[1]

# Print the values
cat("permanova R2 =", R2_value, "\nP-value =", p_value, "\n")

# Perform PCoA analysis using ordinate()
pcoa_res <- ordinate(ps_rarefied, method = "PCoA", distance = "bray")

# Extract PCoA scores and convert them to a data frame
pcoa_data <- data.frame(pcoa_res$vectors[, 1:2])  # Extract the first two axes
pcoa_data$group <- sample_data_df$group  # Add group information to the data frame

# Extract explained variance for each axis
eig_values <- pcoa_res$values$Relative_eig

# Plot PCoA with ggplot2
pcoaplot_cs <- ggplot(pcoa_data, aes(x = Axis.1, y = Axis.2, color = group)) +
  geom_point(size = 9, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, linetype = 2) +  # Add ellipses
  scale_color_manual(values = c("#FEB4A9","#9AC9DB", "#40588F", "#FFF795","#F9C00F")) +
  # Add dashed perpendicular (y-axis) and parallel (x-axis) lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Dashed horizontal line at y = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Dashed vertical line at x = 0
  labs(#title = "Conventional soil",
       x = paste0("PCoA1 (", round(eig_values[1] * 100, 2), "%)"),
       y = paste0("PCoA2 (", round(eig_values[2] * 100, 2), "%)"),
       color = "Time (Days)") +
  # Format x and y axis scales to have one decimal place
  scale_x_continuous(
    labels = number_format(accuracy = 0.1),
    limits = c(-0.5, 0.5)  # ✅ Set x-axis range (adjust as needed)
  ) +
  scale_y_continuous(
    labels = number_format(accuracy = 0.1),
    limits = c(-0.4, 0.6)  # ✅ Set y-axis range (adjust as needed)
  ) +
  # Add ANOSIM annotation
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.8,
           label = bquote("PERMANOVA " ~ R^2 ~ "=" ~ .(R2_value) ~ "\n" ~ italic(p) ~ "=" ~ "< 0.001"),
           size = 5, color = "black", family = "Arial") +
  # Theme adjustments to remove background and grid lines
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add black border
    panel.grid = element_blank(),  # Remove all grid lines
    axis.title.x = element_text(size = 15, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial"),
    axis.text.x = element_text(size = 15, color = "black", family = "Arial"),
    axis.text.y = element_text(size = 15, color = "black", family = "Arial"),
    plot.title = element_text(size = 20, hjust = 0.5, family = "Arial"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),  # Increase the text size of the legend
    legend.key.size = unit(0.4, "lines"), # Corrected syntax, removed misplaced 
    legend.key.height = unit(0.3, "lines"),  # Adjust legend key height
    legend.key.width = unit(0.3, "lines"),  # Adjust legend key width
    legend.position = c(0.85, 0.8),  # Move legend inside the plot box (adjust coordinates as needed)
  ) +
  guides(color = guide_legend(override.aes = list(linetype = NA, shape = 16, size = 5)))  # Remove lines from legend  # Guides should be outside `theme()`

print(pcoaplot_cs)

############Organic soil##########
# Set working directory
setwd(dir = 'C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os')

# Import Qiime2 data
otu <- "featuretable.qza"
rep <- "rep-seqs.qza"
tree <- "rooted-tree.qza"
tax <- "taxonomy.qza"
sample <- "sample-metadata-dnaos.tsv"
ps_dada2 <- import_qiime2(otuqza=otu, taxaqza=tax, refseqqza=rep,
                          mapfilename=sample, treeio=tree)
ps <- ps_dada2

# Rarefy the phyloseq object 
set.seed(123)  # for reproducibility
ps_rarefied <- rarefy_even_depth(ps_dada2, 
                                 sample.size = 1435, 
                                 rngseed = 123, 
                                 replace = FALSE, 
                                 verbose = TRUE)

# Calculate Bray-Curtis distance matrix
dist_matrix <- phyloseq::distance(ps_rarefied, method = "bray")

# Extract sample data
sample_data_df <- data.frame(sample_data(ps_rarefied))

# Set grouping variable with correct factor levels
sample_data_df$group <- gsub("OSn-14", "14", sample_data_df$group)
sample_data_df$group <- gsub("OSn-28", "28", sample_data_df$group)
sample_data_df$group <- gsub("OSn-0", "0", sample_data_df$group)
sample_data_df$group <- gsub("OSc-14", "14(CHX)", sample_data_df$group)
sample_data_df$group <- gsub("OSc-28", "28(CHX)", sample_data_df$group)
sample_data_df$group <- factor(sample_data_df$group, 
                               levels = c("0", "14", "28","14(CHX)","28(CHX)"))

# Update phyloseq object with updated sample data
sample_data(ps_rarefied) <- sample_data(sample_data_df)

# Perform PERMANOVA test
permanova_result <- adonis2(dist_matrix ~ group, data = sample_data_df, permutations = 999)
print(permanova_result)

# Extract R² and p-value
R2_value <- round(permanova_result$R2[1], 3)
p_value <- permanova_result$`Pr(>F)`[1]

# Print the values
cat("R =", R2_value, "\nP-value =", p_value, "\n")

# Perform PCoA analysis using ordinate()
pcoa_res <- ordinate(ps_rarefied, method = "PCoA", distance = "bray")
# Extract PCoA scores and convert them to a data frame
pcoa_data <- data.frame(pcoa_res$vectors[, 1:2])  # Extract the first two axes
pcoa_data$group <- sample_data_df$group  # Add group information to the data frame

# Extract explained variance for each axis
eig_values <- pcoa_res$values$Relative_eig

# Plot PCoA with ggplot2
pcoaplot_os <- ggplot(pcoa_data, aes(x = Axis.1, y = Axis.2, color = group)) +
  geom_point(size = 9, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, linetype = 2) +  # Add ellipses
  scale_color_manual(values = c("#FEB4A9","#9AC9DB", "#40588F", "#FFF795","#F9C00F")) +
  # Add dashed perpendicular (y-axis) and parallel (x-axis) lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Dashed horizontal line at y = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Dashed vertical line at x = 0
  labs(#title = "Organic soil",
       x = paste0("PCoA1 (", round(eig_values[1] * 100, 2), "%)"),
       y = paste0("PCoA2 (", round(eig_values[2] * 100, 2), "%)"))+
  #color = "Time (Days)") +
  # Format x and y axis scales to have one decimal place
  scale_x_continuous(
    labels = number_format(accuracy = 0.1),
    limits = c(-0.5, 0.5)  # ✅ Set x-axis range (adjust as needed)
  ) +
  scale_y_continuous(
    labels = number_format(accuracy = 0.1),
    limits = c(-0.4, 0.6)  # ✅ Set y-axis range (adjust as needed)
  ) +
  # Add PERMANOVA annotation
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.8,
           label = bquote("PERMANOVA " ~ R^2 ~ "=" ~ .(R2_value) ~ "\n" ~ italic(p) ~ "=" ~ "< 0.001"),
           size = 5, color = "black", family = "Arial") +
  # Theme adjustments to remove background and grid lines
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Add black border
    panel.grid = element_blank(),  # Remove all grid lines
    axis.title.x = element_text(size = 15, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial"),
    axis.text.x = element_text(size = 15, color = "black", family = "Arial"),
    axis.text.y = element_text(size = 15, color = "black", family = "Arial"),
    plot.title = element_text(size = 20, hjust = 0.5, family = "Arial"),
    legend.position = "none")

print(pcoaplot_os)



############Code to Merge several Boxplots into One Figure#####################

combined_plot <- ggarrange(
  p_alpha_TH_os+
    theme(
      axis.title.x = element_blank()
    ), 
  p_alpha_TH_cs+
    theme(
      axis.title.x = element_blank()
    ), 
  pcoaplot_os, 
  pcoaplot_cs,
  labels = c("A", "B", "C", "D"),  # Add subplot labels
  ncol = 2, nrow = 2,              # Arrange in a 2x2 grid
  common.legend = FALSE)           # Share a common legend (optional)


#legend = "right"  # Legend position: "top", "bottom", "left", or "right"
# Print the combined plot
print(combined_plot)

# Save the combined plot
ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/merged_PCoA_Shannon-new.jpg", plot = combined_plot, device = "jpeg",
       dpi = 600, width = 25, height = 25, units = "cm")

ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/merged_PCoA_Shannon-new.tiff", plot = combined_plot, device = "tiff",
       dpi = 600, width = 25, height = 25, units = "cm")
