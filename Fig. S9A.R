# Clean and Publication-Ready Version of Code for Protist RA Boxplots (r < 0)

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
rm(list = ls())

# Set directory and load relative abundance data
setwd('C:/Users/User/OneDrive/Desktop/qian/R/cyc+-relative_abundance(lefse)/18dna-os')
RA_data <- read.table("level-8.txt", sep = "\t", header = TRUE, check.names = FALSE)

# Transform to long format and normalize to %
RA_data <- RA_data %>%
  pivot_longer(cols = -group, names_to = "Sample", values_to = "RelativeAbundance") %>%
  mutate(RelativeAbundance = RelativeAbundance * 100,
         group = gsub("st__un_d[0-9]+__", "un_", group),
         group = gsub("d8__", "", group)) %>%
  mutate(Treatment = case_when(
    grepl("^OSc-", Sample) ~ "+CHX",
    grepl("^OSn-", Sample) ~ "-CHX"
  ),
  Time = case_when(
    grepl("-14$", Sample) ~ "Day14",
    grepl("-28$", Sample) ~ "Day28",
    grepl("-0$", Sample) ~ "Day0"
  ))

# Load trophic assignment data
setwd('C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-os')
trophicdata <- read.table("df_filtered-trophic.txt", sep = "\t", header = TRUE, check.names = FALSE) %>%
  rename(name = Rank8) %>%
  mutate(name = gsub("st__un_d[0-9]+__", "un_", name),
         name = gsub("d8__", "", name)) %>%
  distinct(name, .keep_all = TRUE)

# Load E. coli network correlations (r < 0 only)
setwd('C:/Users/User/OneDrive/Desktop/qian/R/network/18os-cyc-ecoli')
ecoli_r <- read.csv("R.ecoli_level8_sig.csv", sep = ",") %>%
  select(name = X, rvalue = E.coli) %>%
  filter(name != "E.coli" & rvalue < 0)

# Merge to retain only consumers
merged_data <- ecoli_r %>%
  left_join(trophicdata %>% select(name, Trophic), by = "name") %>%
  filter(Trophic == "Consumer") %>%
  distinct()

# Define output directory
output_dir <- "C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/key protists RA change_scatterplot_18os/network-cyc_r_less0/Consumer"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Filter and duplicate Day0 -CHX samples for visual comparison
plot_list <- list()
num_cols <- 4

for (i in seq_len(nrow(merged_data))) {
  taxon <- merged_data$name[i]
  data <- RA_data %>% filter(group == taxon)
  
  if (max(data$RelativeAbundance, na.rm = TRUE) >= 1) {
    duplicated <- data %>% filter(Time == "Day0" & Treatment == "-CHX") %>% mutate(Treatment = "+CHX")
    updated_data <- bind_rows(data, duplicated)
    updated_data$Time <- factor(updated_data$Time, levels = c("Day0", "Day14", "Day28"))
    
    p <- ggplot(updated_data, aes(x = Treatment, y = RelativeAbundance, color = Time)) +
      geom_boxplot(position = position_dodge(0.6), width = 0.5, size = 0.8, outlier.shape = NA) +
      geom_jitter(position = position_jitterdodge(0.15, 0.6), size = 3, alpha = 0.7) +
      stat_summary(fun = median, geom = "crossbar", width = 0.4, size = 0.4,
                   position = position_dodge(0.5), aes(color = Time), show.legend = FALSE) +
      labs(title = taxon, x = NULL, y = "Relative abundance (%)") +
      scale_y_continuous(labels = number_format(accuracy = 0.1), limits = c(0, max(updated_data$RelativeAbundance) * 1.1)) +
      scale_color_manual(values = c("Day0" = "#3A3737", "Day14" = "#1f77b4", "Day28" = "#ffcc00")) +
      theme_classic(base_family = "Arial") +
      theme(
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom"
      )
    
    # Save individual plots
    fname <- file.path(output_dir, paste0("taxa_", i, "_", gsub("[^a-zA-Z0-9]", "_", taxon), ".jpg"))
    ggsave(fname, plot = p, width = 9, height = 8, units = "cm", dpi = 600)
    
    # Adjust y-label only for leftmost plot
    if ((i - 1) %% num_cols == 0) {
      p <- p + labs(y = "Relative abundance (%)")
    } else {
      p <- p + labs(y = "")
    }
    
    plot_list[[length(plot_list) + 1]] <- p
  }
}

# Combine plots
combined_plot <- wrap_plots(plotlist = plot_list, ncol = num_cols) +
  plot_annotation(title = "Organic soil",
                  theme = theme(plot.title = element_text(size = 30, hjust = 0.5)))

# Save combined figure
ggsave(file.path(output_dir, "combined_taxa_plots_r_less0_box.jpg"),
       plot = combined_plot, width = 28, height = 11, units = "cm", dpi = 600)

print(combined_plot)
