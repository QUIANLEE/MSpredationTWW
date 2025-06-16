# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
rm(list = ls())

# Load and preprocess relative abundance data
setwd("C:/Users/User/OneDrive/Desktop/qian/R/cyc+-relative_abundance(lefse)/18dna-cs")
RA_data <- read.table("level-8.txt", sep = "\t", header = TRUE, check.names = FALSE)
RA_data <- RA_data %>%
  pivot_longer(cols = -group, names_to = "Sample", values_to = "RelativeAbundance") %>%
  mutate(RelativeAbundance = RelativeAbundance * 100,
         group = gsub("st__un_d[0-9]+__", "un_", group),
         group = gsub("d[0-9]+__", "", group),
         Treatment = case_when(grepl("^CSc-", Sample) ~ "+CHX", grepl("^CSn-", Sample) ~ "-CHX"),
         Time = case_when(grepl("-14$", Sample) ~ "Day14", grepl("-28$", Sample) ~ "Day28", grepl("-0$", Sample) ~ "Day0"))

# Load and preprocess trophic classification
setwd("C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-cs")
trophicdata <- read.table("df_filtered-trophic.txt", sep = "\t", header = TRUE, check.names = FALSE)
colnames(trophicdata)[colnames(trophicdata) == "Rank8"] <- "name"
trophicdata$name <- gsub("st__un_d[0-9]+__", "un_", trophicdata$name)
trophicdata$name <- gsub("d[0-9]+__", "", trophicdata$name)
trophicdata <- trophicdata[, !(colnames(trophicdata) == "" | is.na(colnames(trophicdata)))]
trophicdata <- trophicdata %>% mutate(Trophic = ifelse(name == "Spumella", "Consumer", Trophic))
trophicdata_clean <- trophicdata %>% distinct(name, .keep_all = TRUE)
trophicdata_clean$name <- trimws(trophicdata_clean$name)

# Load network correlations
setwd("C:/Users/User/OneDrive/Desktop/qian/R/network/18cs-cyc-ecoli")
ecoli_r <- read.csv("R.ecoli_level8_sig.csv", sep = ",", header = TRUE)
colnames(ecoli_r) <- trimws(colnames(ecoli_r))
ecoli_r <- ecoli_r %>% dplyr::select(X, `E.coli`) %>% filter(X != "E.coli")
colnames(ecoli_r)[1] <- "name"; colnames(ecoli_r)[2] <- "rvalue"
ecoli_r$name <- gsub("st__un_d[0-9]+__", "un_", ecoli_r$name)
ecoli_r$name <- gsub("d[0-9]+__", "", ecoli_r$name)
ecoli_r_less0 <- ecoli_r %>% filter(rvalue < 0)

# Merge network and trophic classification, keep consumers
merged_trophic_network_r_less0 <- ecoli_r_less0 %>%
  left_join(trophicdata_clean %>% select(name, Trophic), by = "name") %>%
  distinct() %>% filter(Trophic == "Consumer")

# Set output directory
output_dir_r_less0 <- "C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/key protists RA change_scatterplot_18cs/network-cyc_r_less0/Consumer"

dir.create(output_dir_r_less0, showWarnings = FALSE, recursive = TRUE)

# Extract relevant taxa
for (i in 1:nrow(merged_trophic_network_r_less0)) assign(paste0("taxa_less0_", i), merged_trophic_network_r_less0$name[i])

# Extract data for plotting (RA >= 1%)
for (i in 1:nrow(merged_trophic_network_r_less0)) {
  data <- RA_data %>% filter(group == get(paste0("taxa_less0_", i))) %>% select(-group, -Sample)
  if (max(data$RelativeAbundance, na.rm = TRUE) >= 1) assign(paste0("taxa_less0_", i, "_data"), data)
}
remaining_taxa_indices_r_less0 <- which(sapply(1:nrow(merged_trophic_network_r_less0), function(i) exists(paste0("taxa_less0_", i, "_data"))))

# Add simulated Day0 +CHX for plotting
for (i in remaining_taxa_indices_r_less0) {
  original_data <- get(paste0("taxa_less0_", i, "_data"))
  new_rows <- original_data %>% filter(Time == "Day0", Treatment == "-CHX") %>% mutate(Treatment = "+CHX")
  assign(paste0("updated_taxa_less0_", i, "_data"), bind_rows(original_data, new_rows))
}

# Convert Time to factor
for (i in remaining_taxa_indices_r_less0) {
  df <- get(paste0("updated_taxa_less0_", i, "_data"))
  df$Time <- factor(df$Time, levels = c("Day0", "Day14", "Day28"))
  assign(paste0("updated_taxa_less0_", i, "_data"), df)
}

# Generate plots
num_cols <- 3
plot_list <- list()
for (i in remaining_taxa_indices_r_less0) {
  data_to_plot <- get(paste0("updated_taxa_less0_", i, "_data"))
  plot_title <- get(paste0("taxa_less0_", i))
  ylim_max <- max(data_to_plot$RelativeAbundance, na.rm = TRUE) * 1.1
  is_leftmost <- ((which(remaining_taxa_indices_r_less0 == i) - 1) %% num_cols) == 0
  
  plot_for_saving <- ggplot(data_to_plot, aes(x = Treatment, y = RelativeAbundance, color = Time)) +
    geom_boxplot(position = position_dodge(0.6), width = 0.5, size = 0.8, outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6), size = 3, alpha = 0.7) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4, size = 0.4, position = position_dodge(0.5), show.legend = FALSE) +
    labs(title = plot_title, x = "", y = "Relative abundance (%)") +
    scale_y_continuous(labels = number_format(accuracy = 0.1), limits = c(0, ylim_max)) +
    scale_color_manual(values = c("Day0" = "#3A3737", "Day14" = "#1f77b4", "Day28" = "#ffcc00")) +
    theme_classic() +
    theme(axis.text = element_text(size = 20, color = "black", family = "Arial"),
          axis.title = element_text(size = 20, color = "black", family = "Arial"),
          plot.title = element_text(size = 20, hjust = 0.5, family = "Arial"),
          legend.title = element_blank(), legend.position = "bottom")
  
  ggsave(paste0(output_dir_r_less0, "/taxa_less0_", i, "_", gsub("[^a-zA-Z0-9]", "_", plot_title), ".jpg"),
         plot = plot_for_saving, device = "jpg", width = 9, height = 8, units = "cm", dpi = 600)
  
  plot_for_merging <- plot_for_saving + labs(y = ifelse(is_leftmost, "Relative abundance (%)", ""))
  plot_list[[length(plot_list) + 1]] <- plot_for_merging
}

combined_plot_r_less0 <- wrap_plots(plotlist = plot_list, ncol = num_cols) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom") +
  plot_annotation(title = "Conventional soil", theme = theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Arial")))

print(combined_plot_r_less0)
ggsave(file.path(output_dir_r_less0, "combined_taxa_plots_r_less0-box.jpg"),
       plot = combined_plot_r_less0, device = "jpg", width = 28, height = 22, units = "cm", dpi = 600)
