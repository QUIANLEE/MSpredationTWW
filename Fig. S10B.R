# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

rm(list = ls())  # Clear the workspace

############################ 18-CS+cyc ###################################

# Load relative abundance table
setwd("C:/Users/User/OneDrive/Desktop/qian/R/cyc+-relative_abundance(lefse)/18dna-cs")
RA_data <- read.table("level-8.txt", sep = "\t", header = TRUE, check.names = FALSE)

# Convert to long format
RA_data <- RA_data %>%
  pivot_longer(cols = -group, names_to = "Sample", values_to = "RelativeAbundance") %>%
  mutate(RelativeAbundance = RelativeAbundance * 100,
         group = gsub("st__un_d[0-9]+__", "un_", group),
         group = gsub("d[0-9]+__", "", group)) %>%
  mutate(
    Treatment = case_when(
      grepl("^CSc-", Sample) ~ "+CHX",
      grepl("^CSn-", Sample) ~ "-CHX"
    ),
    Time = case_when(
      grepl("-14$", Sample) ~ "Day14",
      grepl("-28$", Sample) ~ "Day28",
      grepl("-0$", Sample)  ~ "Day0"
    )
  )

# Load trophic annotation
trophic_path <- "C:/Users/User/OneDrive/Desktop/qian/R/CYC-filter0.01/18dna-cs"
setwd(trophic_path)
trophicdata <- read.table("df_filtered-trophic.txt", sep = "\t", header = TRUE, check.names = FALSE)
colnames(trophicdata)[colnames(trophicdata) == "Rank8"] <- "name"
trophicdata <- trophicdata %>%
  mutate(name = gsub("st__un_d[0-9]+__", "un_", name),
         name = gsub("d[0-9]+__", "", name),
         name = trimws(name)) %>%
  distinct(name, .keep_all = TRUE)

# Load network correlation data
network_path <- "C:/Users/User/OneDrive/Desktop/qian/R/network/18cs+cyc-ecoli"
setwd(network_path)
ecoli_r <- read.csv("R.ecoli_level8_sig.csv", sep = ",", header = TRUE) %>%
  dplyr::select(X, E.coli) %>%
  filter(X != "E.coli") %>%
  rename(name = X, rvalue = E.coli)

# Clean taxa names
ecoli_r <- ecoli_r %>%
  mutate(name = gsub("st__un_d[0-9]+__", "un_", name),
         name = gsub("d[0-9]+__", "", name))

# Filter taxa negatively correlated with E. coli
ecoli_r_less0 <- ecoli_r %>% filter(rvalue < 0)

# Merge with trophic data and filter Consumers
merged <- ecoli_r_less0 %>%
  left_join(trophicdata %>% select(name, Trophic), by = "name") %>%
  filter(Trophic == "Consumer") %>%
  distinct()

# Output path
output_dir <- "C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/key protists RA change_scatterplot_18cs/network+cyc_r_less0/Consumer"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Extract and filter RA data for selected taxa
for (i in 1:nrow(merged)) {
  taxa <- merged$name[i]
  assign(paste0("taxa_less0_", i), taxa)
  data <- RA_data %>% filter(group == taxa) %>% select(-group, -Sample)
  if (max(data$RelativeAbundance, na.rm = TRUE) >= 1) {
    assign(paste0("taxa_less0_", i, "_data"), data)
  }
}

# Determine valid taxa
valid_indices <- which(sapply(1:nrow(merged), function(i) exists(paste0("taxa_less0_", i, "_data"))))

# Simulate additional +CHX for Day0
for (i in valid_indices) {
  df <- get(paste0("taxa_less0_", i, "_data"))
  dup_rows <- df %>% filter(Time == "Day0", Treatment == "-CHX") %>% mutate(Treatment = "+CHX")
  df_updated <- bind_rows(df, dup_rows)
  df_updated$Time <- factor(df_updated$Time, levels = c("Day0", "Day14", "Day28"))
  assign(paste0("updated_taxa_less0_", i, "_data"), df_updated)
}

# Plotting
num_cols <- 3
plot_list <- list()

for (i in valid_indices) {
  df <- get(paste0("updated_taxa_less0_", i, "_data"))
  title <- get(paste0("taxa_less0_", i))
  ylim_max <- max(df$RelativeAbundance, na.rm = TRUE) * 1.1
  is_left <- ((which(valid_indices == i) - 1) %% num_cols) == 0
  
  p <- ggplot(df, aes(x = Treatment, y = RelativeAbundance, color = Time)) +
    geom_boxplot(position = position_dodge(0.6), width = 0.5, size = 0.8, outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(0.15, 0.6), size = 3, alpha = 0.7) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4, size = 0.4,
                 position = position_dodge(0.5), aes(color = Time), show.legend = FALSE) +
    labs(title = title, x = "", y = ifelse(is_left, "Relative abundance (%)", "")) +
    scale_y_continuous(labels = number_format(accuracy = 0.1), limits = c(0, ylim_max)) +
    scale_color_manual(values = c("Day0" = "#3A3737", "Day14" = "#1f77b4", "Day28" = "#ffcc00")) +
    theme_classic() +
    theme(axis.text = element_text(size = 20, color = "black", family = "Arial"),
          axis.title = element_text(size = 20, color = "black", family = "Arial"),
          plot.title = element_text(size = 20, hjust = 0.5, family = "Arial"),
          legend.title = element_blank(),
          legend.position = "bottom")
  
  ggsave(file.path(output_dir, paste0("taxa_less0_", i, "_", gsub("[^a-zA-Z0-9]", "_", title), ".jpg")),
         plot = p, device = "jpg", width = 9, height = 8, units = "cm", dpi = 600)
  
  plot_list[[length(plot_list) + 1]] <- p
}

# Combine plots
combined <- wrap_plots(plotlist = plot_list, ncol = num_cols) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") +
  plot_annotation(title = "Conventional soil",
                  theme = theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Arial")))

# Save combined figure
ggsave(file.path(output_dir, "combined_taxa_plots_r_less0-box.jpg"),
       plot = combined, device = "jpg", width = 28, height = 11, units = "cm", dpi = 600)

print(combined)
