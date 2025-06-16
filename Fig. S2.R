library(ggplot2)

# Sample data (make sure `count` is correctly named)
data <- data.frame(
  group = c("STWW", "STWW+E.coli"),
  count = c(0, 8.46)
)

# Create the bar plot
a <- ggplot(data, aes(x = group, y = count, fill = group, color = group)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge(width = 0.2)) +  # Explicitly set stat = "identity"
  scale_fill_manual(values = c("STWW" = "#9AC9DB", "STWW+E.coli" = "#40588F")) +  # Fix fill color mapping
  scale_color_manual(values = c("STWW" = "#9AC9DB", "STWW+E.coli" = "#40588F")) + 
  theme_classic() +
  # Fix border color mapping
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 20, color = "black", family = "Arial"),
    axis.text.y = element_text(size = 20, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 20, color = "black", family = "Arial"),
    axis.title.x = element_text(size = 20, color = "black", family = "Arial"),
    plot.title = element_text(size = 20, hjust = 0.5, family = "Arial")
  ) +
  labs(title = "", x = "", y = expression(atop(italic("E. coli")~"concentration", "(log10 CFU/mL)")))


# Print the plot
print(a)

ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/ecoli culture in STWW.tiff", plot =a, device = "tiff", dpi = 600, width = 16, height = 16, units = "cm")
ggsave("C:/Users/User/OneDrive/Desktop/qian/R/paper1_figs/ecoli culture in STWW.jpg", plot =a, device = "jpeg", dpi = 600, width = 16, height = 16, units = "cm")

