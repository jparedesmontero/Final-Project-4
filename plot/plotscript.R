# Load libraries

library(ggplot2)
library(readr)

# Read the file
hotspots <- read_csv("mutation_hotspot_counts.csv")

# Define gene schematic
gene_start <- 0
gene_end <- 1100
gene_y <- -2  # Y-position to draw the gene schematic

# Plot with gene schematic and lollipops
ggplot(hotspots, aes(x = Position, y = Frequency)) +
  # Gene schematic
  geom_segment(aes(x = gene_start, xend = gene_end, y = gene_y, yend = gene_y),
               size = 6, color = "gray80") +
  
  # Lollipop stems
  geom_segment(aes(xend = Position, yend = gene_y), color = "gray40", linewidth = 0.7) +
  
  # Lollipop heads
  geom_point(color = "firebrick", size = 3) +
  
  # Axis and labels
  labs(
    title = "Mutation Hotspots in the Beta-lactamase Gene",
    x = "Position in Gene (bp)",
    y = "Number of Samples with SNP"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  ) +
  ylim(gene_y - 1, max(hotspots$Frequency) + 2)


metadata <- read.csv("metadata.csv")
head(metadata)
snp_counts <- read.csv("snp_counts.csv")
head(snp_counts)
snp_counts$Run <- sub("\\.sorted$", "", snp_counts$Sample)

merged <- merge(snp_counts, metadata, by = "Run")
head(merged)

merged <- merged[, !names(merged) %in% "Sample"]

# Filter out samples with 0 SNPs
filtered <- merged %>% filter(SNP_Count > 0)

# Plot
ggplot(filtered, aes(x = reorder(Run, -SNP_Count), y = SNP_Count, fill = isolation_source)) +
  geom_col() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, size = 12),  # 🔎 Larger x labels
    axis.text.y = element_text(size = 12),              # 🧭 Optional: enlarge y labels too
    axis.title = element_text(size = 14),               # 📐 Axis titles
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  labs(
    title = "Mutational load in betalactamase gene",
    x = "Sample",
    y = "Number of SNPs"
  )

library(ggplot2)
library(dplyr)
library(scales)

# Summarize the SNP status
merged$SNP_Status <- ifelse(merged$SNP_Count > 0, "With SNPs", "No SNPs")

snp_summary <- merged %>%
  count(SNP_Status) %>%
  mutate(
    percent = Count / sum(Count),
    label = paste0(SNP_Status, "\n", percent(percent))
  )

# Make the pie chart
ggplot(snp_summary, aes(x = "", y = Count, fill = SNP_Status)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            size = 6, fontface = "bold", color = "white") +
  scale_fill_manual(values = c("With SNPs" = "#D73027", "No SNPs" = "#999999")) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.position = "none"
  ) +
  labs(title = "Proportion of Samples With and Without SNPs")


library(ggplot2)
library(dplyr)
library(scales)

# Prepare SNP status summary
merged$SNP_Status <- ifelse(merged$SNP_Count > 0, "With SNPs", "No SNPs")

snp_summary <- merged %>%
  count(SNP_Status) %>%
  mutate(
    percent = Count / sum(Count),
    label = paste0(Count, " (", percent(percent), ")")
  )

# Create stacked bar
ggplot(snp_summary, aes(x = "Samples", y = Count, fill = SNP_Status)) +
  geom_bar(stat = "identity", width = 0.6, color = "white") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            color = "white", size = 6, fontface = "bold") +
  scale_fill_manual(values = c("With SNPs" = "#D73027", "No SNPs" = "#999999")) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.position = "none"
  ) +
  labs(title = "Distribution of Samples by SNP Presence")
