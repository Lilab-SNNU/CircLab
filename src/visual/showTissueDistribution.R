library(ggplot2)
library(dplyr)

showTissueDistribution <- function(file_path, output_file) {
  # check
  if (!file.exists(file_path)) {
    stop("The file is not exist：", file_path)
  }
  bed <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (!("tissue" %in% colnames(bed)) || !("circRNA_type" %in% colnames(bed))) {
    stop("The file lack of 'tissue' or 'circRNA_type' column")
  }
  tissue_counts <- bed %>%
    group_by(tissue, circRNA_type) %>%
    summarise(count = n(), .groups = 'drop')  
  colors <- c("#91C9A4", "#E5755B", "#C0F1D3", "#9AC6D9", "#E5A456", "#7C97C7", "#FAECBC", "#A38DBD")
  p <- ggplot(tissue_counts, aes(x = tissue, y = count, fill = circRNA_type)) +
    geom_bar(stat = "identity") +
    labs(x = "Tissue Type", y = "Number of circRNA", title = "Distribution of circRNA in Tissues") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
    scale_fill_manual(values = colors) 
  print(p)
  ggsave(output_file, plot = p, width = 7, height = 7, bg = "white")
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript showTissueDistribution.R <input_file> <output_file>")
}
input_file <- args[1]
output_file <- args[2]

showTissueDistribution(input_file, output_file)