# showCodon.R
suppressPackageStartupMessages({
  library(ggplot2)
  library(Biostrings)
  library(dplyr)
  library(tibble)
})

showCodon <- function(input_file, max_patterns = 20, output_file = "codon_distribution.png") {
  circ <- readDNAStringSet(input_file)
  codon_matrix <- trinucleotideFrequency(circ, step = 1)
  codon_types <- colnames(codon_matrix)
  codon_totals <- colSums(codon_matrix)
  codon_df <- data.frame(Freq = as.vector(codon_totals), type = codon_types)
  codon_df <- codon_df[order(-codon_df$Freq), ]
  codon_df <- head(codon_df, max_patterns)
  max_freq <- max(codon_df$Freq)
  
  ggplot(codon_df) +
    geom_bar(stat = "identity", aes(x = type, y = Freq, fill = Freq)) +
    scale_fill_gradient(low = "#ecf8f3", high = "#097b66fa") +
    coord_flip() +
    scale_x_discrete(limits = rev(codon_df$type)) +
    labs(title = paste("Top", max_patterns, "most frequent 3-nt patterns"), x = "3-nt Pattern", y = "Count") +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.background = element_blank()) +
    scale_y_continuous(limits = c(0, max_freq * 1.1), expand = c(0, 0)) +  
    geom_text(aes(x = type, y = Freq + max_freq * 0.05, label = Freq), size = 3, color = "black", hjust = 0.5)
  
  ggsave(output_file, width = 9, height = 6)
}

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
max_patterns <- as.integer(args[2])
output_file <- args[3]

showCodon(input_file, max_patterns, output_file)
