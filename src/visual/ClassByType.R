suppressPackageStartupMessages({
  library(circlize)
  library(readr)
  library(dplyr)
  library(RColorBrewer)
})

ClassByType <- function(input_file, output_file) {
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  circ <- read_delim(input_file, delim = "\t", show_col_types = FALSE)
  if (nrow(circ) == 0) {
    stop("Input file has no data or cannot be read.")
  } else {
    print("Data loaded successfully.")
  }
  
  circ_data <- circ %>%
    select(chr, start, end, circRNA_type)
  circ_data$chr <- factor(circ_data$chr, levels = unique(circ_data$chr))
  chr_limits <- circ_data %>%
    group_by(chr) %>%
    summarize(start = min(start), end = max(end))
  print("Chromosome limits calculated:")
  print(chr_limits)
  
  chr_limits <- chr_limits %>%
    filter(!is.na(start) & !is.na(end))
  chr_limits <- arrange(chr_limits, chr)
  circ_data <- arrange(circ_data, chr)
  chromosomes <- chr_limits$chr
  xlim_matrix <- as.matrix(chr_limits[, c("start", "end")])
  
  print("Chromosome names and order should match with xlim:")
  print(chromosomes)
  print(paste("Number of rows in xlim_matrix:", nrow(xlim_matrix)))

  print(paste("Creating PDF at:", output_file))
  pdf(output_file, width = 10, height = 10)
  circos.par(cell.padding = c(0.02, 0.02, 0.02, 0.02))
  circos.par(gap.after = rep(5, length(unique(circ_data$chr))))
  circos.clear()
  circos.genomicInitialize(circ_data)
  
  chromosomes <- unique(circ_data$chr)
  num_chromosomes <- length(chromosomes)
  background_colors <- brewer.pal(min(num_chromosomes, 12), "Set3")
  if (num_chromosomes > 12) {
    background_colors <- colorRampPalette(brewer.pal(12, "Set3"))(num_chromosomes)
  }
  
  circos.track(
    ylim = c(0, 1), 
    bg.col = background_colors,
    bg.border = NA,
    track.height = 0.05
  )
  
  color_map <- c(
    "eig-circRNA" = "#5BC0EB",
    "iig-circRNA" = "#FDE74C",
    "dg-circRNA" = "#9BC53D",
    "e-circRNA" = "#E55934",
    "i-circRNA" = "#FA7921",
    "ei-circRNA" = "#3CC453",
    "ig-circRNA" = "#197EEA"
  )
  
  for (type in names(color_map)) {
    circ_data_tmp <- circ_data[circ_data$circRNA_type == type, ]
    if (nrow(circ_data_tmp) > 0) {
      print(paste("Adding circRNA type:", type))
      circos.genomicDensity(
        circ_data_tmp,
        col = color_map[type], 
        track.height = 0.05,
        bg.border = NA
      )
    }
  }
  
  circos.trackPlotRegion(
    factors = circ_data$chr,
    y = rep(1, nrow(circ_data)),
    ylim = c(0, 1),
    bg.col = "lightgrey",
    bg.border = NA,
    track.height = 0.1,
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      circos.text(CELL_META$xcenter, CELL_META$ycenter, chr,
                  facing = "inside", niceFacing = TRUE, cex = 0.6)
    }
  )
  
  legend("topright", legend = names(color_map), fill = color_map, title = "circRNA Type", cex = 0.7)
  dev.off()
  print("PDF file created successfully.")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript classByType.R <input_file> <output_file>")
}

input_file <- args[1]
output_file <- args[2]
ClassByType(input_file, output_file)
