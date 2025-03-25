options(repos = c(CRAN = "https://cran.r-project.org"))

required_list <- c(
  'ggplot2',
  'gtable',
  'grid'
)

install_if_missing <- function(packages) {
  installed_list <- installed.packages()
  
  for (pkg in packages) {
    if (!(pkg %in% installed_list[, "Package"])) {
      print(paste(pkg, 'is not installed. Installing now...'))
      install.packages(pkg)
    } else {
      print(paste(pkg, 'has been installed'))
    }
  }
}

install_if_missing(required_list)
lapply(required_list, library, character.only = TRUE)

showLength2 <- function(file_path, max_length) {
  if (!file.exists(file_path)) {
    stop("The file is not exist, Please check", file_path)
  }
  bed <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  if (!("start" %in% colnames(bed)) || !("end" %in% colnames(bed))) {
    stop("The file is lack of 'start' or 'end' column")
  }
  len <- bed$end - bed$start
  less_max <- len[len < max_length]
  
  if (length(less_max) == 0) {
    stop("CircRNA length less than ", max_length, " not find")
  }
  
  m <- mean(less_max)
  ave <- gsub("%", floor(m), "Average length: % bp")
  ti <- gsub("%", max_length, "length\n(less than %)")
  
  print(paste("Total", length(less_max), "circRNA，average length", floor(m), "bp"))
  
  g1 <- ggplot(as.data.frame(less_max), aes(x = less_max, y = ..density..)) +
    geom_histogram(binwidth = max_length / 30, colour = "black", fill = "lightblue") +
    geom_density(colour = "#660066") +
    geom_vline(aes(xintercept = m), linetype = 5, colour = "darkred") +
    labs(x = ti, title = "Distribution of circRNA length") +
    theme(plot.title = element_text(hjust = 0.5))
  
  g2 <- ggplot(as.data.frame(less_max), aes(x = less_max)) +
    geom_histogram(binwidth = max_length / 30, colour = "black", fill = "lightblue") +
    labs(x = ti, title = "Distribution of circRNA length") +
    theme(plot.title = element_text(hjust = 0.5))
  
  g3 <- ggplot(as.data.frame(less_max), aes(y = less_max)) +
    geom_boxplot(fill = "#b3e3cd") +
    labs(y = "CircRNA Length", title = "Boxplot of CircRNA Length") +
    theme_minimal()
  
  ggplot2.two_y_axis <- function(g1, g2) {
    g1 <- ggplotGrob(g1)
    g2 <- ggplotGrob(g2)
    pp <- c(subset(g1$layout, name == "panel", se = t:r))
    g1 <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)
    g1 <- gtable_add_cols(g1, g2$widths[g2$layout[g2$layout$name == "ylab-l", ]$l], pp$r)
    
    grid.newpage()
    grid.draw(g1)
  }
  
  ggplot2.two_y_axis(g1, g2)
  print(g3)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript showLength_barlot.R <input_file> <max_length> <output_file>")
}

input_file <- args[1]
max_length <- as.numeric(args[2])
output_file <- args[3]
png(output_file)
showLength2(input_file, max_length)
dev.off()
