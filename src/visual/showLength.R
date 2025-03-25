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

# 定义 showLength 函数
showLength <- function(file_path, max_length) {
  if (!file.exists(file_path)) {
    stop("文件不存在，请检查路径：", file_path)
  }
  
  bed <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  if (!("start" %in% colnames(bed)) || !("end" %in% colnames(bed))) {
    stop("文件中缺少 'start' 或 'end' 列")
  }
  
  len <- bed$end - bed$start
  less_max <- len[len < max_length]
  
  if (length(less_max) == 0) {
    stop("没有长度小于 ", max_length, " 的 circRNA")
  }
  
  m <- mean(less_max)
  ave <- gsub("%", floor(m), "Average length: % bp")
  ti <- gsub("%", max_length, "length\n(less than %)")
  
  print(paste("共有", length(less_max), "条 circRNA，平均长度为", floor(m), "bp"))
  
  # 生成第一个图 (带密度线)
  g1 <- ggplot(as.data.frame(less_max), aes(x = less_max, y = after_stat(density))) +
    geom_histogram(binwidth = max_length / 30, colour = "black", fill = "#a1dcc1") +
    geom_line(stat = "density", colour = "#660066") +
    expand_limits(y = 0) +
    geom_vline(aes(xintercept = m), linetype = 5, colour = "darkred") +
    labs(x = ti, title = "Distribution of circRNA length") +
    theme_minimal()
  
  
  g2 <- ggplot(as.data.frame(less_max), aes(x = less_max)) +
    labs(x = ti, title = "Distribution of circRNA length") +
    geom_histogram(binwidth = max_length / 30, colour = "black", fill = "#a1dcc1") +
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
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript showLength_xm.R <input_file> <max_length> <output_file>")
}

input_file <- args[1]
max_length <- as.numeric(args[2])
output_file <- args[3]

png(output_file)

showLength(input_file, max_length)

dev.off()
