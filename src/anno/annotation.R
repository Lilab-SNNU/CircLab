#!/usr/bin/env Rscript
# annotation.R - 自动化功能注释脚本

# 自动安装缺失的包
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing missing package:", pkg))
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE)
    } else {
      install.packages(pkg, repos = 'http://cran.rstudio.com/')
    }
  }
}

# 检查并安装必要的R包
install_if_missing("argparse")
install_if_missing("ggplot2")
install_if_missing("dplyr")
install_if_missing("clusterProfiler", bioc = TRUE)
install_if_missing("org.Hs.eg.db", bioc = TRUE)
install_if_missing("org.Mm.eg.db", bioc = TRUE)
install_if_missing("org.Rn.eg.db", bioc = TRUE)
install_if_missing("org.At.tair.db", bioc = TRUE)
install_if_missing("org.Osativa.eg.db", bioc = TRUE)
install_if_missing("org.Zmays.eg.db", bioc = TRUE)
install_if_missing("org.Dm.eg.db", bioc = TRUE)
install_if_missing("org.Ce.eg.db", bioc = TRUE)

suppressPackageStartupMessages({
  library(argparse)
  library(ggplot2)
  library(dplyr)
  library(clusterProfiler)
})

# 解析命令行参数
parser <- ArgumentParser(description = "Automatic Functional Annotation")
parser$add_argument("--input", type = "character", required = TRUE)
parser$add_argument("--species", type = "integer", required = TRUE)  # 物种编号
parser$add_argument("--output", type = "character", required = TRUE)
args <- parser$parse_args()

# 物种编号对应的数据库和 KEGG 代码
species_list <- list(
  "1" = list("OrgDb" = "org.Hs.eg.db", "KEGG" = "hsa"),  # 人类
  "2" = list("OrgDb" = "org.Mm.eg.db", "KEGG" = "mmu"),  # 小鼠
  "3" = list("OrgDb" = "org.Rn.eg.db", "KEGG" = "rno"),  # 大鼠
  "4" = list("OrgDb" = "org.At.tair.db", "KEGG" = "ath"),  # 拟南芥
  "5" = list("OrgDb" = "org.Osativa.eg.db", "KEGG" = "osa"),  # 水稻
  "6" = list("OrgDb" = "org.Zmays.eg.db", "KEGG" = "zma"),  # 玉米
  "7" = list("OrgDb" = "org.Dm.eg.db", "KEGG" = "dme"),  # 果蝇
  "8" = list("OrgDb" = "org.Ce.eg.db", "KEGG" = "cel")   # 线虫
)

# 验证物种输入
species <- as.character(args$species)
if (!species %in% names(species_list)) {
  stop("ERROR: Invalid species selection. Please select a number between 1-8.")
}

orgdb <- species_list[[species]]$OrgDb
kegg <- species_list[[species]]$KEGG

suppressPackageStartupMessages({
  library(orgdb, character.only = TRUE)
})

circ_data <- read.delim(args$input, header = TRUE, sep = "\t")

# 检查是否包含 circRNA_name 和 geneName 列
if (!all(c("circRNA_name", "geneName") %in% colnames(circ_data))) {
  stop("ERROR: The input file must contain 'circRNA_name' and 'geneName' columns!")
}

# 提取基因名称并去重
genes <- unique(circ_data$geneName)

# GO 富集分析
go_res <- tryCatch({
  enrichGO(
    gene = genes,
    OrgDb = get(orgdb),
    keyType = "SYMBOL",
    ont = "ALL",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
}, error = function(e) {
  warning("GO enrichment failed: ", conditionMessage(e))
  NULL
})

# KEGG 富集分析
entrez_ids <- tryCatch({
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = get(orgdb))$ENTREZID
}, error = function(e) {
  warning("Failed to convert SYMBOL to ENTREZID: ", conditionMessage(e))
  character(0)
})

kegg_res <- tryCatch({
  if (length(entrez_ids) > 0) {
    enrichKEGG(
      gene = entrez_ids,
      organism = kegg,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
  } else {
    NULL
  }
}, error = function(e) {
  warning("KEGG enrichment failed: ", conditionMessage(e))
  NULL
})

# 创建输出目录
dir.create(args$output, showWarnings = FALSE, recursive = TRUE)

# 保存结果
if (!is.null(go_res)) {
  write.table(go_res@result, file.path(args$output, "GO_results.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

if (!is.null(kegg_res)) {
  write.table(kegg_res@result, file.path(args$output, "KEGG_results.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# 生成可视化报告
pdf(file.path(args$output, "Enrichment_Report.pdf"), width = 10, height = 8)

if (!is.null(go_res)) {
  print(
    dotplot(go_res, showCategory = 15, split = "ONTOLOGY") +
      facet_grid(ONTOLOGY ~ ., scales = "free") +
      ggtitle("GO Enrichment Analysis")
  )
}

if (!is.null(kegg_res)) {
  print(
    dotplot(kegg_res, showCategory = 15) +
      ggtitle("KEGG Pathway Enrichment")
  )
}

dev.off()

message("Annotation completed successfully")