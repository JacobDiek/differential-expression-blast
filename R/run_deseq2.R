# Generic DESeq2 step of the pipeline. run_pipeline.py runs this with a JSON
# file of settings taken from pipeline_config.yaml:
#   Rscript R/run_deseq2.R settings.json

suppressPackageStartupMessages({
  library("DESeq2")
  library("ggplot2")
  library("jsonlite")
})

args <- commandArgs(trailingOnly = TRUE)
settings <- fromJSON(args[1], simplifyVector = FALSE)

results_dir <- settings$output_dir
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Make cts table

files <- list.files(path = settings$counts_dir, pattern = "ReadsPerGene\\.out\\.tab$", full.names = TRUE)

sample_sheet <- read.csv(settings$sample_sheet, stringsAsFactors = FALSE)
sample_sheet$sample <- trimws(sample_sheet$sample)
sample_sheet$condition <- trimws(sample_sheet$condition)

exclude_samples <- unlist(settings$exclude_samples)
sample_sheet <- sample_sheet[!(sample_sheet$sample %in% exclude_samples), ]

# Each sample's count file starts with the sample name followed by "_" or "."
find_count_file <- function(sample) {
  matches <- files[startsWith(basename(files), paste0(sample, "_")) |
                     startsWith(basename(files), paste0(sample, "."))]
  if (length(matches) != 1) {
    stop("Expected 1 count file for sample ", sample, ", found ", length(matches))
  }
  matches
}

sample_files <- sapply(sample_sheet$sample, find_count_file)

# STAR ReadsPerGene columns: 2 = unstranded, 3 = forward, 4 = reverse stranded
read_star <- function(file) {
  read.table(file, skip = 4, header = FALSE, stringsAsFactors = FALSE)
}

star_tables <- lapply(sample_files, read_star)
gene_ids <- star_tables[[1]][, 1]

for (table in star_tables) {
  stopifnot(identical(table[, 1], gene_ids))
}

cts <- do.call(cbind, lapply(star_tables, function(table) table[, settings$star_count_column]))
rownames(cts) <- gene_ids
colnames(cts) <- sample_sheet$sample

cts <- cts[rownames(cts) != "MissingGeneID", ]

# Make col data

coldata <- data.frame(
  row.names = sample_sheet$sample,
  condition = factor(sample_sheet$condition)
)

# Running DESeq2

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)

# PCA, for spotting outlier samples to add to exclude_samples

if (nrow(dds) >= 1000) {
  vsd <- vst(dds, blind = TRUE)
} else {
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
}
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
write.csv(pca_data[, c("name", "condition", "PC1", "PC2")],
          file.path(results_dir, "pca_coordinates.csv"), row.names = FALSE)
ggsave(file.path(results_dir, "pca.png"), plotPCA(vsd, intgroup = "condition") + theme_minimal(),
       width = 7, height = 6, dpi = 300)

make_volcano_plot <- function(res_df, title, padj_cutoff, lfc_cutoff) {

  res_df$significant <- !is.na(res_df$padj) &
    res_df$padj < padj_cutoff &
    abs(res_df$log2FoldChange) > lfc_cutoff

  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significant), alpha = 0.4, size = 1) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey40") +
    labs(title = title, x = "log2 fold change", y = "-log10 adjusted p-value") +
    theme_minimal() +
    theme(legend.position = "none")
}

# Extract each comparison

sig_gene_counts <- list()

for (contrast in settings$contrasts) {
  condition <- contrast[[1]]
  reference <- contrast[[2]]
  comparison_name <- paste0(condition, "_vs_", reference)

  res <- results(dds, contrast = c("condition", condition, reference))

  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]

  write.csv(res_df, file.path(results_dir, paste0("all_genes_", comparison_name, ".csv")), row.names = FALSE)

  sig_genes <- res_df[!is.na(res_df$padj) &
                        res_df$padj < settings$padj_cutoff &
                        abs(res_df$log2FoldChange) > settings$log2fc_cutoff, ]
  sig_genes <- sig_genes[order(sig_genes$padj), ]

  sig_gene_counts[[comparison_name]] <- nrow(sig_genes)

  write.csv(sig_genes, file.path(results_dir, paste0("significant_genes_", comparison_name, ".csv")), row.names = FALSE)

  volcano <- make_volcano_plot(res_df, gsub("_", " ", comparison_name),
                               settings$padj_cutoff, settings$log2fc_cutoff)
  ggsave(file.path(results_dir, paste0("volcano_", comparison_name, ".png")), volcano,
         width = 7, height = 6, dpi = 300)
}

print(sig_gene_counts)

# Versions used, for the pipeline's run record
write_json(list(R = R.version.string, DESeq2 = as.character(packageVersion("DESeq2")),
                ggplot2 = as.character(packageVersion("ggplot2"))),
           file.path(results_dir, "versions.json"), auto_unbox = TRUE, pretty = TRUE)
