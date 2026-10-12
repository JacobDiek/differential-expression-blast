library("DESeq2")
library("apeglm")
library("ashr")
library("ggplot2")
library("ggrepel")

# Run from the project root folder, e.g. setwd("path/to/differential-expression-blast")
# All results (CSVs and plots) are written to results/brook_analysis/

results_dir <- file.path("results", "brook_analysis")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Make cts table

path <- file.path("data", "counts", "brook_reads_updated")

files <- list.files(path = path, pattern = "\\.out\\.tab$", full.names = TRUE)

sample_names <- gsub("_S\\d+_L\\d+_ReadsPerGene\\.out\\.tab$", "", basename(files))
names(files) <- sample_names

# STAR ReadsPerGene column 4 = reverse-stranded counts
read_star <- function(file) {
  read.table(file, skip = 4, header = FALSE)[, 4]
}

gene_ids <- read.table(files[1], skip = 4, header = FALSE)[, 1]

cts <- do.call(cbind, lapply(files, read_star))
rownames(cts) <- gene_ids
colnames(cts) <- sample_names

cts <- cts[rownames(cts) != "MissingGeneID", ]

metadata <- read.csv(file.path("data", "25030-01-final-QC-summary.csv"))

# Line metadata rows up with the count columns by sample name
metadata <- metadata[match(colnames(cts), metadata$Sample), ]
stopifnot(identical(metadata$Sample, colnames(cts)))

# Make col data

coldata <- data.frame(
  row.names = colnames(cts),
  customer_id = metadata$CustomerID,
  tissue = factor(ifelse(grepl("^LF", metadata$CustomerID), "liver_female",
                         ifelse(grepl("^LM", metadata$CustomerID), "liver_male",
                                ifelse(grepl("^T", metadata$CustomerID), "testes", "ovary"))))
)

# Running DESeq2

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ tissue)
dds <- DESeq(dds)

# set outlier_sample below and re-run from here if needed

vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "tissue", returnData = TRUE)
print(pca_data[, c("name", "tissue", "PC1", "PC2")])
plotPCA(vsd, intgroup = "tissue")

outlier_sample <- "25030R-01-02"

if (!is.null(outlier_sample)) {
  cts <- cts[, colnames(cts) != outlier_sample]
  coldata <- coldata[rownames(coldata) != outlier_sample, , drop = FALSE]
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ tissue)
  dds <- DESeq(dds)
}

# Extract sex-differential comparisons: testes vs. ovary (gonad vs. gonad,
# sex-matched tissue) and liver_male vs. liver_female (same tissue, sex only)

contrasts_to_run <- list(
  testes_vs_ovary            = c("tissue", "testes", "ovary"),
  liver_male_vs_liver_female = c("tissue", "liver_male", "liver_female")
)

print(resultsNames(dds))

sig_gene_counts <- list()

for (comparison_name in names(contrasts_to_run)) {
  
  res <- results(dds, contrast = contrasts_to_run[[comparison_name]])
  
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  sig_genes <- res_df[!is.na(res_df$padj) &
                        res_df$padj < 0.05 &
                        abs(res_df$log2FoldChange) > 1, ]
  sig_genes <- sig_genes[order(sig_genes$padj), ]
  sig_genes <- sig_genes[, c("gene_id", setdiff(colnames(sig_genes), "gene_id"))]
  
  sig_gene_counts[[comparison_name]] <- nrow(sig_genes)
  
  write.csv(sig_genes, file.path(results_dir, paste0("significant_genes_", comparison_name, ".csv")), row.names = FALSE)
}

print(sig_gene_counts)

# Plot making

# Genes of interest to label on volcano plots, by comparison
genes_of_interest_liver <- c(
  egapxtmp_003028 = "estrogen receptor 1",
  egapxtmp_025941 = "vitellogenin-2-like",
  egapxtmp_027432 = "vitellogenin 3",
  egapxtmp_018468 = "zona pellucida 3-like",
  egapxtmp_004929 = "zona pellucida 3-like",
  egapxtmp_028631 = "nothepsin",
  egapxtmp_009392 = "cytosolic sulfotransferase 3-like"
)

genes_of_interest_gonad <- c(
  egapxtmp_033296 = "insulin-like 3",
  egapxtmp_015645 = "cyp11c1",
  egapxtmp_014751 = "SUN domain-containing 3-like",
  egapxtmp_002091 = "claudin-3-like",
  egapxtmp_033281 = "zona pellucida 3-like",
  egapxtmp_031446 = "high choriolytic enzyme 1",
  egapxtmp_018019 = "protein Bouncer-like"
)

make_volcano_plot <- function(res_df, genes_of_interest, title, padj_cutoff = 0.05, lfc_cutoff = 1) {
  
  res_df$significant <- !is.na(res_df$padj) &
    res_df$padj < padj_cutoff &
    abs(res_df$log2FoldChange) > lfc_cutoff
  
  res_df$label <- genes_of_interest[res_df$gene_id]
  
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significant), alpha = 0.4, size = 1) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey40") +
    geom_point(data = subset(res_df, !is.na(label)), color = "blue", size = 2) +
    geom_text_repel(data = subset(res_df, !is.na(label)), aes(label = label),
                    size = 3.5, max.overlaps = Inf, box.padding = 0.5) +
    labs(title = title, x = "log2 fold change", y = "-log10 adjusted p-value") +
    theme_minimal() +
    theme(legend.position = "none")
}

make_ma_plot <- function(shrunk_res, title) {
  res_df <- as.data.frame(shrunk_res)
  res_df$significant <- !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
  
  ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
    geom_point(aes(color = significant), alpha = 0.4, size = 1) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    labs(title = title, x = "log10 mean expression", y = "shrunken log2 fold change") +
    theme_minimal() +
    theme(legend.position = "none")
}

# --- Testes vs. Ovary ---

res_gonad <- results(dds, contrast = c("tissue", "testes", "ovary"))
res_gonad_df <- as.data.frame(res_gonad)
res_gonad_df$gene_id <- rownames(res_gonad_df)

volcano_gonad <- make_volcano_plot(res_gonad_df, genes_of_interest_gonad, "Testes vs. Ovary")
ggsave(file.path(results_dir, "volcano_testes_vs_ovary.png"), volcano_gonad, width = 7, height = 6, dpi = 300)

shrunk_gonad <- lfcShrink(dds, contrast = c("tissue", "testes", "ovary"), type = "ashr")
ma_gonad <- make_ma_plot(shrunk_gonad, "Testes vs. Ovary")
ggsave(file.path(results_dir, "ma_testes_vs_ovary.png"), ma_gonad, width = 7, height = 6, dpi = 300)

# --- Liver Male vs. Liver Female ---

res_liver <- results(dds, contrast = c("tissue", "liver_male", "liver_female"))
res_liver_df <- as.data.frame(res_liver)
res_liver_df$gene_id <- rownames(res_liver_df)

volcano_liver <- make_volcano_plot(res_liver_df, genes_of_interest_liver, "Liver: Male vs. Female")
ggsave(file.path(results_dir, "volcano_liver_male_vs_liver_female.png"), volcano_liver, width = 7, height = 6, dpi = 300)

shrunk_liver <- lfcShrink(dds, contrast = c("tissue", "liver_male", "liver_female"), type = "ashr")
ma_liver <- make_ma_plot(shrunk_liver, "Liver: Male vs. Female")
ggsave(file.path(results_dir, "ma_liver_male_vs_liver_female.png"), ma_liver, width = 7, height = 6, dpi = 300)


# AmhY egapxtmp_015627
# Wt1b egapxtmp_013970
# Cyp19a1a egapxtmp_013395
# gsdf egapxtmp_030606

wt1b_data <- plotCounts(dds,
                        gene = "egapxtmp_013970",
                        intgroup = "tissue",
                        transform = TRUE,
                        normalized = FALSE,
                        returnData = FALSE,
                        main = "Wt1b")

gsdf_data <- plotCounts(dds,
                        gene = "egapxtmp_030606",
                        intgroup = "tissue",
                        transform = TRUE,
                        normalized = FALSE,
                        returnData = FALSE,
                        main = "Gsdf")

cyp19a1a_data <- plotCounts(dds,
                            gene = "egapxtmp_013395",
                            intgroup = "tissue",
                            transform = TRUE,
                            normalized = FALSE,
                            returnData = FALSE,
                            main = "Cyp19a1a")

amhy_data <- plotCounts(dds,
                            gene = "egapxtmp_015627",
                            intgroup = "tissue",
                            transform = TRUE,
                            normalized = FALSE,
                            returnData = FALSE,
                            main = "AmhY")

# === Correlation testing ===

genes_of_interest <- c("egapxtmp_015627",
                       "egapxtmp_013970",
                       "egapxtmp_013395",
                       "egapxtmp_030606")

norm_counts <- counts(dds, normalized = TRUE)
sub_counts <- norm_counts[genes_of_interest, ]
rownames(sub_counts) <- c("AmhY", "Wt1b", "Cyp19a1a", "Gsdf")

# Correlation across all samples
cor(t(sub_counts), method = "pearson")

# Correlation within testes only
testes_samples <- rownames(coldata)[coldata$tissue == "testes"]
sub_counts_testes <- sub_counts[, testes_samples]
cor(t(sub_counts_testes), method = "pearson")

# Correlation within ovary only
ovary_samples <- rownames(coldata)[coldata$tissue == "ovary"]
sub_counts_ovary <- sub_counts[, ovary_samples]
cor(t(sub_counts_ovary), method = "pearson")
