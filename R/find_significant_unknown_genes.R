library("DESeq2")
library("apeglm")

# Make cts table

path <- "/Users/jacobdiek/Coding/projects/differential-expression-blast/R/brook_reads_updated"

files <- list.files(path = path, pattern = "\\.out\\.tab$", full.names = TRUE)

sample_names <- gsub("_S\\d+_L\\d+_ReadsPerGene\\.out\\.tab$", "", basename(files))
names(files) <- sample_names

read_star <- function(file) {
  read.table(file, skip = 4, header = FALSE)[, 4]
}

gene_ids <- read.table(files[1], skip = 4, header = FALSE)[, 1]

cts <- do.call(cbind, lapply(files, read_star))
rownames(cts) <- gene_ids
colnames(cts) <- sample_names

cts <- cts[rownames(cts) != "MissingGeneID", ]

metadata <- read.csv("/Users/jacobdiek/Coding/projects/differential-expression-blast/R/25030-01-final-QC-summary.csv")

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

# PCA-based outlier check — inspect the plot and sample labels, then
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
  
  write.csv(sig_genes, paste0("significant_genes_", comparison_name, ".csv"), row.names = FALSE)
}

print(sig_gene_counts)
