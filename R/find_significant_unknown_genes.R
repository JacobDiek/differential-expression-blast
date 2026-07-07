library("DESeq2")
library('apeglm')

# Make cts table

path <- "/Users/jacobdiek/Downloads/Coding/R/6.25.genetables.namedgenes"

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

metadata <- read.csv("/Users/jacobdiek/Downloads/Coding/R/25030-01-final-QC-summary.csv")

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
res <- results(dds)
summary(res)