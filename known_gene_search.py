import pandas as pd


def parse_attributes(attr_string):
    d = {}
    for pair in attr_string.split(";"):
        if "=" in pair:
            key, value = pair.split("=", 1)
            d[key] = value
    return d


def load_gene_descriptions(gff_path):
    gff_columns = [
        "seqid", "source", "type", "start", "end",
        "score", "strand", "phase", "attributes",
    ]
    df = pd.read_csv(gff_path, sep="\t", comment="#", header=None, names=gff_columns)
    genes = df[df["type"] == "gene"]

    parsed = genes["attributes"].apply(parse_attributes)
    attr_df = pd.DataFrame(list(parsed))

    return attr_df[["locus_tag", "description"]].dropna(subset=["locus_tag"])


def annotate_top_genes(sig_csv_path, gff_path, top_n=30, rank_by="padj", output_path=None):
    sig = pd.read_csv(sig_csv_path)
    descriptions = load_gene_descriptions(gff_path)

    if rank_by == "abs_log2fc":
        sig = sig.reindex(sig["log2FoldChange"].abs().sort_values(ascending=False).index)
    else:
        sig = sig.sort_values("padj")

    top = sig.head(top_n)
    annotated = top.merge(descriptions, left_on="gene_id", right_on="locus_tag", how="left")
    annotated = annotated[["gene_id", "log2FoldChange", "padj", "description"]]

    if output_path:
        annotated.to_csv(output_path, index=False)

    return annotated


# Liver male vs. liver female — rank by padj (values aren't saturated at 0 here)
liver_top = annotate_top_genes(
    sig_csv_path="/Users/jacobdiek/Coding/projects/differential-expression-blast/R/significant_genes_liver_male_vs_liver_female.csv",
    gff_path="/Users/jacobdiek/Coding/projects/differential-expression-blast/brook_genome/brook_genome_hap1.gff",
    top_n=100,
    rank_by="padj",
    output_path="annotated_top_liver_male_vs_liver_female.csv",
)
print(liver_top)

# Testes vs. ovary — rank by |log2FC| since many padj values are saturated at 0
gonad_top = annotate_top_genes(
    sig_csv_path="/Users/jacobdiek/Coding/projects/differential-expression-blast/R/significant_genes_testes_vs_ovary.csv",
    gff_path="/Users/jacobdiek/Coding/projects/differential-expression-blast/brook_genome/brook_genome_hap1.gff",
    top_n=100,
    rank_by="abs_log2fc",
    output_path="annotated_top_testes_vs_ovary.csv",
)
print(gonad_top)