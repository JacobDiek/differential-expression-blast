import pandas as pd


def parse_attributes(attr_string):
    d = {}
    for pair in attr_string.split(";"):
        if "=" in pair:
            key, value = pair.split("=", 1)
            d[key] = value
    return d


def read_gff_genes(gff_path):
    # Every "gene" row of the GFF, with its attributes split into columns
    gff_columns = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]

    df = pd.read_csv(
        gff_path,
        sep="\t",
        comment="#",
        header=None,
        names=gff_columns,
    )

    genes = df[df["type"] == "gene"]

    parsed = genes["attributes"].apply(parse_attributes)
    attr_df = pd.DataFrame(list(parsed))

    return pd.concat(
        [
            genes[["seqid", "start", "end", "strand"]].reset_index(drop=True),
            attr_df.reset_index(drop=True),
        ],
        axis=1,
    )


def make_gff_genes(gff_path, gene_id_attribute="locus_tag"):
    # gene_id_attribute is the GFF attribute that holds the gene IDs used in the
    # count tables (locus_tag for the brook stickleback annotation)
    genes = read_gff_genes(gff_path)

    result = genes[["seqid", "start", "end", "strand", gene_id_attribute]]
    result = result.rename(columns={gene_id_attribute: "locus_tag"})

    return result
