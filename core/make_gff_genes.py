import pandas as pd


def parse_attributes(attr_string):
    d = {}
    for pair in attr_string.split(";"):
        if "=" in pair:
            key, value = pair.split("=", 1)
            d[key] = value
    return d


def make_gff_genes():
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
        "brook_genome/brook_genome_hap1.gff",
        sep="\t",
        comment="#",
        header=None,
        names=gff_columns,
    )

    genes = df[df["type"] == "gene"]

    parsed = genes["attributes"].apply(parse_attributes)
    attr_df = pd.DataFrame(list(parsed))

    result = pd.concat(
        [
            genes[["seqid", "start", "end", "strand"]].reset_index(drop=True),
            attr_df[["locus_tag"]].reset_index(drop=True),
        ],
        axis=1,
    )

    return result
