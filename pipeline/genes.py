import re
from urllib.parse import unquote

import pandas as pd

from core.make_gff_genes import read_gff_genes
from pipeline.deseq2 import significant_genes_path
from pipeline.errors import ConfigError


def load_gene_table(settings):
    # One row per GFF gene: ID, name, description, biotype, length, and whether it is "unknown"
    id_attribute = settings["gene_id_attribute"]
    genes = read_gff_genes(settings["genome_gff"])

    if id_attribute not in genes.columns:
        raise ConfigError(f"gene_id_attribute '{id_attribute}' is not an attribute of the GFF genes. "
                          f"Some attributes it does have: {[c for c in genes.columns if c not in ('seqid', 'start', 'end', 'strand')][:10]}")

    table = pd.DataFrame({
        "gene_id": genes[id_attribute],
        "gene_name": genes.get("Name"),
        "gff_description": genes.get("description"),
        "gene_biotype": genes.get("gene_biotype"),
        "length": genes["end"] - genes["start"] + 1,
    }).dropna(subset=["gene_id"]).drop_duplicates(subset="gene_id")

    # GFF values are URL-encoded (e.g. "%2C" for ",")
    table["gff_description"] = table["gff_description"].map(lambda d: unquote(d) if isinstance(d, str) else d)

    # Unknown: no real name, i.e. Name is missing, just the gene ID, or LOC-style
    patterns = [re.compile(p) for p in settings["unknown_name_patterns"]]
    no_real_name = (table["gene_name"].isna()
                    | (table["gene_name"] == table["gene_id"])
                    | table["gene_name"].map(lambda n: isinstance(n, str) and any(p.fullmatch(n) for p in patterns)))
    table["unknown_in_gff"] = no_real_name

    return table.set_index("gene_id")


def load_significant_genes(settings, deseq2_dir):
    # {comparison: DataFrame of significant genes}
    return {tuple(c): pd.read_csv(significant_genes_path(deseq2_dir, c)) for c in settings["contrasts"]}


def check_gene_ids_match(gene_table, significant):
    sig_ids = set().union(*(set(df["gene_id"]) for df in significant.values()))
    if sig_ids and not sig_ids & set(gene_table.index):
        raise ConfigError("None of the significant gene IDs from the count tables are in the GFF. "
                          "Check that gene_id_attribute matches the IDs in the count tables.")


def select_genes_to_blast(settings, gene_table, significant):
    # Significant in any comparison + unknown in the GFF (+ allowed biotype)
    sig_ids = set().union(*(set(df["gene_id"]) for df in significant.values()))
    candidates = gene_table.loc[gene_table.index.isin(sig_ids) & gene_table["unknown_in_gff"]]

    if settings["blast_biotypes"]:
        candidates = candidates[candidates["gene_biotype"].isin(settings["blast_biotypes"])]

    too_long = candidates["length"] > settings["max_gene_length"]
    return sorted(candidates.index[~too_long]), sorted(candidates.index[too_long])
