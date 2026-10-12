import os

import pandas as pd
from Bio import SeqIO

from core.blast_search import blast_search
from core.make_gff_genes import make_gff_genes
from core.nucleotide_search import nucleotide_search
from ncbi_tools.ncbi_gene_identifier import gene_identifier
from ncbi_tools.ncbi_mass_parser import ncbi_mass_parser
from pipeline.errors import ConfigError


def get_ncbi_email(email=None) -> str:
    # NCBI asks for an email with every online BLAST request. Pass --email, set it once with
    # `export NCBI_EMAIL=you@example.com`, or type it in when asked.
    email = (email or os.environ.get("NCBI_EMAIL", "")).strip()
    while not email:
        email = input("Enter your email for NCBI BLAST: ").strip()
    return email


def collect_sequences(settings, gene_ids, fasta_path):
    # {name: sequence} from gene IDs (looked up in the genome) and/or a FASTA file
    sequences = {}
    if gene_ids:
        for key in ("genome_fasta", "genome_gff"):
            if not settings.get(key) or not settings[key].is_file():
                raise ConfigError(f"Looking up gene IDs needs {key} in the config (not found: {settings.get(key)})")

        # Check the IDs first, so a typo fails before the slow genome loading
        known_ids = set(make_gff_genes(settings["genome_gff"], settings["gene_id_attribute"])["locus_tag"])
        unknown = [g for g in gene_ids if g not in known_ids]
        if unknown:
            raise ConfigError(f"Not found in the GFF ({settings['gene_id_attribute']}): {', '.join(unknown)}")

        print(f"Looking up {len(gene_ids)} gene sequences (loading the genome takes a minute)...")
        found = nucleotide_search(gene_ids, genome_fasta=settings["genome_fasta"], gff_path=settings["genome_gff"],
                                  gene_id_attribute=settings["gene_id_attribute"])
        sequences.update({gene: str(seq) for gene, seq in found.items()})

    if fasta_path:
        if not os.path.isfile(fasta_path):
            raise ConfigError(f"FASTA file not found: {fasta_path}")
        for record in SeqIO.parse(fasta_path, "fasta"):
            sequences[record.id] = str(record.seq)

    return sequences


def run_quick_blast(settings, gene_ids, fasta_path=None, email=None):
    sequences = collect_sequences(settings, gene_ids, fasta_path)
    if not sequences:
        raise ConfigError("Nothing to BLAST: give gene IDs and/or --fasta with sequences")

    out_dir = settings["output_dir"] / "quick_blast"
    for folder in ("raw_results", "parsed_results", "identified_genes"):
        (out_dir / folder).mkdir(parents=True, exist_ok=True)

    blast_search(email=get_ncbi_email(email), gene_list=sequences, file_location=out_dir)

    # Same parsing and naming as the bulk pipeline
    ncbi_mass_parser(out_dir / "raw_results", out_dir / "parsed_results")
    for old_file in (out_dir / "identified_genes").glob("*.txt"):
        old_file.unlink()
    summary = pd.DataFrame(gene_identifier(out_dir / "parsed_results", out_dir / "identified_genes"))
    summary = summary.set_index("gene_id").sort_index()
    summary.to_csv(out_dir / "summary.csv")

    this_run = summary.loc[[g for g in sequences if g in summary.index]]
    print()
    for gene, row in this_run.iterrows():
        if row["blast_status"] == "identified":
            print(f"{gene}: {row['blast_top_hit']} ({row['name_support']:.0%} support; "
                  f"best single hit: {row['best_single_hit']})")
        else:
            print(f"{gene}: unknown ({row['informative_hits']} informative, "
                  f"{row['uninformative_hits']} uninformative hits)")
    print(f"\nDetails in {out_dir / 'identified_genes'}, all quick-blast results in {out_dir / 'summary.csv'}")
    return this_run
