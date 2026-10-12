import os

import pandas as pd

from ncbi_tools.fasta_file_splitter import fasta_file_splitter
from ncbi_tools.ncbi_fasta_maker import ncbi_fasta_maker
from ncbi_tools.ncbi_gene_identifier import extract_locus_tag, gene_identifier
from ncbi_tools.ncbi_mass_parser import ncbi_mass_parser
from ncbi_tools.ncbi_xml_splitter import split_blast_xml

INSTRUCTIONS = """\
BLAST batch {batch_name}: {n_genes} genes in {n_chunks} chunk files

1. Open NCBI Nucleotide BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn
2. For each chunk file in:
     {chunks_dir}
   a. Click "Choose File" under "Enter Query Sequence" and pick the chunk file.
   b. Database: "Core nucleotide database (core_nt)". Program: "Highly similar
      sequences (megablast)". Leave the other settings at their defaults.
   c. Click BLAST and wait for the results page.
   d. Click "Download All" and choose "XML" (not "Single-file XML2").
   e. Move the downloaded "...-Alignment.xml" file into:
        {xml_dir}
3. When every chunk is downloaded, run the pipeline again: python run_pipeline.py
"""


class BlastFolders:
    def __init__(self, ncbi_dir):
        self.ncbi_dir = ncbi_dir
        self.xml_dir = ncbi_dir / "ncbi_files_large_alignment"   # downloaded BLAST XMLs go here
        self.split_dir = ncbi_dir / "ncbi_files_raw_output"
        self.parsed_dir = ncbi_dir / "ncbi_files_parsed_output"
        self.identified_dir = ncbi_dir / "identified_genes"
        self.batches_dir = ncbi_dir / "blast_batches"

    def make(self):
        for folder in (self.xml_dir, self.split_dir, self.parsed_dir, self.identified_dir, self.batches_dir):
            folder.mkdir(parents=True, exist_ok=True)


def process_blast_results(folders):
    # Split any newly downloaded XMLs per gene, then parse and name every gene with results
    split_blast_xml(folders.xml_dir, folders.split_dir, skip_already_split=True)
    ncbi_mass_parser(folders.split_dir, folders.parsed_dir)

    # identified_genes/ is rebuilt from scratch every run
    for old_file in folders.identified_dir.glob("*.txt"):
        old_file.unlink()
    summary = gene_identifier(folders.parsed_dir, folders.identified_dir)

    summary = pd.DataFrame(summary, columns=["gene_id", "blast_status", "blast_top_hit", "name_support",
                                             "best_single_hit", "informative_hits", "uninformative_hits"])
    summary.sort_values("gene_id").to_csv(folders.ncbi_dir / "blast_summary.csv", index=False)
    return summary.set_index("gene_id")


def genes_with_results(folders):
    return {extract_locus_tag(f) for f in os.listdir(folders.parsed_dir)}


def read_batches(folders):
    # {batch folder: list of genes in that batch}
    batches = {}
    for batch_dir in sorted(folders.batches_dir.glob("batch_*"), key=lambda p: int(p.name.split("_")[1])):
        genes_file = batch_dir / "genes.txt"
        if genes_file.exists():
            batches[batch_dir] = genes_file.read_text().split()
    return batches


def make_batch(settings, folders, genes):
    existing = [int(p.name.split("_")[1]) for p in folders.batches_dir.glob("batch_*")]
    batch_dir = folders.batches_dir / f"batch_{max(existing, default=0) + 1}"
    chunks_dir = batch_dir / "chunks"
    batch_dir.mkdir(parents=True)

    genes_file = batch_dir / "genes.txt"
    genes_file.write_text("\n".join(genes) + "\n")

    print(f"Writing sequences for {len(genes)} genes to {batch_dir.name}/ (loading the genome takes a minute)...")
    ncbi_fasta_maker(genes_file, batch_dir / "genes.fasta",
                     genome_fasta=settings["genome_fasta"],
                     gff_path=settings["genome_gff"],
                     gene_id_attribute=settings["gene_id_attribute"])
    fasta_file_splitter(batch_dir / "genes.fasta", chunks_dir,
                        max_bases_per_chunk=settings["chunk_max_bases"],
                        max_seqs_per_chunk=settings["chunk_max_seqs"])

    write_instructions(batch_dir, folders, genes)
    return batch_dir


def write_instructions(batch_dir, folders, genes):
    instructions = INSTRUCTIONS.format(
        batch_name=batch_dir.name,
        n_genes=len(genes),
        n_chunks=len(list((batch_dir / "chunks").glob("*.fasta"))),
        chunks_dir=batch_dir / "chunks",
        xml_dir=folders.xml_dir,
    )
    (batch_dir / "INSTRUCTIONS.txt").write_text(instructions)


def batch_status(folders, batches, have_results):
    # Pending = none of its genes have results yet (still waiting to be BLASTed)
    # Incomplete = some genes have results but others are missing
    pending, incomplete = [], {}
    for batch_dir, genes in batches.items():
        missing = [g for g in genes if g not in have_results]
        if len(missing) == len(genes):
            pending.append(batch_dir)
        elif missing:
            incomplete[batch_dir] = missing
            (batch_dir / "missing_results.txt").write_text("\n".join(missing) + "\n")
        else:
            (batch_dir / "missing_results.txt").unlink(missing_ok=True)
    return pending, incomplete
