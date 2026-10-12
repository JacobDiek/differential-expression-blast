from pathlib import Path
from time import sleep

from Bio import Blast


def blast_search(email: str, gene_list: dict, file_location: Path, database="core_nt", megablast=True):
    # Same database and program as the NCBI website defaults used for the bulk BLAST step
    Blast.email = email
    raw_file_location = file_location / "raw_results"

    for gene in gene_list:
        # FASTA header so the results are labeled with the gene name
        query = f">{gene}\n{gene_list[gene]}\n"

        print(f"Submitting BLAST search for {gene} (usually takes a minute or two)...")
        # NCBI's URL API expects MEGABLAST=on
        result_stream = Blast.qblast("blastn", database, query, megablast="on" if megablast else None)
        print("Search complete, saving results...")

        file_name = f"blast_{gene}.xml"
        file_path = raw_file_location / file_name

        with open(file_path, "wb") as out_stream:
            out_stream.write(result_stream.read())
        result_stream.close()

        sleep(5)
