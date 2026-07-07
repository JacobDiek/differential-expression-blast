from pathlib import Path
from time import sleep

from Bio import Blast
from Bio.Seq import Seq


def blast_search(email: str, gene_list: dict, file_location: Path):
    Blast.email = email
    raw_file_location = file_location / "raw_results"

    for gene in gene_list:
        nucleotide_string = Seq(gene_list[gene])

        print("Submitting BLAST search...")
        result_stream = Blast.qblast("blastn", "nt", nucleotide_string)
        print("Search complete, saving results...")

        file_name = f"blast_{gene}.xml"
        file_path = raw_file_location / file_name

        with open(file_path, "wb") as out_stream:
            out_stream.write(result_stream.read())
        result_stream.close()

        sleep(5)
