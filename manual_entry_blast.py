from pathlib import Path

from core.blast_parse import blast_parse
from core.blast_search import blast_search
from core.gene_input import manual_gene_input


def manual_entry_blast():

    file_location = Path(
        "/Users/jacobdiek/Downloads/Coding/Python Projects/Python Blast/Results"
    )
    (file_location / "raw_results").mkdir(parents=True, exist_ok=True)
    (file_location / "parsed_results").mkdir(parents=True, exist_ok=True)

    gene_input_dictionary = manual_gene_input()

    blast_search(
        email="diek.jacob@gmail.com",
        gene_list=gene_input_dictionary,
        file_location=file_location,
    )

    blast_parse(gene_list=gene_input_dictionary, file_location=file_location)
