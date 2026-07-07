import json

from core.nucleotide_search import nucleotide_search


def manual_gene_input() -> dict:
    new_list = (
        input("Do you want to start a new list of genes? [Y/N]: ").strip().upper()
    )
    if new_list == "Y":
        gene_input_dictionary = {}

    else:
        try:
            with open("genes_input.json") as json_file:
                gene_input_dictionary = json.load(json_file)

        except json.decoder.JSONDecodeError:
            print("You do not have a current list of genes")
            gene_input_dictionary = {}

        except FileNotFoundError:
            print("No existing gene list found, starting a new one.")
            gene_input_dictionary = {}

    add_genes = input("Do you want to add a gene to your list [Y/N]: ").strip().upper()
    if add_genes == "Y":
        while True:
            gene_name = input("Enter gene name: ").strip()
            gene_sequence = input("Enter gene sequence: ").strip().upper()

            gene_input_dictionary[gene_name] = gene_sequence

            repeat_gene_entry = input(
                "Do you have another gene to enter [Y/N]: "
            ).upper()
            if repeat_gene_entry != "Y":
                break

        save_list = (
            input("Do you want to save your list of genes? [Y/N]: ").strip().upper()
        )
        if save_list != "N":
            with open("genes_input.json", "w") as f:
                json.dump(gene_input_dictionary, f, indent=4)

    return gene_input_dictionary


def automatic_gene_input() -> dict:
    locus_tag_input = []

    while True:
        locus_tag = input("Enter locus tag: ").strip().lower()
        locus_tag_input.append(locus_tag)

        repeat_gene_entry = input("Do you have another gene to enter [Y/N]: ").upper()
        if repeat_gene_entry != "Y":
            break

    return nucleotide_search(locus_tag_input)
