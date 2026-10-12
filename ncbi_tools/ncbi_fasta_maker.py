from core.nucleotide_search import nucleotide_search


def ncbi_fasta_maker(input_file, output_file, **genome_settings):
    # genome_settings (genome_fasta, gff_path, gene_id_attribute) are passed to nucleotide_search
    with open(input_file) as significant_unknown_genes:
        significant_gene_list = significant_unknown_genes.read().splitlines()

    unknown_gene_dict = nucleotide_search(significant_gene_list, **genome_settings)

    x = 0

    with open(output_file, "w") as file:
        for gene in unknown_gene_dict:
            file.write(f">{gene}\n")
            file.write(f"{unknown_gene_dict[gene]}\n\n")
            x += 1

            if x % 100 == 0:
                print(f"{x} genes done!")
