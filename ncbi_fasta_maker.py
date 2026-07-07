from core.nucleotide_search import nucleotide_search

with open("significant_genes.txt") as significant_unknown_genes:
    significant_gene_list = significant_unknown_genes.read().splitlines()

unknown_gene_dict = nucleotide_search(significant_gene_list)

x = 0

with open("ncbi_unknown.fasta", "w") as file:
    for gene in unknown_gene_dict:
        file.write(f">{gene}\n")
        file.write(f"{unknown_gene_dict[gene]}\n\n")
        x += 1

        if x % 100 == 0:
            print(f"{x} genes done!")
