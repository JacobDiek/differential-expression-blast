from Bio import SeqIO

input_file = "ncbi_unknown.fasta"
output_file = "ncbi_unknown_filter.fasta"
max_bases_per_gene = 20_000

genes = {}

with open(input_file, "r") as file:
    gene_name = None
    sequence = []

    for line in file:
        line = line.strip()

        if not line:
            continue  # Skip blank lines

        if line.startswith(">"):
            # Save the previous gene
            if gene_name is not None:
                genes[gene_name] = "".join(sequence)

            # Start a new gene
            gene_name = line[1:]  # Remove the ">"
            sequence = []
        else:
            sequence.append(line)

    # Save the last gene
    if gene_name is not None:
        genes[gene_name] = "".join(sequence)

filtered_gene_dict = {k:genes[k] for k in genes if len(genes[k]) <= max_bases_per_gene}

print(len(genes))
print(len(filtered_gene_dict))

with open(output_file, "w") as file:
    for gene in filtered_gene_dict:
        file.write(f">{gene}\n")
        file.write(f"{filtered_gene_dict[gene]}\n\n")
