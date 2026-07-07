from Bio import SeqIO

from core.make_gff_genes import make_gff_genes


def nucleotide_search(gene_list: list) -> dict:

    gff_table = make_gff_genes()

    sequence_dict = {}
    scaffold_dict = SeqIO.to_dict(
        SeqIO.parse("brook_genome/brook_genome_hap1_v1.fasta", "fasta")
    )

    lookup = gff_table.set_index("locus_tag")
    for locus_tag in gene_list:
        series = lookup.loc[locus_tag]

        seqid = series["seqid"]
        start = series["start"]
        end = series["end"]
        strand = series["strand"]

        record = scaffold_dict[seqid]
        nucleotide_sequence = record.seq[start - 1 : end]

        if strand == "-":
            nucleotide_sequence = nucleotide_sequence.reverse_complement()

        sequence_dict[locus_tag] = nucleotide_sequence

    return sequence_dict
