from pathlib import Path

from Bio import SeqIO

from core.make_gff_genes import make_gff_genes
from core.nucleotide_search import nucleotide_search
from pipeline.report import sheet_name

TEST_DATA = Path(__file__).resolve().parent / "data"


def test_sequences_come_from_the_right_place_and_strand():
    genes = make_gff_genes(TEST_DATA / "genome.gff", "locus_tag").set_index("locus_tag")
    genome = SeqIO.to_dict(SeqIO.parse(TEST_DATA / "genome.fasta", "fasta"))
    plus, minus = "testgene_0000", "testgene_0001"
    assert genes.at[plus, "strand"] == "+" and genes.at[minus, "strand"] == "-"

    sequences = nucleotide_search([plus, minus], genome_fasta=TEST_DATA / "genome.fasta",
                                  gff_path=TEST_DATA / "genome.gff")

    for gene, sequence in sequences.items():
        row = genes.loc[gene]
        forward = genome[row["seqid"]].seq[row["start"] - 1:row["end"]]
        expected = forward if row["strand"] == "+" else forward.reverse_complement()
        assert str(sequence) == str(expected)
        assert len(sequence) == row["end"] - row["start"] + 1


def test_excel_sheet_names_are_valid_and_unique():
    used = {"summary"}
    long_name = "liver_male_vs_liver_female_and_more_text"
    first, second = sheet_name(long_name, used), sheet_name(long_name, used)
    assert len(first) <= 31 and len(second) <= 31
    assert first != second
    assert sheet_name("a/b:c", used) == "a_b_c"
