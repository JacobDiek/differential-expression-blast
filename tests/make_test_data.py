"""Builds the small synthetic test dataset in tests/data/ (already committed; rerun only to change it).

    python tests/make_test_data.py

Everything is made up: a random genome, a GFF with named and unknown genes, STAR count
files for 7 samples with known differentially expressed genes, and a BLAST XML "download"
whose hits exercise the gene-naming rules. expected.json records what the pipeline should find.
"""
import json
import random
import shutil
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent / "data"
SEED = 7

N_GENES = 120
N_DE = 30                 # genes 0-29 are differentially expressed
LONG_GENE = 4             # unknown + DE, but longer than max_gene_length in the test config
MISSING_FROM_BLAST = 25   # unknown + DE, left out of the BLAST "download"
SAMPLES = [("sampleC1", "control"), ("sampleC2", "control"), ("sampleC3", "control"),
           ("sampleT1", "treated"), ("sampleT2", "treated"), ("sampleT3", "treated"),
           ("sampleT4", "treated")]   # sampleT4 is an outlier, excluded in the test config
OUTLIER = "sampleT4"


def gene_id(i):
    return f"testgene_{i:04d}"


def gene_kind(i):
    # Two out of four genes are unknown: Name = gene ID, or a LOC-style name
    return {0: "id_name", 1: "loc_name"}.get(i % 4, "named")


def make_genome(rng):
    genes, scaffolds = [], {}
    for s in range(3):
        sequence, position = [], 0
        for i in range(s * N_GENES // 3, (s + 1) * N_GENES // 3):
            gap = "".join(rng.choice("ACGT") for _ in range(300))
            length = 6000 if i == LONG_GENE else rng.randint(400, 900)
            body = "".join(rng.choice("ACGT") for _ in range(length))
            sequence += [gap, body]
            start = position + len(gap) + 1
            genes.append({"i": i, "seqid": f"scaffold_{s + 1}", "start": start, "end": start + length - 1,
                          "strand": "-" if i % 2 else "+"})
            position += len(gap) + length
        scaffolds[f"scaffold_{s + 1}"] = "".join(sequence)
    return scaffolds, genes


def gff_line(gene):
    i = gene["i"]
    attributes = [f"ID=gene-{gene_id(i)}"]
    kind = gene_kind(i)
    if kind == "id_name":
        attributes += [f"Name={gene_id(i)}", "gene_biotype=" + ("protein_coding" if i % 8 == 0 else "lncRNA")]
    elif kind == "loc_name":
        attributes += [f"Name=LOC{100000 + i}", f"description=uncharacterized protein LOC{100000 + i}",
                       "gene_biotype=protein_coding"]
    else:
        attributes += [f"Name=tg{i}", f"description=test protein {i}%2C made up", "gene_biotype=protein_coding"]
    attributes.append(f"locus_tag={gene_id(i)}")
    return "\t".join([gene["seqid"], "Test", "gene", str(gene["start"]), str(gene["end"]), ".",
                      gene["strand"], ".", ";".join(attributes)])


def make_counts(np_rng):
    base = np_rng.uniform(200, 2000, N_GENES)
    fold = np.ones(N_GENES)
    fold[:N_DE] = np.where(np.arange(N_DE) % 2 == 0, 8.0, 1 / 8)   # up / down in "treated"
    dispersion = 0.05
    counts = {}
    for sample, condition in SAMPLES:
        mean = base * (fold if condition == "treated" else 1)
        if sample == OUTLIER:
            mean = mean * np_rng.uniform(0.05, 20, N_GENES)
        n = 1 / dispersion
        counts[sample] = np_rng.negative_binomial(n, n / (n + mean))
    return counts


def write_star_file(path, counts, np_rng):
    # Reverse-stranded library: column 4 has the reads, column 3 only ~2% wrong-strand reads
    wrong_strand = np_rng.binomial(counts, 0.02)
    lines = ["N_unmapped\t5000\t5000\t5000", "N_multimapping\t800\t800\t800",
             f"N_noFeature\t300\t{int(counts.sum())}\t{int(wrong_strand.sum())}", "N_ambiguous\t50\t20\t30"]
    lines += [f"{gene_id(i)}\t{c + w}\t{w}\t{c}" for i, (c, w) in enumerate(zip(counts, wrong_strand))]
    path.write_text("\n".join(lines) + "\n")


# BLAST hit templates: (hit definitions with e-values, expected status, expected name)
SPECIES = ["Gasterosteus aculeatus", "Pungitius pungitius", "Larimichthys crocea", "Oryzias latipes",
           "Danio rerio", "Esox lucius"]


def hits_majority(n):
    # Best single hit says titin, but most of the weight says zona pellucida
    hits = [(f"PREDICTED: Pungitius pungitius titin-like (LOC{n}1), mRNA", "0")]
    hits += [(f"PREDICTED: {sp} zona pellucida sperm-binding protein 4-like (LOC{n}{k}), mRNA", f"1e-{150 - k}")
             for k, sp in enumerate(SPECIES[1:], start=2)]
    return hits, "identified", "zona pellucida sperm binding protein 4 like"


def hits_clean(n):
    hits = [(f"PREDICTED: {sp} solute carrier family 1 member 2 (slc1a2), transcript variant X{k}, mRNA", f"1e-{120 - k}")
            for k, sp in enumerate(SPECIES)]
    return hits, "identified", "solute carrier family 1 member 2"


def hits_noise(n):
    hits = [(f"{sp} DNA, chromosome {k}, sequence", "1e-90") for k, sp in enumerate(SPECIES[:3])]
    hits += [(f"Danio rerio clone CH211-{n}, complete sequence", "1e-80"),
             (f"Esox lucius voucher KU{n} ultra-conserved element locus uce-308 genomic sequence", "1e-70")]
    return hits, "unknown", None


def hits_uninformative(n):
    hits = [(f"PREDICTED: {sp} uncharacterized LOC{n}{k} (LOC{n}{k}), ncRNA", "1e-60") for k, sp in enumerate(SPECIES)]
    hits.append((f"PREDICTED: Danio rerio keratin 5 (krt5), mRNA", "1e-20"))
    return hits, "unknown", None


def hits_none(n):
    return [], "unknown", None


def hits_subspecies(n):
    # "Seriola lalandi dorsalis" is a three-part organism name; it must not end up in the gene name
    hits = [(f"PREDICTED: Seriola lalandi dorsalis claudin-4 (LOC{n}{k}), mRNA", "0") for k in range(3)]
    hits += [(f"PREDICTED: {sp} claudin-4 (cldn4), mRNA", "1e-100") for sp in SPECIES[:3]]
    return hits, "identified", "claudin 4"


TEMPLATES = [hits_majority, hits_clean, hits_noise, hits_uninformative, hits_none, hits_subspecies]


def blast_xml(iterations):
    parts = ['<?xml version="1.0"?>',
             '<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">',
             "<BlastOutput>", "  <BlastOutput_program>blastn</BlastOutput_program>",
             "  <BlastOutput_version>BLASTN 2.17.0+</BlastOutput_version>",
             "  <BlastOutput_db>core_nt</BlastOutput_db>", "<BlastOutput_iterations>"]
    for number, (query, hits) in enumerate(iterations, start=1):
        parts += ["<Iteration>", f"  <Iteration_iter-num>{number}</Iteration_iter-num>",
                  f"  <Iteration_query-def>{query}</Iteration_query-def>", "<Iteration_hits>"]
        for hit_number, (hit_def, evalue) in enumerate(hits, start=1):
            hit_def = hit_def.replace("&", "&amp;").replace("<", "&lt;")
            parts += ["<Hit>", f"  <Hit_num>{hit_number}</Hit_num>", f"  <Hit_id>gi|{hit_number}|ref|XM_TEST.1|</Hit_id>",
                      f"  <Hit_def>{hit_def}</Hit_def>",
                      f"  <Hit_hsps><Hsp><Hsp_evalue>{evalue}</Hsp_evalue></Hsp></Hit_hsps>", "</Hit>"]
        parts += ["</Iteration_hits>", "</Iteration>"]
    parts += ["</BlastOutput_iterations>", "</BlastOutput>"]
    return "\n".join(parts) + "\n"


TEST_CONFIG = f"""\
# Config for the synthetic test dataset (made by tests/make_test_data.py).
# Try the pipeline on it with: python run_pipeline.py tests/data/test_config.yaml
genome_fasta: genome.fasta
genome_gff: genome.gff
counts_dir: counts
sample_sheet: samples.csv
output_dir: output

gene_id_attribute: locus_tag
star_count_column: auto
exclude_samples:
  - {OUTLIER}
contrasts:
  - [treated, control]
padj_cutoff: 0.05
log2fc_cutoff: 1

unknown_name_patterns:
  - "^LOC\\\\d+$"
blast_biotypes: []
max_gene_length: 5000
chunk_max_bases: 1000000
chunk_max_seqs: 5
"""


def main():
    rng = random.Random(SEED)
    np_rng = np.random.default_rng(SEED)

    if DATA.exists():
        shutil.rmtree(DATA)
    (DATA / "counts").mkdir(parents=True)
    (DATA / "blast_results").mkdir()

    scaffolds, genes = make_genome(rng)
    with open(DATA / "genome.fasta", "w") as f:
        for name, sequence in scaffolds.items():
            f.write(f">{name}\n")
            f.writelines(sequence[i:i + 80] + "\n" for i in range(0, len(sequence), 80))
    (DATA / "genome.gff").write_text("##gff-version 3\n" + "\n".join(gff_line(g) for g in genes) + "\n")

    counts = make_counts(np_rng)
    for lane, (sample, _) in enumerate(SAMPLES, start=1):
        write_star_file(DATA / "counts" / f"{sample}_S{lane}_L001_ReadsPerGene.out.tab", counts[sample], np_rng)
    (DATA / "samples.csv").write_text("sample,condition\n" + "".join(f"{s},{c}\n" for s, c in SAMPLES))

    unknown_de = [i for i in range(N_DE) if gene_kind(i) != "named"]
    to_blast = [i for i in unknown_de if i != LONG_GENE]
    iterations, expected_genes = [], {}
    for n, i in enumerate(to_blast):
        hits, status, name = TEMPLATES[n % len(TEMPLATES)](i)
        if i == MISSING_FROM_BLAST:
            expected_genes[gene_id(i)] = {"status": "waiting for BLAST results", "name": None}
            continue
        iterations.append((gene_id(i), hits))
        expected_genes[gene_id(i)] = {"status": status, "name": name}
    (DATA / "blast_results" / "TEST0001-Alignment.xml").write_text(blast_xml(iterations))

    expected = {
        "significant_genes": [gene_id(i) for i in range(N_DE)],
        "genes_to_blast": [gene_id(i) for i in to_blast],
        "too_long": [gene_id(LONG_GENE)],
        "missing_from_blast": [gene_id(MISSING_FROM_BLAST)],
        "blast": expected_genes,
    }
    (DATA / "expected.json").write_text(json.dumps(expected, indent=2) + "\n")
    (DATA / "test_config.yaml").write_text(TEST_CONFIG)
    print(f"Wrote test data to {DATA}")


if __name__ == "__main__":
    main()
