import pytest

from ncbi_tools.ncbi_gene_identifier import (MAX_HIT_WEIGHT, Hit, find_multi_word_species, gene_identifier,
                                             hit_weight, vote_key, weighted_vote)


@pytest.mark.parametrize("hit_def, species, name", [
    ("PREDICTED: Gasterosteus aculeatus solute carrier organic anion transporter family, member 1E1 (slco1e1), "
     "transcript variant X2, mRNA",
     "Gasterosteus aculeatus", "solute carrier organic anion transporter family, member 1e1"),
    ("TPA_asm: Oryzias latipes cab protein (cab), mRNA", "Oryzias latipes", "cab protein"),
    ("PREDICTED: Rhamphochromis sp. 'chilingali' tetraspanin-4-like (LOC1), mRNA",
     "Rhamphochromis sp. 'chilingali'", "tetraspanin 4 like"),
    ("PREDICTED: Xiphophorus maculatus x Xiphophorus helleri claudin 4 (cldn4), mRNA",
     "Xiphophorus maculatus x Xiphophorus helleri", "claudin 4"),
    ("[Epinephelus coioides] Epinephelus coioides EcSCYA113 mRNA for CC chemokine, complete cds",
     "Epinephelus coioides", "ecscya113 mrna for cc chemokine, complete cds"),
])
def test_species_is_removed_from_hit_name(hit_def, species, name):
    hit = Hit(hit_def)
    assert hit.species == species
    assert hit.normalized_name == name


def test_three_part_species_names_are_learned_from_the_hits():
    hit_defs = [
        "PREDICTED: Seriola lalandi dorsalis claudin-4 (LOC1), mRNA",
        "PREDICTED: Seriola lalandi dorsalis keratin 5 (LOC2), mRNA",
        "PREDICTED: Gasterosteus aculeatus claudin-4 (cldn4), mRNA",
        "PREDICTED: Danio rerio keratin 5 (krt5), mRNA",
        # "neural" starts real gene names, so it must not become part of an organism name
        "PREDICTED: Gasterosteus aculeatus neural cell adhesion molecule 1 (ncam1), mRNA",
        "PREDICTED: Danio rerio cell adhesion molecule 1 (cadm1), mRNA",
        "PREDICTED: Danio rerio neural cell adhesion molecule 2 (ncam2), mRNA",
    ]
    multi_word_species = find_multi_word_species(hit_defs)
    assert multi_word_species == {"Seriola lalandi dorsalis"}

    hit = Hit(hit_defs[0], multi_word_species)
    assert hit.species == "Seriola lalandi dorsalis"
    assert hit.normalized_name == "claudin 4"
    assert Hit(hit_defs[4], multi_word_species).normalized_name == "neural cell adhesion molecule 1"


@pytest.mark.parametrize("evalue, weight", [
    ("0", MAX_HIT_WEIGHT), ("0.0", MAX_HIT_WEIGHT), ("1e-50", 50), ("1e-300", MAX_HIT_WEIGHT),
    ("0.5", 1), ("5", 1), ("not a number", 1), (None, 1),
])
def test_hit_weight(evalue, weight):
    assert hit_weight(evalue) == pytest.approx(weight)


def test_vote_key_groups_versions_of_the_same_name():
    same = ["zona pellucida sperm binding protein 4", "zona pellucida sperm binding protein 4 like",
            "zona pellucida sperm binding protein 4 like transcript variant x3",
            "zona pellucida sperm binding protein 4, mrna", "zona pellucida sperm binding protein 4 homolog"]
    assert len({vote_key(name) for name in same}) == 1
    assert vote_key("claudin 4") != vote_key("claudin 5")


def test_weighted_vote_prefers_the_majority_over_the_single_best_hit():
    hits = [("titin like", 180)] + [("zona pellucida sperm binding protein 4 like", 140)] * 4
    name, support, _ = weighted_vote(hits)
    assert name == "zona pellucida sperm binding protein 4 like"
    assert support == pytest.approx(560 / 740)


def test_weighted_vote_ties_go_to_the_better_ranked_hit():
    name, support, _ = weighted_vote([("gene a", 50), ("gene b", 50)])
    assert name == "gene a"
    assert support == pytest.approx(0.5)


def write_parsed(path, gene, hits):
    # Same format as ncbi_mass_parser output
    text = f"Query: {gene}\n\n"
    for number, (hit_def, evalue) in enumerate(hits, start=1):
        text += f"Hit Num: {number}\ngi|{number}|ref|XM_1|\n{hit_def}\nEvalue: {evalue}\n\n"
    (path / f"blast_output_{gene}.txt").write_text(text)


def test_gene_identifier(tmp_path):
    parsed, identified = tmp_path / "parsed", tmp_path / "identified"
    parsed.mkdir()
    lectin = [("PREDICTED: Danio rerio lectin 1 (lec1), mRNA", "1e-80")] * 3
    write_parsed(parsed, "gene_a", lectin)
    write_parsed(parsed, "gene_b", lectin)   # same name as gene_a: must not overwrite its file
    write_parsed(parsed, "gene_noise", [("Danio rerio DNA, chromosome 3, sequence", "0"),
                                        ("Danio rerio clone XYZ, complete sequence", "0")])
    write_parsed(parsed, "gene_uncharacterized", [("PREDICTED: Danio rerio uncharacterized LOC1 (LOC1), ncRNA", "0")] * 3
                 + [("PREDICTED: Danio rerio keratin 5 (krt5), mRNA", "1e-10")])
    write_parsed(parsed, "gene_no_hits", [])

    summary = {row["gene_id"]: row for row in gene_identifier(parsed, identified)}

    assert summary["gene_a"]["blast_status"] == "identified"
    assert summary["gene_a"]["blast_top_hit"] == "lectin 1"
    assert summary["gene_a"]["name_support"] == 1.0
    for gene in ("gene_noise", "gene_uncharacterized", "gene_no_hits"):
        assert summary[gene]["blast_status"] == "unknown"
        assert summary[gene]["blast_top_hit"] is None

    files = sorted(f.name for f in identified.iterdir())
    assert "lectin_1__gene_a.txt" in files and "lectin_1__gene_b.txt" in files
    assert "UNKNOWN_gene_no_hits.txt" in files
    assert len(files) == 5
