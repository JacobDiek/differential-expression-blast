"""quick-blast with NCBI replaced by a fake that returns prepared BLAST XML (no network)."""
import io
import subprocess
import sys
from pathlib import Path

import pytest

import core.blast_search
from core.nucleotide_search import nucleotide_search
from pipeline.errors import ConfigError
from pipeline.quick_blast import run_quick_blast
from pipeline.settings import read_settings
from tests.make_test_data import blast_xml, hits_majority, hits_noise

PROJECT_ROOT = Path(__file__).resolve().parent.parent
TEST_DATA = Path(__file__).resolve().parent / "data"
FAKE_RESULTS = {"testgene_0000": hits_majority(0)[0], "my_sequence": hits_noise(1)[0]}


@pytest.fixture
def fake_ncbi(monkeypatch):
    calls = []

    def fake_qblast(program, database, query, megablast=None):
        name = query.splitlines()[0].removeprefix(">")
        calls.append({"program": program, "database": database, "megablast": megablast,
                      "name": name, "sequence": "".join(query.splitlines()[1:])})
        return io.BytesIO(blast_xml([(name, FAKE_RESULTS[name])]).encode())

    monkeypatch.setattr(core.blast_search.Blast, "qblast", fake_qblast)
    monkeypatch.setattr(core.blast_search, "sleep", lambda seconds: None)
    return calls


@pytest.fixture
def settings(tmp_path):
    settings = read_settings(TEST_DATA / "test_config.yaml")
    settings["output_dir"] = tmp_path
    return settings


def test_quick_blast_gene_ids_and_fasta(fake_ncbi, settings, tmp_path):
    fasta = tmp_path / "mine.fasta"
    fasta.write_text(">my_sequence\nACGTACGTACGTACGT\n")

    results = run_quick_blast(settings, ["testgene_0000"], fasta_path=fasta, email="test@example.com")

    assert [c["name"] for c in fake_ncbi] == ["testgene_0000", "my_sequence"]
    assert all(c["program"] == "blastn" and c["database"] == "core_nt" and c["megablast"] == "on"
               for c in fake_ncbi)
    genome_sequence = nucleotide_search(["testgene_0000"], genome_fasta=TEST_DATA / "genome.fasta",
                                        gff_path=TEST_DATA / "genome.gff")["testgene_0000"]
    assert fake_ncbi[0]["sequence"] == str(genome_sequence)
    assert fake_ncbi[1]["sequence"] == "ACGTACGTACGTACGT"

    assert results.at["testgene_0000", "blast_status"] == "identified"
    assert results.at["testgene_0000", "blast_top_hit"] == "zona pellucida sperm binding protein 4 like"
    assert results.at["my_sequence", "blast_status"] == "unknown"
    assert (tmp_path / "quick_blast" / "summary.csv").exists()


def test_quick_blast_rejects_unknown_gene_ids_before_blasting(fake_ncbi, settings):
    with pytest.raises(ConfigError, match="Not found in the GFF \\(locus_tag\\): testgene_9999"):
        run_quick_blast(settings, ["testgene_0000", "testgene_9999"], email="test@example.com")
    assert fake_ncbi == []


def test_quick_blast_needs_genes_or_fasta():
    result = subprocess.run([sys.executable, "run_pipeline.py", "quick-blast"], cwd=PROJECT_ROOT,
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "give gene IDs and/or --fasta" in result.stderr
