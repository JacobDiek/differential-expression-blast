import shutil
from pathlib import Path

import pytest
import yaml

from pipeline.errors import ConfigError
from pipeline.settings import load_settings
from pipeline.strandedness import detect_star_count_column

TEST_DATA = Path(__file__).resolve().parent / "data"


def write_star(path, forward, reverse):
    lines = ["N_unmapped\t0\t0\t0", "N_multimapping\t0\t0\t0", "N_noFeature\t0\t0\t0", "N_ambiguous\t0\t0\t0"]
    lines += [f"gene_{i}\t{f + r}\t{f}\t{r}" for i, (f, r) in enumerate(zip(forward, reverse))]
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.mark.parametrize("forward, reverse, column", [
    ([5, 3, 2], [500, 300, 200], 4),     # reverse stranded
    ([500, 300, 200], [5, 3, 2], 3),     # forward stranded
    ([100, 90, 110], [95, 105, 100], 2), # unstranded
    ([100, 100, 100], [200, 200, 200], None),  # 2x: unclear
])
def test_detect_star_count_column(tmp_path, forward, reverse, column):
    star_file = write_star(tmp_path / "s_ReadsPerGene.out.tab", forward, reverse)
    assert detect_star_count_column([star_file])[0] == column


@pytest.fixture
def test_config(tmp_path):
    # A copy of the test dataset's config that can be edited per test
    config = yaml.safe_load((TEST_DATA / "test_config.yaml").read_text())
    for key in ("genome_fasta", "genome_gff", "counts_dir", "sample_sheet"):
        config[key] = str(TEST_DATA / config[key])
    config["output_dir"] = str(tmp_path / "output")

    def write(**changes):
        path = tmp_path / "config.yaml"
        path.write_text(yaml.safe_dump({**config, **changes}))
        return path
    return write


def test_valid_config_loads_and_detects_strandedness(test_config, capsys):
    settings = load_settings(test_config())
    assert settings["star_count_column"] == 4
    assert "reverse stranded" in capsys.readouterr().out


def test_wrong_strandedness_setting_warns(test_config, capsys):
    settings = load_settings(test_config(star_count_column=3))
    assert settings["star_count_column"] == 3
    assert "WARNING" in capsys.readouterr().out


def test_config_problems_are_all_reported_together(test_config):
    with pytest.raises(ConfigError) as error:
        load_settings(test_config(genome_fasta="missing.fasta", star_count_column=7,
                                  exclude_samples=["nope"], contrasts=[["treated", "ctrl"]]))
    message = str(error.value)
    assert "genome_fasta: file not found" in message
    assert "star_count_column must be" in message
    assert "exclude_samples not in sample_sheet: ['nope']" in message
    assert "no samples with condition 'ctrl'" in message


def test_missing_count_file_is_reported(test_config, tmp_path):
    counts = tmp_path / "counts"
    shutil.copytree(TEST_DATA / "counts", counts)
    next(counts.glob("sampleC1_*")).unlink()
    with pytest.raises(ConfigError, match="sample sampleC1: expected 1 ReadsPerGene.out.tab file"):
        load_settings(test_config(counts_dir=str(counts)))


def test_missing_required_setting(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text("genome_fasta: x.fasta\n")
    with pytest.raises(ConfigError, match="missing: genome_gff, counts_dir, sample_sheet, contrasts"):
        load_settings(path)
