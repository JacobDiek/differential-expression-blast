import os
import shutil
from pathlib import Path

import pandas as pd
import yaml

from pipeline.errors import ConfigError
from pipeline.strandedness import COLUMN_NAMES, detect_star_count_column


DEFAULTS = {
    "output_dir": "results",
    "gene_id_attribute": "locus_tag",
    "star_count_column": 4,
    "exclude_samples": [],
    "padj_cutoff": 0.05,
    "log2fc_cutoff": 1,
    "unknown_name_patterns": [r"^LOC\d+$"],
    "blast_biotypes": [],
    "max_gene_length": 20_000,
    "chunk_max_bases": 1_000_000,
    "chunk_max_seqs": 30,
}

REQUIRED = ["genome_fasta", "genome_gff", "counts_dir", "sample_sheet", "contrasts"]
PATH_SETTINGS = ["genome_fasta", "genome_gff", "counts_dir", "sample_sheet", "output_dir"]


def load_settings(config_path):
    # Settings for the full pipeline, with every input checked
    settings = read_settings(config_path, required=REQUIRED)
    validate_settings(settings)
    return settings


def read_settings(config_path, required=()):
    # Config + defaults, with paths resolved; only the `required` settings are checked
    config_path = Path(config_path).resolve()
    if not config_path.exists():
        raise ConfigError(f"Config file not found: {config_path}")

    with open(config_path) as f:
        user_settings = yaml.safe_load(f) or {}

    missing = [key for key in required if not user_settings.get(key)]
    if missing:
        raise ConfigError(f"{config_path.name} is missing: {', '.join(missing)}")

    settings = {**DEFAULTS, **user_settings}
    settings["exclude_samples"] = settings["exclude_samples"] or []
    settings["blast_biotypes"] = settings["blast_biotypes"] or []
    settings["unknown_name_patterns"] = settings["unknown_name_patterns"] or []

    # Paths in the config are relative to the config file's folder
    for key in PATH_SETTINGS:
        if settings.get(key):
            settings[key] = (config_path.parent / settings[key]).resolve()

    return settings


def find_count_file(sample, count_files):
    # Same rule as R/run_deseq2.R: file name starts with the sample name followed by "_" or "."
    return [f for f in count_files if f.startswith(f"{sample}_") or f.startswith(f"{sample}.")]


def validate_settings(settings):
    problems = []

    for key in ["genome_fasta", "genome_gff", "sample_sheet"]:
        if not settings[key].is_file():
            problems.append(f"{key}: file not found ({settings[key]})")
    if not settings["counts_dir"].is_dir():
        problems.append(f"counts_dir: folder not found ({settings['counts_dir']})")

    if settings["star_count_column"] not in (2, 3, 4, "auto"):
        problems.append("star_count_column must be 2 (unstranded), 3 (forward), 4 (reverse) or auto")

    if shutil.which("Rscript") is None:
        problems.append("Rscript was not found. Install R and make sure Rscript is on your PATH")

    contrasts = settings["contrasts"]
    if not isinstance(contrasts, list) or not all(isinstance(c, list) and len(c) == 2 for c in contrasts):
        problems.append("contrasts must be a list of [condition, reference] pairs")
        contrasts = []

    if settings["sample_sheet"].is_file() and settings["counts_dir"].is_dir():
        sample_sheet = pd.read_csv(settings["sample_sheet"], dtype=str)
        sample_sheet.columns = sample_sheet.columns.str.strip()

        if not {"sample", "condition"} <= set(sample_sheet.columns):
            problems.append("sample_sheet must have 'sample' and 'condition' columns")
        else:
            samples = sample_sheet["sample"].str.strip()
            count_files = [f for f in os.listdir(settings["counts_dir"]) if f.endswith("ReadsPerGene.out.tab")]

            if samples.duplicated().any():
                problems.append(f"sample_sheet has duplicate samples: {sorted(set(samples[samples.duplicated()]))}")

            for sample in samples:
                matches = find_count_file(sample, count_files)
                if len(matches) != 1:
                    problems.append(f"sample {sample}: expected 1 ReadsPerGene.out.tab file in counts_dir, found {len(matches)}")

            unknown_excludes = set(settings["exclude_samples"]) - set(samples)
            if unknown_excludes:
                problems.append(f"exclude_samples not in sample_sheet: {sorted(unknown_excludes)}")

            kept = sample_sheet[~samples.isin(settings["exclude_samples"])]
            kept_files = [settings["counts_dir"] / matches[0] for sample in kept["sample"].str.strip()
                          if len(matches := find_count_file(sample, count_files)) == 1]
            if kept_files and settings["star_count_column"] in (2, 3, 4, "auto"):
                problems += check_strandedness(settings, kept_files)
            conditions = set(kept["condition"].str.strip())
            for condition, reference in contrasts:
                for name in (condition, reference):
                    if name not in conditions:
                        problems.append(f"contrast [{condition}, {reference}]: no samples with condition '{name}'")

    if problems:
        raise ConfigError("Problems with the config:\n  - " + "\n  - ".join(problems))


def check_strandedness(settings, count_files):
    # The wrong STAR column silently gives bad counts, so compare the setting with what the
    # count files look like. Returns config problems; a mismatch is only a warning.
    detected, forward, reverse = detect_star_count_column(count_files)
    reads = f"{forward:,} reads on genes forward vs {reverse:,} reverse"
    setting = settings["star_count_column"]

    if setting == "auto":
        if detected is None:
            return [f"star_count_column: auto could not tell the strandedness ({reads}). Set it to 2, 3 or 4"]
        settings["star_count_column"] = detected
        print(f"Count files look {COLUMN_NAMES[detected]} ({reads}), using STAR column {detected}.")
    elif detected is None:
        print(f"WARNING: could not confirm star_count_column {setting} from the count files ({reads}).")
    elif detected != setting:
        print(f"WARNING: star_count_column is {setting} ({COLUMN_NAMES[setting]}), but the count files look "
              f"{COLUMN_NAMES[detected]} ({reads}). Check your library prep; column {detected} may be correct.")
    return []
