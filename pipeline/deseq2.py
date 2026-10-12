import json
import subprocess
import sys
from pathlib import Path

R_SCRIPT = Path(__file__).resolve().parent.parent / "R" / "run_deseq2.R"


def comparison_name(contrast):
    return f"{contrast[0]}_vs_{contrast[1]}"


def significant_genes_path(deseq2_dir, contrast):
    return deseq2_dir / f"significant_genes_{comparison_name(contrast)}.csv"


def run_deseq2(settings, deseq2_dir):
    deseq2_dir.mkdir(parents=True, exist_ok=True)

    r_settings = {
        "counts_dir": str(settings["counts_dir"]),
        "sample_sheet": str(settings["sample_sheet"]),
        "star_count_column": settings["star_count_column"],
        "exclude_samples": settings["exclude_samples"],
        "contrasts": settings["contrasts"],
        "padj_cutoff": settings["padj_cutoff"],
        "log2fc_cutoff": settings["log2fc_cutoff"],
        "output_dir": str(deseq2_dir),
    }

    # Skip if DESeq2 already ran with exactly these settings
    settings_path = deseq2_dir / "settings_used.json"
    outputs_exist = all(significant_genes_path(deseq2_dir, c).exists() for c in settings["contrasts"])
    if outputs_exist and settings_path.exists() and json.loads(settings_path.read_text()) == r_settings:
        print("DESeq2 results are up to date, skipping.")
        return

    # Written only after a successful run, so a failed run is retried next time
    settings_path.unlink(missing_ok=True)
    pending_path = deseq2_dir / "settings_pending.json"
    pending_path.write_text(json.dumps(r_settings, indent=2))

    print("Running DESeq2 (this takes a few minutes)...")
    # Flush first so our messages print before the R output, not after it
    sys.stdout.flush()
    result = subprocess.run(["Rscript", str(R_SCRIPT), str(pending_path)])
    if result.returncode != 0:
        raise RuntimeError("The DESeq2 R script failed, see the R output above.")

    pending_path.rename(settings_path)
