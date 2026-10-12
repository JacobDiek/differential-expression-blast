import hashlib
import json
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime
from importlib import metadata
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PYTHON_PACKAGES = ["biopython", "pandas", "numpy", "PyYAML", "openpyxl"]


def package_versions():
    versions = {"python": platform.python_version()}
    for package in PYTHON_PACKAGES:
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def code_version():
    # Git commit of this tool, and whether it had uncommitted changes
    def git(*args):
        result = subprocess.run(["git", "-C", str(PROJECT_ROOT), *args], capture_output=True, text=True)
        return result.stdout.strip() if result.returncode == 0 else None

    if shutil.which("git") is None:
        return None
    commit = git("rev-parse", "HEAD")
    if commit is None:
        return None
    return {"commit": commit, "uncommitted_changes": bool(git("status", "--porcelain"))}


def file_info(path):
    path = Path(path)
    if path.is_dir():
        files = [f for f in path.iterdir() if f.is_file()]
        return {"path": str(path), "files": len(files)}
    stat = path.stat()
    info = {"path": str(path), "bytes": stat.st_size,
            "modified": datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds")}
    if stat.st_size < 10_000_000:
        info["md5"] = hashlib.md5(path.read_bytes()).hexdigest()
    return info


def header_field(header, tag):
    match = re.search(f"<BlastOutput_{tag}>([^<]*)<", header)
    return match.group(1) if match else None


def blast_info(xml_dir):
    # Program, version and database from the header of every downloaded BLAST XML
    searches = {}
    dates = []
    for xml_file in Path(xml_dir).glob("*.xml"):
        with open(xml_file, errors="replace") as f:
            header = f.read(2000)
        fields = tuple(header_field(header, tag) for tag in ("program", "version", "db"))
        searches[fields] = searches.get(fields, 0) + 1
        dates.append(xml_file.stat().st_mtime)

    if not dates:
        return {"downloaded_xml_files": 0}
    return {
        "downloaded_xml_files": len(dates),
        "searches": [{"program": p, "version": v, "database": db, "files": n}
                     for (p, v, db), n in searches.items()],
        "files_dated": [datetime.fromtimestamp(min(dates)).date().isoformat(),
                        datetime.fromtimestamp(max(dates)).date().isoformat()],
    }


def write_run_record(started, config_path, settings, deseq2_dir, xml_dir, counts, status):
    run_info_dir = settings["output_dir"] / "run_info"
    run_info_dir.mkdir(parents=True, exist_ok=True)
    stamp = started.strftime("%Y-%m-%d_%H%M%S")

    deseq2_versions_file = deseq2_dir / "versions.json"
    record = {
        "started": started.isoformat(timespec="seconds"),
        "finished": datetime.now().isoformat(timespec="seconds"),
        "status": status,
        "command": " ".join([Path(sys.executable).name] + sys.argv),
        "config_file": str(Path(config_path).resolve()),
        "code_version": code_version(),
        "software": package_versions(),
        "deseq2_run_with": json.loads(deseq2_versions_file.read_text()) if deseq2_versions_file.exists() else None,
        "settings": {k: str(v) if isinstance(v, Path) else v for k, v in settings.items()},
        "inputs": {key: file_info(settings[key])
                   for key in ("genome_fasta", "genome_gff", "sample_sheet", "counts_dir")},
        "blast": blast_info(xml_dir),
        "results": counts,
    }

    record_path = run_info_dir / f"run_{stamp}.json"
    record_path.write_text(json.dumps(record, indent=2, default=str))
    shutil.copy(record_path, run_info_dir / "latest_run.json")
    shutil.copy(config_path, run_info_dir / f"config_{stamp}.yaml")
    return record_path
