import importlib.util
import shutil
import subprocess
import sys

MIN_PYTHON = (3, 11)
# import name -> pip package name
PYTHON_PACKAGES = {"Bio": "biopython", "pandas": "pandas", "numpy": "numpy",
                   "yaml": "pyyaml", "openpyxl": "openpyxl"}
R_PACKAGES = ["DESeq2", "ggplot2", "jsonlite"]


def check_install():
    # Prints what is installed and what is missing; returns True if everything is there
    ok = True

    def report(passed, label, fix=""):
        nonlocal ok
        ok = ok and passed
        print(f"  [{'OK' if passed else 'MISSING'}] {label}" + ("" if passed or not fix else f"\n            -> {fix}"))

    print("Python:")
    version = ".".join(map(str, sys.version_info[:3]))
    report(sys.version_info >= MIN_PYTHON, f"Python {version} (need {MIN_PYTHON[0]}.{MIN_PYTHON[1]}+)",
           "install a newer Python, then rerun ./setup.sh")
    for module, package in PYTHON_PACKAGES.items():
        report(importlib.util.find_spec(module) is not None, package, "pip install -r requirements.txt")

    print("R:")
    rscript = shutil.which("Rscript")
    report(rscript is not None, "Rscript", "install R from https://cran.r-project.org and make sure Rscript is on your PATH")
    if rscript:
        packages = ", ".join(f'"{p}"' for p in R_PACKAGES)
        result = subprocess.run(
            [rscript, "-e", f'for (p in c({packages})) cat(p, requireNamespace(p, quietly = TRUE), "\\n")'],
            capture_output=True, text=True)
        installed = dict(line.split() for line in result.stdout.splitlines() if len(line.split()) == 2)
        for package in R_PACKAGES:
            report(installed.get(package) == "TRUE", f"R package {package}", "Rscript R/install_packages.R")

    print("\nEverything is installed." if ok else "\nSome requirements are missing (see above).")
    return ok
