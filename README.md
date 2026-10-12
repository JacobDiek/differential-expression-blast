# Differential Expression BLAST

Annotating unknown genes that come up as significant in a differential expression analysis. Built for the brook stickleback (*Culaea inconstans*) genome, but works for any organism with the same kinds of input files.

Genome annotations often leave many genes without a real name (placeholder IDs like `egapxtmp_XXXXXX` or `LOC12345`). The pipeline:

1. Runs DESeq2 on RNA-seq counts to find differentially expressed genes for the comparisons you choose.
2. Picks out the significant genes that have no real name in the GFF.
3. Pulls each one's nucleotide sequence out of the genome and splits them into chunks for NCBI BLAST.
4. Reads the BLAST results and names each gene from its hits, or flags it as still unknown.
5. Writes one table per comparison with every significant gene, its statistics and its best annotation.

## Quick start

```bash
python run_pipeline.py
```

The first run does everything up to BLAST, then stops and prints instructions for uploading the chunk files to the NCBI BLAST website. After downloading the results into the folder it names, run the same command again to finish. Steps that are already done are skipped, so it is always safe to rerun.

Final tables end up in `results/final_tables/`: one CSV per comparison, plus `annotated_results.xlsx` with all of them.

## Setup

You need **Python 3.11+** (tested with 3.14) and **R** (tested with 4.5, from [CRAN](https://cran.r-project.org)) installed. Then, on macOS or Linux:

```bash
./setup.sh
source venv/bin/activate
```

`setup.sh` does three things:

1. Creates a Python environment in `venv/` and installs the Python packages.
2. Installs any missing R packages (DESeq2, ggplot2, jsonlite).
3. Checks that everything works.

You can rerun it at any time.

To check an existing install, run `python run_pipeline.py --check`. It lists anything missing and how to fix it.

<details>
<summary>Manual setup (e.g. on Windows)</summary>

```bash
python -m venv venv
venv\Scripts\activate          # macOS/Linux: source venv/bin/activate
pip install -r requirements.txt
Rscript R/install_packages.R
python run_pipeline.py --check
```

`Rscript R/install_packages.R --analysis` also installs the extra packages for `R/brook_stickleback_analysis.R`.
</details>

### Test that it works

```bash
python -m pytest
```

This runs about 30 tests in under 10 seconds, including the whole pipeline on a small made-up dataset in `tests/data/`. You can also try the pipeline on that dataset yourself:

```bash
python run_pipeline.py tests/data/test_config.yaml
```

Then copy `tests/data/blast_results/TEST0001-Alignment.xml` into `tests/data/output/ncbi/ncbi_files_large_alignment/` and run it again to see the second half.

## Inputs

| Input | Description |
|---|---|
| Genome FASTA | The genome assembly |
| Genome GFF | Its annotation. Gene rows need an ID attribute matching the count tables (e.g. `locus_tag`), and ideally `Name`, `description` and `gene_biotype` |
| Count tables | One STAR `<sample>_..._ReadsPerGene.out.tab` file per sample, all in one folder |
| Sample sheet | CSV with `sample` and `condition` columns. Each sample's count file must start with its sample name followed by `_` or `.` |

Example sample sheet:

```csv
sample,condition
25030R-01-01,liver_female
25030R-01-06,liver_male
25030R-01-11,testes
```

## Configuration

Everything is set in [`pipeline_config.yaml`](pipeline_config.yaml). Paths are relative to the config file. To run a different organism or analysis, copy it, edit it and run `python run_pipeline.py my_config.yaml`. Give each run its own `output_dir`.

| Setting | Meaning |
|---|---|
| `genome_fasta`, `genome_gff`, `counts_dir`, `sample_sheet` | Input files (see above) |
| `output_dir` | Where all results go (default `results`) |
| `gene_id_attribute` | GFF gene attribute holding the IDs used in the count tables (`locus_tag` for brook stickleback) |
| `star_count_column` | STAR column to use: `2` unstranded, `3` forward, `4` reverse stranded, or `auto`. The pipeline checks the count files and warns if the number looks wrong for your library |
| `exclude_samples` | Samples left out of DESeq2, e.g. outliers on `deseq2/pca.png` |
| `contrasts` | Comparisons as `[condition, reference]`. Positive log2FC = higher in the first condition |
| `padj_cutoff`, `log2fc_cutoff` | Significance cutoffs (`padj <` and `abs(log2FC) >`) |
| `unknown_name_patterns` | A gene counts as unknown if its GFF `Name` is missing, equal to its ID, or matches one of these patterns (default `^LOC\d+$`) |
| `blast_biotypes` | Only BLAST these `gene_biotype` values, e.g. `[protein_coding]`. Empty = all |
| `max_gene_length` | Genes longer than this (bp) are not BLASTed (default 20,000) |
| `chunk_max_bases`, `chunk_max_seqs` | Size limits for each file uploaded to NCBI |

## The BLAST step

Thousands of genes are too many to send to NCBI one by one, so they are BLASTed through the website in batches:

1. The pipeline writes a numbered batch to `results/ncbi/blast_batches/batch_N/`, with the chunk files, the gene list and an `INSTRUCTIONS.txt`.
2. Upload each chunk to [NCBI Nucleotide BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn) (database **core_nt**, program **megablast**, the defaults) and download the results as **XML**.
3. Put the downloaded `...-Alignment.xml` files in `results/ncbi/ncbi_files_large_alignment/` and rerun the pipeline.

Genes that already have BLAST results are never BLASTed again, so changing the comparisons or cutoffs later only creates a batch for the new genes. If some genes in a batch get no results (e.g. a chunk was missed), the pipeline lists them in that batch's `missing_results.txt`. Delete the batch folder and rerun to BLAST them again.

**How genes are named:**

- **Organism names are removed** from each hit's description. This includes prefixes like `PREDICTED:`, "sp." names, hybrids, and three-part subspecies names such as *Seriola lalandi dorsalis*, which are learned from the hits themselves.
- **Noise:** hits from genomic DNA records (chromosome, scaffold, clone, contig, ultra-conserved elements, etc.) are skipped.
- **Uninformative:** hits saying "uncharacterized", "hypothetical" or "unknown" are counted as uninformative.
- **Unknown:** if there are no informative hits, or more than a quarter of the remaining hits are uninformative, the gene stays unknown.
- **Weighted vote:** otherwise the informative hits vote on the name.
  - Each hit's vote is weighted by -log10(e-value), capped at 180, so strong hits count much more than weak ones.
  - Versions of the same name are counted together ("... transcript variant X2", "... -like", "... homolog").
  - The group with the most weight wins, and is named after its best-ranked hit.
- **Name support:** the winner's share of the total weight. Low support (e.g. under 50%) means the hits disagree, so the name deserves a manual look.
- **Per-gene files:** each file in `identified_genes/` shows the full vote breakdown and the single best hit for comparison.

## BLASTing a few genes (quick-blast)

To check a handful of genes without the website upload, `quick-blast` sends them straight to NCBI BLAST:

```bash
python run_pipeline.py quick-blast egapxtmp_018614 egapxtmp_009144   # gene IDs from the GFF
python run_pipeline.py quick-blast --fasta my_sequences.fasta         # your own sequences
```

- **Settings:** it uses the same database (core_nt) and program (megablast) as the website step, and the same weighted-vote naming. It prints each gene's name and support.
- **Config:** gene IDs are looked up with the genome, GFF and `gene_id_attribute` from `pipeline_config.yaml` (use `--config` for another file). The count tables and sample sheet aren't needed.
- **Output:** results go to `results/quick_blast/`, with `summary.csv` and the same per-gene files as the main pipeline.
- **Speed:** each gene is a separate NCBI request, which takes a minute or two, so this is for a few genes, not thousands.
- **Email:** NCBI asks for an email address. Pass `--email you@example.com`, set `export NCBI_EMAIL=you@example.com`, or type it when asked.

## Outputs

```
results/
├── deseq2/
│   ├── significant_genes_<comparison>.csv   # padj and log2FC cutoffs applied
│   ├── all_genes_<comparison>.csv           # full DESeq2 results
│   ├── pca.png, pca_coordinates.csv         # for spotting outlier samples
│   └── volcano_<comparison>.png
├── ncbi/
│   ├── genes_to_blast.txt, genes_too_long.txt
│   ├── blast_batches/batch_N/               # chunks to upload + INSTRUCTIONS.txt
│   ├── ncbi_files_large_alignment/          # downloaded BLAST XMLs go here
│   ├── ncbi_files_raw_output/               # one XML per gene
│   ├── ncbi_files_parsed_output/            # readable hit list per gene
│   ├── identified_genes/                    # <name>__<gene>.txt (with vote breakdown) or UNKNOWN_<gene>.txt
│   └── blast_summary.csv                    # one row per BLASTed gene
├── final_tables/
│   ├── annotated_<comparison>.csv
│   └── annotated_results.xlsx               # all comparisons + a summary sheet
└── run_info/
    ├── run_<date_time>.json                 # details of each run
    ├── config_<date_time>.yaml              # the config used for that run
    └── latest_run.json
```

Each final table has every significant gene sorted by padj, with:

- `best_annotation` and `annotation_source`: the GFF description if it is informative, otherwise the BLAST name. "Uncharacterized protein LOC…"-style descriptions don't count as informative.
- The DESeq2 statistics.
- The GFF name, description and biotype.
- `blast_status`: `identified`, `unknown`, `waiting for BLAST results`, `not needed (named gene)` or `not BLASTed (too long / biotype excluded)`.
- `blast_top_hit` (the weighted-vote name), `blast_name_support` (0–1) and `blast_best_single_hit`.

The Excel file has the same tables, one sheet per comparison. Each sheet has filters and a frozen header row, and names with under 50% support are highlighted. A Summary sheet lists the counts per comparison and explains each column.

**Run details:** every run saves a record in `run_info/`, so results can be traced back and reported. It includes:

- the config and resolved settings;
- the Python, package, R and DESeq2 versions;
- the code's git commit;
- input file sizes, dates and checksums;
- the BLAST program, version and database read from the downloaded XML files (e.g. `core_nt`), with their download dates;
- the result counts.

## Project layout

```
differential-expression-blast/
├── run_pipeline.py            # the tool: full pipeline, quick-blast and --check
├── pipeline_config.yaml       # pipeline settings (set up for the brook stickleback data)
├── setup.sh                   # one-time install (macOS/Linux)
├── requirements.txt           # Python packages
├── pipeline/                  # pipeline steps: settings, DESeq2, gene selection, BLAST batches,
│                              #   tables, run record, quick-blast, install check
├── ncbi_tools/                # BLAST result handling: FASTA writing, chunking, XML splitting,
│                              #   parsing and gene naming
├── core/                      # GFF/genome lookup and online BLAST requests
├── R/
│   ├── run_deseq2.R           # DESeq2 step used by the pipeline
│   ├── install_packages.R     # installs the R packages
│   └── brook_stickleback_analysis.R   # this project's exploratory analysis (see below)
├── tests/                     # automated tests + small synthetic test dataset
├── data/                      # inputs (not in git)
└── results/                   # outputs (not in git)
```

`data/` and `results/` are too large for git (the brook stickleback genome alone is ~1 GB), so they are git-ignored and need to be shared separately.
