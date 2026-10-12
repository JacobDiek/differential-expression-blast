"""Run the full pipeline: DESeq2 -> pick unknown significant genes -> BLAST -> annotated tables.

    python run_pipeline.py                      # uses pipeline_config.yaml
    python run_pipeline.py my_other_config.yaml
    python run_pipeline.py --check               # check Python/R requirements are installed

    python run_pipeline.py quick-blast GENE_ID [GENE_ID ...]   # BLAST a few genes straight away
    python run_pipeline.py quick-blast --fasta sequences.fasta

The NCBI BLAST step is done on the NCBI website, so the first run stops after
writing the sequences to BLAST, with instructions. Run it again once the BLAST
results are downloaded to finish. Steps that are already done are skipped.
"""
import argparse
import sys
from datetime import datetime

from pipeline.errors import ConfigError
from pipeline.install_check import check_install


def step(number, title):
    print(f"\n=== Step {number}/5: {title} ===")


def run(config_path):
    # Imported here so --check still works when packages are missing
    from pipeline.blast import (BlastFolders, batch_status, genes_with_results, make_batch,
                                process_blast_results, read_batches)
    from pipeline.deseq2 import run_deseq2
    from pipeline.genes import (check_gene_ids_match, load_gene_table, load_significant_genes,
                                select_genes_to_blast)
    from pipeline.report import write_final_tables
    from pipeline.run_record import write_run_record
    from pipeline.settings import load_settings

    started = datetime.now()
    settings = load_settings(config_path)
    output_dir = settings["output_dir"]
    deseq2_dir = output_dir / "deseq2"
    final_dir = output_dir / "final_tables"
    folders = BlastFolders(output_dir / "ncbi")
    folders.make()

    step(1, "Differential expression (DESeq2)")
    run_deseq2(settings, deseq2_dir)

    step(2, "Finding significant genes with unknown names")
    print("Reading the GFF...")
    gene_table = load_gene_table(settings)
    significant = load_significant_genes(settings, deseq2_dir)
    check_gene_ids_match(gene_table, significant)
    genes_to_blast, too_long = select_genes_to_blast(settings, gene_table, significant)
    (folders.ncbi_dir / "genes_to_blast.txt").write_text("\n".join(genes_to_blast) + "\n")
    (folders.ncbi_dir / "genes_too_long.txt").write_text("\n".join(too_long) + "\n")
    print(f"{len(genes_to_blast)} unknown significant genes to BLAST "
          f"({len(too_long)} more skipped for being over {settings['max_gene_length']} bp)")

    step(3, "Reading downloaded BLAST results")
    blast_summary = process_blast_results(folders)
    have_results = genes_with_results(folders)

    step(4, "Preparing genes that still need BLAST")
    batches = read_batches(folders)
    in_a_batch = set().union(*batches.values()) if batches else set()
    need_batch = [g for g in genes_to_blast if g not in have_results and g not in in_a_batch]
    if need_batch:
        make_batch(settings, folders, need_batch)
        batches = read_batches(folders)
    else:
        print("Every gene to BLAST already has results or is in a batch.")
    pending, incomplete = batch_status(folders, batches, have_results)

    step(5, "Writing annotated tables")
    counts = write_final_tables(settings, significant, gene_table, blast_summary,
                                genes_to_blast, too_long, final_dir)
    for comparison, numbers in counts.items():
        print(f"\n{comparison}:")
        for label, value in numbers.items():
            print(f"  {label:<28}{value}")
    print(f"\nTables written to {final_dir} (CSV files + annotated_results.xlsx)")

    status = "waiting for BLAST results" if pending else "complete"
    record_path = write_run_record(started, config_path, settings, deseq2_dir, folders.xml_dir, counts, status)
    print(f"Run details saved to {record_path}")

    for batch_dir, missing in incomplete.items():
        print(f"\nNOTE: no BLAST results for {len(missing)} of the genes in {batch_dir.name} "
              f"(listed in {batch_dir / 'missing_results.txt'}).\n"
              f"Their chunk may not have been downloaded. To BLAST them again, delete the "
              f"{batch_dir.name} folder and rerun the pipeline.")

    if pending:
        print("\n" + "=" * 70)
        print("ACTION NEEDED: BLAST these genes on the NCBI website, then run the pipeline again.")
        print("=" * 70)
        for batch_dir in pending:
            print()
            print((batch_dir / "INSTRUCTIONS.txt").read_text())
    else:
        print("\nPipeline complete.")


def quick_blast_main(argv):
    parser = argparse.ArgumentParser(
        prog="run_pipeline.py quick-blast",
        description="BLAST a few genes or sequences directly through NCBI (one request each, "
                    "a minute or two per gene) and name them the same way as the pipeline.")
    parser.add_argument("genes", nargs="*", help="gene IDs to look up in the genome from the config")
    parser.add_argument("--fasta", help="FASTA file of your own sequences to BLAST")
    parser.add_argument("--config", default="pipeline_config.yaml",
                        help="config with the genome, GFF and output_dir (default: pipeline_config.yaml)")
    parser.add_argument("--email", help="email address for NCBI (or set NCBI_EMAIL)")
    args = parser.parse_args(argv)
    if not args.genes and not args.fasta:
        parser.error("give gene IDs and/or --fasta")

    from pipeline.quick_blast import run_quick_blast
    from pipeline.settings import read_settings
    run_quick_blast(read_settings(args.config), args.genes, args.fasta, args.email)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "quick-blast":
        try:
            quick_blast_main(sys.argv[2:])
        except ConfigError as error:
            print(f"\n{error}", file=sys.stderr)
            sys.exit(1)
        return

    parser = argparse.ArgumentParser(
        description="DESeq2 + BLAST annotation pipeline",
        epilog="To BLAST just a few genes: python run_pipeline.py quick-blast --help")
    parser.add_argument("config", nargs="?", default="pipeline_config.yaml",
                        help="pipeline config file (default: pipeline_config.yaml)")
    parser.add_argument("--check", action="store_true",
                        help="only check that the Python and R requirements are installed")
    args = parser.parse_args()

    if args.check:
        sys.exit(0 if check_install() else 1)

    try:
        run(args.config)
    except ImportError as error:
        print(f"\n{error}\nRun 'python run_pipeline.py --check' to see what is missing.", file=sys.stderr)
        sys.exit(1)
    except ConfigError as error:
        print(f"\n{error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
