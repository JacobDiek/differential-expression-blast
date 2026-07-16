import os
import xml.etree.ElementTree as ET
from pathlib import Path


def ncbi_mass_parser(input_dir, output_dir):
    raw_file_location = Path(input_dir)
    parsed_file_location = Path(output_dir)
    parsed_file_location.mkdir(parents=True, exist_ok=True)

    ncbi_raw_dir = os.listdir(raw_file_location)

    for file in ncbi_raw_dir:
        if not file.endswith(".xml"):
            continue

        full_file_path = raw_file_location / file
        tree = ET.parse(full_file_path)
        root = tree.getroot()

        file_path = parsed_file_location / f"blast_output_{file}.txt"

        with open(file_path, "w", encoding="utf-8") as f:
            for hit in root.iter("Hit"):
                hit_num = hit.findtext("./Hit_num")
                hit_id = hit.findtext("./Hit_id")
                hit_def = hit.findtext("./Hit_def")
                hsp_evalue = hit.findtext("./Hit_hsps/Hsp/Hsp_evalue")

                f.write(f"Hit Num: {hit_num}\n")
                f.write(f"{hit_id}\n")
                f.write(f"{hit_def}\n")
                f.write(f"Evalue: {hsp_evalue}\n\n")

ncbi_mass_parser(input_dir='/Users/jacobdiek/Coding/projects/differential-expression-blast/ncbi_files_results/ncbi_files_raw_output',
                 output_dir='/Users/jacobdiek/Coding/projects/differential-expression-blast/ncbi_files_results/ncbi_files_parsed_output')