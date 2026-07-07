import xml.etree.ElementTree as ET
from pathlib import Path


def blast_parse(gene_list: dict, file_location: Path):

    raw_file_location = file_location / "raw_results"
    parsed_file_location = file_location / "parsed_results"

    for gene in gene_list:

        file_path_tree = raw_file_location / f"blast_{gene}.xml"
        tree = ET.parse(file_path_tree)
        root = tree.getroot()

        file_path = parsed_file_location / f"blast_output_{gene}.txt"
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
