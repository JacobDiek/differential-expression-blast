import os
import re


class Hit:
    def __init__(self, hit_def):
        self.hit_def = hit_def
        
        self.clean_def = self.hit_def.removeprefix("PREDICTED: ")

        self.parenthetical = self._extract_parenthetical()
        self.classification = self._classify()
        self.description = self._extract_description()
        self.species = self._extract_species()
        self.normalized_name = self._normalize_description()

    def _extract_parenthetical(self):
        if "(" not in self.clean_def:
            return None
        
        open_pos = self.clean_def.find("(")
        close_pos = self.clean_def.find(")")

        return self.clean_def[open_pos+1:close_pos]
    
    def _classify(self):
        if self.parenthetical is None:
            return "no_tag"
        elif "LOC" in self.parenthetical:
            return "locus_tag"
        else:
            return "gene_tag"
        
    def _extract_description(self):
        if self.classification == "no_tag":
            return self.clean_def
        return self.clean_def[:(self.clean_def.find("(")-1)]
    
    def _extract_species(self):
        words = self.clean_def.split()
        if len(words) >= 2:
            return f"{words[0]} {words[1]}"
        return "Unknown"

    def _normalize_description(self):
        name = self.description

        if self.species != "Unknown":
            name = name.replace(self.species, "")

        name = name.lower()
        name = name.replace("-", " ")

        name = " ".join(name.split())

        return name


def sanitize_filename(name, max_length=100):
    name = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
    if not name:
        name = "unnamed"
    return name[:max_length]


def extract_locus_tag(input_filename):
    name = os.path.splitext(input_filename)[0]
    if name.startswith("blast_output_"):
        name = name[len("blast_output_"):]
    return name


def gene_identifier(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    input_files = os.listdir(input_dir)

    # Hard Skips: Genomic Noise
    genomic_noise = ["chromosome", "haplotype", "scaffold", "assembly", "clone"]
    
    # Soft Flags: Uknown but Real Genes
    uninformative = ["uncharacterized", "hypothetical", "unknown"]

    for input_file in input_files:
        input_file_path = os.path.join(input_dir, input_file)
        
        with open(input_file_path, 'r', errors="replace") as f:
            lines = [line.rstrip("\n") for line in f]

        all_hit_defs = lines[4::5]

        locus_tag = extract_locus_tag(input_file)

        file_hits = []
        file_flags = []
        for line in all_hit_defs:
            individual_hit = Hit(hit_def=line)
            normalized_hit = individual_hit.normalized_name
            
            # Skip genomic noise entirely
            if any(noise in normalized_hit for noise in genomic_noise):
                continue
                
            # Flags number of unknown genes
            if any(flag in normalized_hit for flag in uninformative):
                file_flags.append(normalized_hit)
                continue
                
            file_hits.append(normalized_hit)

        if len(file_flags) * 3 > len(file_hits):
            output_filename = f"UNKNOWN_{sanitize_filename(locus_tag)}.txt"
            output_content = (
                f"Locus tag: {locus_tag}\n"
                f"Status: Unknown Gene\n"
                f"Flagged hits ({len(file_flags)}): {sorted(set(file_flags))}\n"
                f"Informative hits ({len(file_hits)}): {sorted(set(file_hits))}\n"
            )
        else:
            unique_hits = sorted(set(file_hits))
            top_hit = unique_hits[0] if unique_hits else "unnamed"
            output_filename = f"{sanitize_filename(top_hit)}_{sanitize_filename(locus_tag)}.txt"
            output_content = (
                f"Locus tag: {locus_tag}\n"
                f"Status: Identified\n"
                f"Hits: {unique_hits}\n"
            )

        output_path = os.path.join(output_dir, output_filename)
        with open(output_path, "w", encoding="utf-8") as out_f:
            out_f.write(output_content)


gene_identifier(input_dir="/Users/jacobdiek/Coding/projects/differential-expression-blast/ncbi_files_results/ncbi_files_parsed_output",
                output_dir="/Users/jacobdiek/Coding/projects/differential-expression-blast/ncbi_files_results/output_dir")