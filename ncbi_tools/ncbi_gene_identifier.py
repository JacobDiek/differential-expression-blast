import math
import os
import re


# Record-type prefixes NCBI puts in front of the organism name
DEF_PREFIXES = ("PREDICTED: ", "TPA_asm: ", "TPA_inf: ", "TPA_exp: ", "TPA: ")

# Words that mean the next word is also part of the organism name
SPECIES_JOINERS = {"sp.", "spp.", "cf.", "aff.", "subsp.", "var.", "x"}


def strip_def_prefixes(hit_def):
    hit_def = hit_def.strip()
    for prefix in DEF_PREFIXES:
        hit_def = hit_def.removeprefix(prefix)
    # e.g. "[Epinephelus coioides] Epinephelus coioides EcSCYA113 mRNA ..."
    return re.sub(r"^\[[^\]]*\]\s*", "", hit_def)


def species_words(words, multi_word_species=frozenset()):
    # Number of leading words that make up the organism name: "Genus species",
    # plus e.g. "sp. 'strain'", "x Genus species" (hybrids) or a known subspecies
    if len(words) < 2 or not words[0][:1].isupper():
        return 0

    n = 2
    while n < len(words):
        if words[n - 1] in SPECIES_JOINERS:
            n += 1
        elif words[n] in SPECIES_JOINERS:
            n += 1 if words[n] != "x" else 3
        elif " ".join(words[:n + 1]) in multi_word_species:
            n += 1
        else:
            break
    return min(n, len(words))


def simplify_description(text):
    # Description without gene symbol, transcript variant, mRNA etc., for comparing between organisms
    text = re.split(r"\s\(", text)[0].lower().replace("-", " ")
    text = re.sub(r",.*$", "", text)
    return " ".join(text.split())


def find_multi_word_species(hit_defs):
    # Some organisms have three-part names (e.g. the subspecies "Seriola lalandi dorsalis").
    # A third word counts as part of the organism name when it follows that organism in
    # most of its hits, and removing it leaves a gene description that other organisms'
    # hits also have (while keeping it doesn't). Ordinary description words like
    # "putative" or "genomic", which start other organisms' descriptions, never count.
    organisms_with = {}   # simplified description -> organisms whose hits have it
    hits_per_organism = {}
    first_word_organisms = {}  # first description word -> organisms it follows
    candidates = {}       # "Genus species word" -> descriptions that follow it
    for hit_def in set(hit_defs):
        words = strip_def_prefixes(hit_def).split()
        n = species_words(words)
        if n == 0 or len(words) < n + 2:
            continue
        organism = " ".join(words[:n])
        hits_per_organism[organism] = hits_per_organism.get(organism, 0) + 1
        organisms_with.setdefault(simplify_description(" ".join(words[n:])), set()).add(organism)
        third = words[n]
        first_word_organisms.setdefault(third, set()).add(organism)
        if third.isalpha() and third.islower():
            candidates.setdefault(f"{organism} {third}", []).append(" ".join(words[n + 1:]))

    multi_word_species = set()
    for name, rests in candidates.items():
        organism, third = name.rsplit(" ", 1)
        if first_word_organisms[third] - {organism}:
            continue
        fits = 0
        for rest in rests:
            without_word = organisms_with.get(simplify_description(rest), set()) - {organism}
            with_word = organisms_with.get(simplify_description(f"{third} {rest}"), set()) - {organism}
            if without_word and not with_word:
                fits += 1
        most_of_its_hits = len(rests) * 2 >= hits_per_organism[organism]
        if len(rests) >= 2 and most_of_its_hits and fits * 2 >= len(rests):
            multi_word_species.add(name)
    return multi_word_species


class Hit:
    def __init__(self, hit_def, multi_word_species=frozenset()):
        self.hit_def = hit_def
        self.multi_word_species = multi_word_species

        self.clean_def = strip_def_prefixes(self.hit_def)

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
        n = species_words(words, self.multi_word_species)
        if n:
            return " ".join(words[:n])
        return "Unknown"

    def _normalize_description(self):
        name = self.description

        if self.species != "Unknown":
            name = " ".join(name.split()).removeprefix(self.species)

        name = name.lower()
        name = name.replace("-", " ")

        name = " ".join(name.split())

        return name


def sanitize_filename(name, max_length=100):
    name = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
    if not name:
        name = "unnamed"
    return name[:max_length]


def unique_in_order(items):
    return list(dict.fromkeys(items))


# E-values of 0 mean "smaller than about 1e-180", so weights are capped there
MAX_HIT_WEIGHT = 180


def hit_weight(evalue_text):
    # Stronger hits (lower e-value) get a bigger vote: weight = -log10(evalue), at least 1
    try:
        evalue = float(evalue_text)
    except (TypeError, ValueError):
        return 1.0
    if evalue <= 0:
        return MAX_HIT_WEIGHT
    return min(max(-math.log10(evalue), 1.0), MAX_HIT_WEIGHT)


def vote_key(name):
    # Groups versions of the same gene name so they vote together, e.g.
    # "... like transcript variant x2", "... homolog" and "... like, mrna" -> "..."
    name = re.sub(r"\btranscript variant \w+", "", name)
    name = re.sub(r"\b(partial mrna|mrna|complete cds|partial cds)\b", "", name)
    name = " ".join(name.replace(",", " ").split())
    name = re.sub(r"( like| homolog)+$", "", name)
    return name


def weighted_vote(hits):
    # hits: (normalized name, weight) in BLAST order (best first).
    # The group of names with the most total weight wins; it is named after
    # its best-ranked hit. Support = the winner's share of all the weight.
    group_weight = {}
    group_name = {}
    for name, weight in hits:
        key = vote_key(name)
        group_weight[key] = group_weight.get(key, 0) + weight
        group_name.setdefault(key, name)

    # max() keeps the first key on ties, i.e. the group with the best-ranked hit
    winner = max(group_weight, key=group_weight.get)
    support = group_weight[winner] / sum(group_weight.values())
    return group_name[winner], support, group_weight


def extract_locus_tag(input_filename):
    name = os.path.splitext(input_filename)[0]
    if name.startswith("blast_output_"):
        name = name[len("blast_output_"):]
    return name


def read_hits(parsed_file):
    with open(parsed_file, 'r', errors="replace") as f:
        lines = [line.rstrip("\n") for line in f]

    # Each hit is "Hit Num: N", hit id, hit definition, "Evalue: x", blank line
    return [(lines[i + 2], lines[i + 3].removeprefix("Evalue: ") if i + 3 < len(lines) else None)
            for i, line in enumerate(lines)
            if line.startswith("Hit Num:") and i + 2 < len(lines)]


def gene_identifier(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    input_files = os.listdir(input_dir)

    # Hard Skips: Genomic Noise
    genomic_noise = ["chromosome", "haplotype", "scaffold", "assembly", "clone",
                     "genomic sequence", "ultra conserved element", "contig"]
    
    # Soft Flags: Uknown but Real Genes
    uninformative = ["uncharacterized", "hypothetical", "unknown"]

    # First pass: read every file, so organism names can be learned from all hits together
    hits_by_file = {input_file: read_hits(os.path.join(input_dir, input_file)) for input_file in input_files}
    multi_word_species = find_multi_word_species(
        hit_def for all_hits in hits_by_file.values() for hit_def, _ in all_hits)

    summary = []
    for input_file, all_hits in hits_by_file.items():
        locus_tag = extract_locus_tag(input_file)

        file_hits = []
        file_flags = []
        hit_votes = []
        for line, evalue in all_hits:
            individual_hit = Hit(hit_def=line, multi_word_species=multi_word_species)
            normalized_hit = individual_hit.normalized_name
            
            # Skip genomic noise
            if any(noise in normalized_hit for noise in genomic_noise):
                continue
                
            # Flags number of unknown genes
            if any(flag in normalized_hit for flag in uninformative):
                file_flags.append(normalized_hit)
                continue
                
            file_hits.append(normalized_hit)
            hit_votes.append((normalized_hit, hit_weight(evalue)))

        # Genes with no informative hits at all (e.g. only genomic noise) are also unknown
        if not file_hits or len(file_flags) * 3 > len(file_hits):
            top_hit = None
            best_single_hit = None
            support = None
            output_filename = f"UNKNOWN_{sanitize_filename(locus_tag)}.txt"
            output_content = (
                f"Locus tag: {locus_tag}\n"
                f"Status: Unknown Gene\n"
                f"Flagged hits ({len(file_flags)}): {unique_in_order(file_flags)}\n"
                f"Informative hits ({len(file_hits)}): {unique_in_order(file_hits)}\n"
            )
        else:
            # BLAST lists hits best (lowest evalue) first, so keep that order
            unique_hits = unique_in_order(file_hits)
            best_single_hit = unique_hits[0]
            top_hit, support, group_weight = weighted_vote(hit_votes)
            votes = sorted(group_weight.items(), key=lambda item: -item[1])
            vote_lines = "".join(f"  {weight / sum(group_weight.values()):6.1%}  {key}\n" for key, weight in votes)
            # Locus tag in the filename so genes with the same top hit don't overwrite each other
            output_filename = f"{sanitize_filename(top_hit)}__{locus_tag}.txt"
            output_content = (
                f"Locus tag: {locus_tag}\n"
                f"Status: Identified\n"
                f"Name: {top_hit} (weighted vote, {support:.0%} support)\n"
                f"Best single hit: {best_single_hit}\n"
                f"Votes (share of -log10 evalue weight):\n{vote_lines}"
                f"Hits: {unique_hits}\n"
            )

        output_path = os.path.join(output_dir, output_filename)
        with open(output_path, "w", encoding="utf-8") as out_f:
            out_f.write(output_content)

        summary.append({
            "gene_id": locus_tag,
            "blast_status": "unknown" if top_hit is None else "identified",
            "blast_top_hit": top_hit,
            "name_support": None if support is None else round(support, 3),
            "best_single_hit": best_single_hit,
            "informative_hits": len(file_hits),
            "uninformative_hits": len(file_flags),
        })

    # One row per gene, used by the pipeline for the final tables
    return summary
