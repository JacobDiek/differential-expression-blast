import pandas as pd

COLUMN_NAMES = {2: "unstranded", 3: "forward stranded", 4: "reverse stranded"}

# Reads on genes in the forward (column 3) vs reverse (column 4) STAR columns:
# a stranded library has far more in one of them, an unstranded one about the same
STRANDED_RATIO = 3
UNSTRANDED_RATIO = 1.5


def detect_star_count_column(count_files, max_files=3):
    # Returns (column or None if unclear, forward reads, reverse reads)
    forward = reverse = 0
    for count_file in count_files[:max_files]:
        table = pd.read_csv(count_file, sep="\t", header=None, skiprows=4)
        forward += table[2].sum()
        reverse += table[3].sum()

    low, high = sorted([forward, reverse])
    if low == 0:
        return (3 if forward > reverse else 4) if high > 0 else None, forward, reverse

    ratio = high / low
    if ratio >= STRANDED_RATIO:
        return (3 if forward > reverse else 4), forward, reverse
    if ratio <= UNSTRANDED_RATIO:
        return 2, forward, reverse
    return None, forward, reverse
