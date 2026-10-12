import os

from Bio import SeqIO

output_prefix = "unknown_chunk"


def write_chunk(records, num, output_dir):
    filename = os.path.join(output_dir, f"{output_prefix}_{num}.fasta")
    with open(filename, "w") as f:
        for record in records:
            f.write(f">{record.id}\n")
            f.write(f"{str(record.seq)}\n")
    total_bp = sum(len(r.seq) for r in records)
    print(f"Wrote {filename}: {len(records)} sequences, {total_bp} bp")


def fasta_file_splitter(input_file, output_dir, max_bases_per_chunk=1_000_000, max_seqs_per_chunk=30):
    os.makedirs(output_dir, exist_ok=True)

    records = list(SeqIO.parse(input_file, "fasta"))

    chunk = []
    chunk_bases = 0
    chunk_num = 1

    for record in records:
        seq_len = len(record.seq)

        if seq_len > max_bases_per_chunk:
            print(f"WARNING: {record.id} is {seq_len} bp, exceeds chunk limit on its own.")

        # Flush current chunk if adding this record would exceed EITHER limit
        would_exceed_bases = chunk_bases + seq_len > max_bases_per_chunk
        would_exceed_seqs = len(chunk) + 1 > max_seqs_per_chunk

        if (would_exceed_bases or would_exceed_seqs) and chunk:
            write_chunk(chunk, chunk_num, output_dir)
            chunk_num += 1
            chunk = []
            chunk_bases = 0

        chunk.append(record)
        chunk_bases += seq_len

    if chunk:
        write_chunk(chunk, chunk_num, output_dir)

    print(f"\nDone. Created {chunk_num} chunk files in '{output_dir}/'.")
