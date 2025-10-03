#!/usr/bin/env python
import argparse
import pandas as pd
import pysam
from tqdm import tqdm
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count
import time
import threading

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Generate FASTA file")
parser.add_argument("--bam_file", required=True, help="Path to BAM file")
parser.add_argument("--modkit_path", required=True, help="Path to modkit TSV file")
parser.add_argument("--ref_genome_path", required=True, help="Path to reference genome")
parser.add_argument("output_fasta", help="Output FASTA file name")
args = parser.parse_args()

bam_file = args.bam_file
modkit_path = args.modkit_path
ref_genome_path = args.ref_genome_path
output_fasta = args.output_fasta

# Debug: print the parsed arguments
print("BAM file:", bam_file)
print("Modkit TSV file:", modkit_path)
print("Reference genome:", ref_genome_path)
print("Output FASTA:", output_fasta)

# Read and process the modkit TSV file
modkitfile = pd.read_csv(modkit_path, sep='\t')
modkitfile['ref_kmer'] = modkitfile['ref_kmer'].str.upper()
modkitfile["fail"] = modkitfile["fail"].astype(str).str.strip().str.lower()
modkitfile = modkitfile[modkitfile["fail"] == "false"]
print("Unique values in 'fail' column after conversion:", modkitfile["fail"].unique())
print(f"Total rows after filtering: {len(modkitfile)}")

modkit_by_read = modkitfile.groupby('read_id')

reference_genome = pysam.FastaFile(ref_genome_path)

mismatch_lock = threading.Lock()

# Add global counters for excluded reads
excluded_counts = {
    "unmapped": 0,
    "secondary": 0,
    "supplementary": 0,
    "qcfail": 0,
    "low_quality": 0
}
excluded_lock = threading.Lock()

def process_read(read):
    """
    Process a single read:
    - Checks for quality and filtering criteria.
    - Uses pre-grouped modkit data for this read.
    - Changes methylated C bases (or G on reverse strand) to 'm'.
    - Returns a tuple (read_id, processed_sequence).
    """
    # Filtering checks
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        with excluded_lock:
            if read.is_unmapped:
                excluded_counts["unmapped"] += 1
            if read.is_secondary:
                excluded_counts["secondary"] += 1
            if read.is_supplementary:
                excluded_counts["supplementary"] += 1
        return None
    if read.is_qcfail or read.mapping_quality < 50:
        with excluded_lock:
            if read.is_qcfail:
                excluded_counts["qcfail"] += 1
            if read.mapping_quality < 50:
                excluded_counts["low_quality"] += 1
        return None

    try:
        # Get the reference sequence for the read
        seq = list(read.get_reference_sequence().upper())
    except Exception as e:
        print(f"Error getting reference sequence for read {read.query_name}: {e}")
        return None

    start_read = read.reference_start
    id_of_read = read.query_name

    # Check if the read has modkit data
    if id_of_read not in modkit_by_read.groups:
        # If no modkit data, include the read without methylation modifications
        seq = ''.join(seq)
        if read.is_reverse:
            reverse_DNA = str.maketrans('ACGT', 'TGCA')
            seq = seq.translate(reverse_DNA)[::-1]
        return (id_of_read, seq)

    # Process modkit data for the read
    modkitfile_read = modkit_by_read.get_group(id_of_read)
    methylation_called = False
    for _, CpG_in_read in modkitfile_read.iterrows():
        position_m = CpG_in_read.ref_position
        relative_position = position_m - start_read
        call_code = CpG_in_read.call_code

        if call_code != 'm':
            continue

        if relative_position < 0 or relative_position >= len(seq):
            continue

        # Added quality checks from original script:
        if relative_position >= len(read.query_alignment_qualities):
            continue
        if read.query_alignment_qualities[relative_position] < 10:
            continue

        # Check and modify the sequence according to strand orientation
        if ((seq[relative_position] == 'C') and (not read.is_reverse)) or \
           ((seq[relative_position] == 'G') and (read.is_reverse)):
            seq[relative_position] = 'm'
            methylation_called = True
        else:
            # Print mismatch messages
            with mismatch_lock:
                if (seq[relative_position] == 'G') and (not read.is_reverse):
                    print(f'G mismatch in read {id_of_read}, mapping quality: {read.mapping_quality}')
                elif (seq[relative_position] == 'C') and (read.is_reverse):
                    print(f'C mismatch in read {id_of_read}, mapping quality: {read.mapping_quality}')
                else:
                    print(f"Something is wrong in read {id_of_read} at relative position {relative_position}")
                    print(f"cigar string: {read.cigarstring} D:{read.cigarstring.find('D')} I:{read.cigarstring.find('I')}")
                    segment_start = max(0, relative_position - 3)
                    segment_end = min(len(seq), relative_position + 3)
                    print(f"  seq segment: {''.join(seq[segment_start:segment_end])}")
                    print(f"  modkit data: {modkitfile_read[['ref_kmer', 'query_kmer']]}")
                    print(f"read.is_reverse: {read.is_reverse}")
                    print(f"  base at position: {seq[relative_position]}")

    # Reassemble the sequence
    seq = ''.join(seq)
    if read.is_reverse:
        reverse_DNA = str.maketrans('ACGTm', 'TGCAm')
        seq = seq.translate(reverse_DNA)[::-1]

    # Include the read even if no methylation was called
    return (id_of_read, seq)

def write_fasta(records, output_file, mode="a"):
    """Writes a list of FASTA records (tuples of (id, seq)) to the output file."""
    with open(output_file, mode) as fasta_file:
        for read_id, seq in records:
            fasta_file.write(f">{read_id}\n")
            fasta_file.write(f"{seq}\n")

def read_batches(bam, batch_size=1000):
    """Generator that yields batches of reads from the BAM file."""
    batch = []
    for read in bam:
        batch.append(read)
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch

def main():
    start_time = time.time()

    # Initialize output file (clear previous content)
    with open(output_fasta, "w") as f:
        f.write("")

    # Create a thread pool once for all batches
    pool = ThreadPool(cpu_count())

    # Process the BAM file in batches
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for batch in tqdm(read_batches(bam, batch_size=1000), desc="Processing batches"):
            results = list(tqdm(pool.imap(process_read, batch), total=len(batch), desc="Processing reads in batch"))
            fasta_records = [record for record in results if record is not None]
            write_fasta(fasta_records, output_fasta, mode="a")

    pool.close()
    pool.join()

    # Print excluded read counts
    print("\nExcluded reads summary:")
    for reason, count in excluded_counts.items():
        print(f"{reason}: {count}")

    print(f"FASTA file saved to {output_fasta}")
    print(f"Total time taken: {time.time() - start_time:.2f} seconds")

if __name__ == '__main__':
    main()
