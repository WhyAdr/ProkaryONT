#!/usr/bin/env python3
# ==============================================================================
# generate_mock_environment.py — Setup/cleanup mock binaries and sequences for testing
# ==============================================================================
# Usage:
#   python utils/generate_mock_environment.py          # Setup mock environment
#   python utils/generate_mock_environment.py --clean  # Remove mock environment
# ==============================================================================

import os
import sys
import stat
import gzip
import shutil

TOOLS = [
    "NanoPlot", "lrge", "raven", "meryl", "fastcat", "porechop_abi",
    "snikt.R", "filtlong", "seqkit", "fastplong", "autocycler",
    "parallel", "dorado", "samtools", "dnaapler", "necat",
    "canu", "flye", "metaMDBG", "miniasm", "minipolish", "minimap2",
    "plassembler", "hifiasm", "myloasm", "nextdenovo", "wtdbg2",
    "rasusa", "lja", "Ilesta", "racon",
    "nucmer", "mummerplot"
]

FILES = [
    "sequencing_summary.txt",
    "polished_assembly.fasta"
]

DIRS = [
    "pod5_dir",
    "autocycler_out",
    "plassembler_db",
]

def setup_mock():
    print("Setting up mock testing environment...")
    
    # 1. Create mock_bin directory and executables
    os.makedirs("mock_bin", exist_ok=True)
    for tool in TOOLS:
        path = os.path.join("mock_bin", tool)
        if tool == "seqkit":
            with open(path, "w", newline="\n") as f:
                f.write("""#!/bin/sh
if [ "$1" = "stats" ]; then
    printf "file\\tformat\\ttype\\tnum_seqs\\tsum_len\\tmin_len\\tavg_len\\tmax_len\\tQ20\\tQ30\\tN50\\tG/C\\n"
    case "$*" in
        *filtered_input.fastq.gz*)
            printf "other.fastq.gz\\tFASTQ\\tDNA\\t50\\t800000000\\t200\\t500\\t1000\\t0.0\\t0.0\\t500\\t0.0\\n"
            ;;
        *input.fastq.gz*)
            printf "input.fastq.gz\\tFASTQ\\tDNA\\t100\\t1000000000\\t200\\t500\\t1000\\t0.0\\t0.0\\t500\\t0.0\\n"
            ;;
        *)
            printf "other.fastq.gz\\tFASTQ\\tDNA\\t50\\t800000000\\t200\\t500\\t1000\\t0.0\\t0.0\\t500\\t0.0\\n"
            ;;
    esac
    exit 0
fi
exit 0
""")
        elif tool == "fastcat":
            with open(path, "w", newline="\n") as f:
                f.write("""#!/bin/sh
if [ "$1" = "lint" ]; then
    printf '@mock_zero_sdust\\nACTG\\n+\\nIIII\\n'
    printf 'Read mock_sdust_positive masked fraction 0.99 exceeds threshold 0.00, skipping.\\n' >&2
fi
exit 0
""")
        else:
            with open(path, "w", newline="\n") as f:
                f.write("#!/bin/sh\nexit 0\n")
        try:
            st = os.stat(path)
            os.chmod(path, st.st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
        except Exception:
            pass
    print("  Created mock binaries in mock_bin/")

    # 2. Create directories
    for d in DIRS:
        os.makedirs(d, exist_ok=True)
    print("  Created directories: " + ", ".join(DIRS))

    # 3. Touch files
    for file in FILES:
        with open(file, "w") as f:
            pass
    print("  Created blank files: " + ", ".join(FILES))

    # 4. Write dummy gzip FASTQ records
    fastq_data = b"@read1\nACTG\n+\nIIII\n@read2\nCCCC\n+\nIIII\n"
    with gzip.open("input.fastq.gz", "wb") as f:
        f.write(fastq_data)
    with gzip.open("filtered_input.fastq.gz", "wb") as f:
        f.write(fastq_data)
    print("  Created gzipped FASTQ inputs: input.fastq.gz, filtered_input.fastq.gz")

    # 5. Create autocycler_out/consensus_assembly.fasta
    with open(os.path.join("autocycler_out", "consensus_assembly.fasta"), "w") as f:
        pass
    print("  Created mock autocycler consensus assembly.")
    print("\nMock environment setup successfully!")
    print("Run dry-runs using: export PATH=\"$(pwd)/mock_bin:$PATH\"")

def clean_mock():
    print("Cleaning up mock environment...")
    # Delete directories
    for d in ["mock_bin", "pod5_dir", "autocycler_out", "plassembler_db", "subsampled_reads", "01_qc", "02_genome_size", "assemblies"]:
        if os.path.exists(d):
            shutil.rmtree(d, ignore_errors=True)
    # Delete files
    to_delete = [
        "input.fastq.gz", "filtered_input.fastq.gz", "sequencing_summary.txt",
        "polished_assembly.fasta", "pipeline.log", "dnaapler_reoriented.fasta",
        ".success_assembly", ".success_cluster", ".success_trim", "contig_depths.tsv",
        "metrics.tsv", "contig_characteristics.tsv", "plassembler_summary.tsv",
        "rasusa_subsample.tsv", "subsample.yaml", "no_keep.log", "keep90.log"
    ]
    for file in to_delete:
        if os.path.exists(file):
            try:
                os.remove(file)
            except Exception:
                pass
    print("Cleaned up successfully.")

def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--clean":
        clean_mock()
    else:
        setup_mock()

if __name__ == "__main__":
    main()
