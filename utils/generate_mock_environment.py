#!/usr/bin/env python3
# ==============================================================================
# generate_mock_environment.py — Setup/cleanup mock binaries and sequences for testing
# ==============================================================================
# Usage:
#   python utils/generate_mock_environment.py          # Setup mock environment
#   python utils/generate_mock_environment.py --clean  # Remove mock environment
# ==============================================================================

import gzip
import os
import shutil
import stat
import sys

TOOLS = [
    "NanoPlot",
    "lrge",
    "raven",
    "meryl",
    "fastcat",
    "porechop_abi",
    "snikt.R",
    "chopper",
    "filtlong",
    "seqkit",
    "fastplong",
    "autocycler",
    "parallel",
    "dorado",
    "samtools",
    "dnaapler",
    "necat",
    "canu",
    "flye",
    "metaMDBG",
    "miniasm",
    "minipolish",
    "minimap2",
    "plassembler",
    "hifiasm",
    "myloasm",
    "nextdenovo",
    "wtdbg2",
    "rasusa",
    "lja",
    "Ilesta",
    "racon",
    "nucmer",
    "mummerplot",
]

FILES = ["sequencing_summary.txt", "polished_assembly.fasta"]

DIRS = [
    "pod5_dir",
    "autocycler_out",
    "plassembler_db",
    "mock_contig_policy",
]


MOCK_TOOL_SCRIPT = r"""#!/bin/sh
tool=$(basename "$0")

if [ "${1:-}" = "--version" ] || [ "${1:-}" = "-V" ]; then
    printf 'mock-%s 1.0.1\n' "$tool"
    exit 0
fi

case "$tool" in
    NanoPlot)
        out=''
        while [ "$#" -gt 0 ]; do
            case "$1" in
                -o|--outdir) out=$2; shift 2 ;;
                *) shift ;;
            esac
        done
        [ -n "$out" ] || exit 2
        mkdir -p "$out"
        printf '<html>mock NanoPlot</html>\n' > "$out/NanoPlot-report.html"
        ;;
    fastcat)
        subcommand=${1:-}
        shift || true
        case "$subcommand" in
            fastq)
                if [ "${1:-}" = "--help" ]; then
                    printf '%s\n' 'fastcat fastq --force-error --output DIRECTORY'
                    exit 0
                fi
                out=''
                input=''
                while [ "$#" -gt 0 ]; do
                    case "$1" in
                        --force-error) shift ;;
                        --output|-o) out=$2; shift 2 ;;
                        --sample|-s) shift 2 ;;
                        --*) exit 2 ;;
                        *) input=$1; shift ;;
                    esac
                done
                [ -n "$out" ] && [ ! -e "$out" ] || exit 2
                mkdir -p "$out/histograms"
                printf 'filename\tn_seqs\tn_bases\nmock.fastq.gz\t2\t8\n' > "$out/per_file.txt"
                printf 'lower\tupper\tcount\n0\t10\t2\n' > "$out/histograms/length.hist"
                printf 'lower\tupper\tcount\n0\t20\t2\n' > "$out/histograms/quality.hist"
                printf 'filename\tbasecaller\tcount\nmock.fastq.gz\tmock\t2\n' > "$out/basecallers.txt"
                printf 'filename\trun_id\tcount\nmock.fastq.gz\tmock-run\t2\n' > "$out/runids.tsv"
                if gzip -t "$input" 2>/dev/null; then gzip -dc "$input"; else cat "$input"; fi
                ;;
            lint)
                if [ "${1:-}" = "--help" ]; then
                    printf '%s\n' 'fastcat lint --threshold N --window N --max-proportion P'
                    exit 0
                fi
                maximum=0.95
                input=''
                while [ "$#" -gt 0 ]; do
                    case "$1" in
                        --threshold|--window) shift 2 ;;
                        --max-proportion) maximum=$2; shift 2 ;;
                        --*) exit 2 ;;
                        *) input=$1; shift ;;
                    esac
                done
                printf '@mock_zero_sdust\nACTG\n+\nIIII\n'
                awk -v m="$maximum" 'BEGIN { exit (0.99 > m) ? 0 : 1 }' && \
                    printf 'Read mock_sdust_positive masked fraction 0.99 exceeds threshold %s, skipping.\n' "$maximum" >&2
                ;;
            *) exit 2 ;;
        esac
        ;;
    seqkit)
        subcommand=${1:-}
        shift || true
        case "$subcommand" in
            stats)
                printf 'file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tQ1\tQ2\tQ3\tN50\tN50_num\tQ20\tQ30\tAvgQual\tGC\n'
                case "$*" in
                    *filtered_input.fastq.gz*|*.tmp.*)
                        printf 'filtered.fastq.gz\tFASTQ\tDNA\t50\t800000000\t200\t16000000\t30000000\t1000\t5000\t10000\t12000\t20\t90.0\t80.0\t12.5\t50.0\n'
                        ;;
                    *input.fastq.gz*)
                        printf 'input.fastq.gz\tFASTQ\tDNA\t100\t1000000000\t100\t10000000\t30000000\t900\t4500\t9000\t11000\t25\t85.0\t75.0\t11.0\t50.0\n'
                        ;;
                    *)
                        printf 'other.fastq\tFASTQ\tDNA\t75\t900000000\t200\t12000000\t30000000\t950\t4750\t9500\t11500\t22\t88.0\t78.0\t12.0\t50.0\n'
                        ;;
                esac
                ;;
            seq)
                input=''
                output=''
                while [ "$#" -gt 0 ]; do
                    case "$1" in
                        -o|--out-file) output=$2; shift 2 ;;
                        --min-qual|--min-len|-m) shift 2 ;;
                        -*) shift ;;
                        *) input=$1; shift ;;
                    esac
                done
                [ -n "$input" ] && [ -n "$output" ] || exit 2
                if gzip -t "$input" 2>/dev/null; then gzip -dc "$input" > "$output"; else cat "$input" > "$output"; fi
                ;;
            *) exit 0 ;;
        esac
        ;;
    lrge)
        printf '5000000\n'
        ;;
    raven)
        graph=''
        while [ "$#" -gt 0 ]; do
            case "$1" in
                --graphical-fragment-assembly) graph=$2; shift 2 ;;
                --threads|-p) shift 2 ;;
                *) shift ;;
            esac
        done
        [ -z "$graph" ] || printf 'H\tVN:Z:1.0\n' > "$graph"
        printf '>mock_raven\nACGTACGT\n'
        ;;
    meryl)
        if [ "${1:-}" = "count" ]; then
            out=''
            previous=''
            for argument in "$@"; do
                [ "$previous" = "output" ] && out=$argument
                previous=$argument
            done
            [ -n "$out" ] || exit 2
            mkdir -p "$out"
            printf 'mock meryl db\n' > "$out/merylIndex"
        elif [ "${1:-}" = "histogram" ]; then
            printf '1\t10\n2\t5\n'
        fi
        ;;
    porechop_abi)
        input=''
        output=''
        scan=false
        while [ "$#" -gt 0 ]; do
            case "$1" in
                -i) input=$2; shift 2 ;;
                -o) output=$2; shift 2 ;;
                --guess_adapter_only) scan=true; shift ;;
                --threads) shift 2 ;;
                *) shift ;;
            esac
        done
        if [ "$scan" = true ]; then
            printf 'mock adapter diagnostic\n'
        elif [ -n "$output" ]; then
            if gzip -t "$input" 2>/dev/null; then gzip -dc "$input" > "$output"; else cat "$input" > "$output"; fi
        fi
        ;;
    snikt.R)
        printf 'mock SNIKT diagnostic\n' > snikt_report.txt
        ;;
    chopper)
        if [ "${1:-}" = "--help" ]; then
            printf '%s\n' 'chopper --trim-approach MODE --split-window N --input FILE'
            exit 0
        fi
        input=''
        while [ "$#" -gt 0 ]; do
            case "$1" in
                --input|-i) input=$2; shift 2 ;;
                --trim-approach|--headcrop|--tailcrop|--cutoff|--split-window|--minlength|--threads) shift 2 ;;
                *) shift ;;
            esac
        done
        [ -n "$input" ] || exit 2
        if gzip -t "$input" 2>/dev/null; then gzip -dc "$input"; else cat "$input"; fi
        ;;
    filtlong)
        input=''
        while [ "$#" -gt 0 ]; do
            case "$1" in
                --keep_percent|--min_length) shift 2 ;;
                --*) shift ;;
                *) input=$1; shift ;;
            esac
        done
        [ -n "$input" ] || exit 2
        if gzip -t "$input" 2>/dev/null; then gzip -dc "$input"; else cat "$input"; fi
        ;;
    fastplong)
        html=''
        json=''
        while [ "$#" -gt 0 ]; do
            case "$1" in
                --html) html=$2; shift 2 ;;
                --json) json=$2; shift 2 ;;
                -i|--thread) shift 2 ;;
                *) shift ;;
            esac
        done
        [ -z "$html" ] || printf '<html>mock Fastplong</html>\n' > "$html"
        [ -z "$json" ] || printf '{"mock": true}\n' > "$json"
        ;;
    minimap2)
        case " $* " in
            *" -x asm5 "*)
                query=''
                target=''
                for argument in "$@"; do target="$query"; query="$argument"; done
                if [ "$query" = "$target" ] && grep -q '^>contained$' "$query" 2>/dev/null \
                   && grep -q '^>target$' "$query" 2>/dev/null; then
                    printf 'contained\t100\t0\t100\t+\ttarget\t120\t0\t100\t100\t100\t60\n'
                fi
                ;;
        esac
        ;;
    *)
        exit 0
        ;;
esac
exit 0
"""


def setup_mock():
    print("Setting up mock testing environment...")

    # 1. Create mock_bin directory and executables
    os.makedirs("mock_bin", exist_ok=True)
    for tool in TOOLS:
        path = os.path.join("mock_bin", tool)
        with open(path, "w", newline="\n") as f:
            f.write(MOCK_TOOL_SCRIPT)
        try:
            st = os.stat(path)
            os.chmod(path, st.st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
        except OSError:
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

    # 6. Create deterministic post-assembly assessment fixtures.
    with open(
        os.path.join("mock_contig_policy", "normal.fasta"), "w", newline="\n"
    ) as f:
        f.write(">normal_1\nACGTACGT\n>normal_2\nAACCGGTT\n")
    fragmented_records = [("target", "A" * 120), ("contained", "C" * 100)]
    for index in range(24):
        length = 121 + index
        sequence = ("ACGT" * ((length // 4) + 1))[:length]
        sequence = sequence[:-4] + f"{index:04d}".translate(
            str.maketrans("0123456789", "ACGTACGTAC")
        )
        fragmented_records.append((f"unique_{index:02d}", sequence))
    with open(
        os.path.join("mock_contig_policy", "fragmented.fasta"), "w", newline="\n"
    ) as f:
        for identifier, sequence in fragmented_records:
            f.write(f">{identifier}\n{sequence}\n")
    with open(
        os.path.join("mock_contig_policy", "contained.paf"), "w", newline="\n"
    ) as f:
        f.write("contained\t100\t0\t100\t+\ttarget\t120\t0\t100\t100\t100\t60\n")
    print("  Created mock post-assembly contig-policy FASTA/PAF fixtures.")
    print("\nMock environment setup successfully!")
    print('Run dry-runs using: export PATH="$(pwd)/mock_bin:$PATH"')


def clean_mock():
    print("Cleaning up mock environment...")
    # Delete directories
    for d in [
        "mock_bin",
        "pod5_dir",
        "autocycler_out",
        "plassembler_db",
        "subsampled_reads",
        "01_qc",
        "02_genome_size",
        "assemblies",
        "assemblies_prepared",
        "assembly_assessment",
        "mock_contig_policy",
        "provenance",
        ".prokaryont_tmp",
        ".test-tmp",
    ]:
        if os.path.exists(d):
            shutil.rmtree(d, ignore_errors=True)
    # Delete files
    to_delete = [
        "input.fastq.gz",
        "filtered_input.fastq.gz",
        "sequencing_summary.txt",
        "polished_assembly.fasta",
        "pipeline.log",
        "dnaapler_reoriented.fasta",
        ".success_assembly",
        ".success_cluster",
        ".success_trim",
        "contig_depths.tsv",
        "metrics.tsv",
        "contig_characteristics.tsv",
        "plassembler_summary.tsv",
        "rasusa_subsample.tsv",
        "subsample.yaml",
        "no_keep.log",
        "keep90.log",
        ".prokaryont_stage1.complete.tsv",
        ".prokaryont_stage2.complete.tsv",
        "02_filtering_metrics.tsv",
    ]
    for file in to_delete:
        if os.path.exists(file):
            try:
                os.remove(file)
            except OSError:
                pass
    print("Cleaned up successfully.")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--clean":
        clean_mock()
    else:
        setup_mock()


if __name__ == "__main__":
    main()
