#!/usr/bin/env bash
# ==============================================================================
# 01_qc_estimate.sh — Read QC, genome size estimation, k-mer profiling
# ==============================================================================
# Usage:
#   bash 01_qc_estimate.sh --input-fastq reads.fq.gz --sequencing-summary summary.txt
#
# Optional flags:
#   --config FILE           Path to pipeline.conf (values override defaults)
#   --threads N             Number of threads (default: 128)
#   --memory N              Memory in GB for Meryl (default: auto from system)
#   --kmer-size N           K-mer size for Meryl (default: 21)
#   --nanoplot-color C      Color for NanoPlot graphs (default: green)
#   --dry-run               Print commands without executing
#   --help                  Show this help
# ==============================================================================

source "$(dirname "$0")/00_setup.sh"

# --- Defaults ----------------------------------------------------------------
threads="${threads:-128}"
input_fastq="${input_fastq:-}"
sequencing_summary="${sequencing_summary:-}"
meryl_memory="${meryl_memory:-}"
meryl_kmer_size="${meryl_kmer_size:-21}"
nanoplot_color="${nanoplot_color:-green}"
config_file="${config_file:-}"

# --- Usage -------------------------------------------------------------------
usage() {
    echo "Usage: $(basename "$0") --input-fastq FILE --sequencing-summary FILE [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --input-fastq FILE          Raw FASTQ file (gzipped ok)"
    echo "  --sequencing-summary FILE   Sequencing summary from MinKNOW"
    echo ""
    echo "Optional:"
    echo "  --config FILE               Path to pipeline.conf"
    echo "  --threads N                 Number of threads (default: 128)"
    echo "  --memory N                  Memory in GB for Meryl (default: auto)"
    echo "  --kmer-size N               K-mer size for Meryl (default: 21)"
    echo "  --nanoplot-color C          Color for NanoPlot graphs (default: green)"
    echo "  --dry-run                   Print commands without executing"
    echo "  --help                      Show this help"
    exit 0
}

# --- Parse arguments ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-fastq)        input_fastq="$2"; shift 2 ;;
        --sequencing-summary) sequencing_summary="$2"; shift 2 ;;
        --config)             config_file="$2"; shift 2 ;;
        --threads)            threads="$2"; shift 2 ;;
        --memory)             meryl_memory="$2"; shift 2 ;;
        --kmer-size)          meryl_kmer_size="$2"; shift 2 ;;
        --nanoplot-color)     nanoplot_color="$2"; shift 2 ;;
        --dry-run)            dry_run=true; shift ;;
        --help|-h)            usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

# --- Load config (CLI flags parsed above take priority) ----------------------
[[ -n "${config_file}" ]] && load_config "${config_file}"

# --- Validate ----------------------------------------------------------------
require_arg "--input-fastq" "${input_fastq}"
require_arg "--sequencing-summary" "${sequencing_summary}"
require_file "${input_fastq}" ""
require_file "${sequencing_summary}" ""

require_tool NanoPlot
require_tool lrge
require_tool raven
require_tool meryl

# --- Derived paths -----------------------------------------------------------
qc_dir="$(pwd)/01_qc"
genome_size_dir="$(pwd)/02_genome_size"

# ==============================================================================
# STEP 1 — NanoPlot QC
# ==============================================================================

log_step "Step 1: NanoPlot QC"
mkdir -p "${qc_dir}"

log_info "Running NanoPlot on sequencing summary..."
run_cmd NanoPlot --threads "${threads}" \
    --summary "${sequencing_summary}" \
    --loglength \
    -o "${qc_dir}/run_summary_profile"

log_info "Running NanoPlot on raw FASTQ..."
run_cmd NanoPlot --threads "${threads}" \
    -c "${nanoplot_color}" --fastq "${input_fastq}" \
    -o "${qc_dir}/NanoPlot_sample" \
    --loglength --plots dot --N50

log_info ">>> CHECK: Open ${qc_dir}/NanoPlot_sample/NanoPlot-report.html"

# ==============================================================================
# STEP 3 — Genome Size Estimation
# ==============================================================================

log_step "Step 3: Genome size estimation"
mkdir -p "${genome_size_dir}"

log_info "Running LRGE..."
run_cmd bash -c 'lrge -P ont -t "$1" -s 123 "$2" > "$3"' \
    _ "${threads}" "${input_fastq}" "${genome_size_dir}/lrge_output.txt"
if [[ -z "${dry_run:-}" && -f "${genome_size_dir}/lrge_output.txt" ]]; then
    lrge_size=$(cat "${genome_size_dir}/lrge_output.txt")
    log_info "LRGE estimated genome size: ${lrge_size}"
fi

log_info "Running Raven quick assembly (0 polish)..."
run_cmd bash -c 'raven --threads "$1" -p 0 \
    --graphical-fragment-assembly "$2" \
    "$3" > "$4"' \
    _ "${threads}" "${genome_size_dir}/assemblyGraph.gfa" "${input_fastq}" "${genome_size_dir}/assembly.fasta"

if [[ -f "${genome_size_dir}/assembly.fasta" ]]; then
    raven_size=$(grep -v '^>' "${genome_size_dir}/assembly.fasta" | tr -d '\n' | wc -c)
    log_info "Raven assembly size: ${raven_size} bp"
fi

if [[ -n "${lrge_size:-}" && -n "${raven_size:-}" ]]; then
    mean_genome_size=$(( (2 * raven_size + lrge_size) / 3 ))
    echo "${mean_genome_size}" > "${genome_size_dir}/mean_genome_size.txt"
    log_info "Weighted Mean Genome Size ((2*Raven + LRGE) / 3): ${mean_genome_size} bp"
fi

log_info "Running Meryl k-mer counting..."
# Default meryl memory to available system GB if not set
if [[ -z "${meryl_memory}" ]]; then
    meryl_memory=$(free -g 2>/dev/null | awk '/^Mem:/{print $7}')
    meryl_memory="${meryl_memory:-16}"
fi
run_cmd meryl count k="${meryl_kmer_size}" memory="${meryl_memory}" threads="${threads}" compress \
    output "${genome_size_dir}/genome.meryl" "${input_fastq}"
run_cmd bash -c 'meryl histogram "$1" > "$2"' \
    _ "${genome_size_dir}/genome.meryl" "${genome_size_dir}/genome.hist"

# --- Summary -----------------------------------------------------------------
log_step "QC & estimation complete."
log_info "LRGE estimate:  ${lrge_size:-N/A}"
log_info "Raven size:     ${raven_size:-N/A} bp"
log_info "Mean (Weight):  ${mean_genome_size:-N/A} bp"
log_info "NanoPlot:       ${qc_dir}/"
log_info "Meryl db:       ${genome_size_dir}/genome.meryl"
