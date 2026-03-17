#!/usr/bin/env bash
# ==============================================================================
# 03_polish_orient.sh — Dorado polishing & Dnaapler reorientation
# ==============================================================================
# Usage:
#   bash 03_polish_orient.sh --assembly consensus.fasta --pod5-dir /data/pod5/
#
# Optional flags:
#   --config FILE           Path to pipeline.conf (values override defaults)
#   --threads N             Number of threads (default: 128)
#   --sample-name NAME      Sample name for dnaapler prefix (default: MyBacteria)
#   --dorado-model MODEL    Dorado basecaller model (default: sup)
#   --min-qscore N          Minimum quality score filter (default: 7)
#   --cleanup-bam           Remove intermediate BAM files after polishing
#   --dry-run               Print commands without executing
#   --help                  Show this help
# ==============================================================================

source "$(dirname "$0")/00_setup.sh"

# --- Defaults ----------------------------------------------------------------
threads="${threads:-128}"
assembly="${assembly:-}"
pod5_dir="${pod5_dir:-}"
sample_name="${sample_name:-MyBacteria}"
dorado_model="${dorado_model:-sup}"
dorado_min_qscore="${dorado_min_qscore:-7}"
cleanup_bam="${cleanup_bam:-}"
config_file="${config_file:-}"

# --- Usage -------------------------------------------------------------------
usage() {
    echo "Usage: $(basename "$0") --assembly FILE --pod5-dir DIR [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --assembly FILE             Draft assembly FASTA (e.g. autocycler_consensus.fasta)"
    echo "  --pod5-dir DIR              Directory containing POD5 files"
    echo ""
    echo "Optional:"
    echo "  --config FILE               Path to pipeline.conf"
    echo "  --threads N                 Number of threads (default: 128)"
    echo "  --sample-name NAME          Sample name for dnaapler prefix (default: MyBacteria)"
    echo "  --dorado-model MODEL        Basecaller model: fast|hac|sup (default: sup)"
    echo "  --min-qscore N              Min quality score for basecalling (default: 7)"
    echo "  --cleanup-bam               Remove intermediate BAM files after polishing"
    echo "  --dry-run                   Print commands without executing"
    echo "  --help                      Show this help"
    exit 0
}

# --- Parse arguments ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --assembly)      assembly="$2"; shift 2 ;;
        --pod5-dir)      pod5_dir="$2"; shift 2 ;;
        --config)        config_file="$2"; shift 2 ;;
        --threads)       threads="$2"; shift 2 ;;
        --sample-name)   sample_name="$2"; shift 2 ;;
        --dorado-model)  dorado_model="$2"; shift 2 ;;
        --min-qscore)    dorado_min_qscore="$2"; shift 2 ;;
        --cleanup-bam)   cleanup_bam=true; shift ;;
        --dry-run)       dry_run=true; shift ;;
        --help|-h)       usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

# --- Load config (CLI flags parsed above take priority) ----------------------
[[ -n "${config_file}" ]] && load_config "${config_file}"

# --- Validate ----------------------------------------------------------------
require_arg "--assembly" "${assembly}"
require_arg "--pod5-dir" "${pod5_dir}"
require_file "${assembly}" "02_filter_assemble.sh"
require_dir "${pod5_dir}" ""

require_tool dorado
require_tool samtools
require_tool dnaapler

# --- Derived paths -----------------------------------------------------------
polished_assembly="$(pwd)/polished_assembly.fasta"
reoriented_assembly="$(pwd)/dnaapler_reoriented.fasta"

# ==============================================================================
# STEP 5 — Dorado Polishing
# ==============================================================================

log_step "Step 5: Dorado polishing (model=${dorado_model}, min-qscore=${dorado_min_qscore})"

log_info "5a. Re-basecalling with move tables..."
run_cmd bash -c 'dorado basecaller "$1" "$2" \
    --emit-moves --min-qscore "$3" \
    > all_reads_w_moves.bam' \
    _ "${dorado_model}" "${pod5_dir}" "${dorado_min_qscore}"

# --- Validate Dorado basecalling output ---
if [[ -z "${dry_run:-}" ]]; then
    log_info "Validating Dorado basecalling output..."
    if ! samtools quickcheck all_reads_w_moves.bam; then
        log_error "all_reads_w_moves.bam is corrupt or truncated. Dorado may have crashed."
    fi
    read_count=$(samtools view -c all_reads_w_moves.bam)
    if [[ "${read_count}" -eq 0 ]]; then
        log_error "Dorado produced 0 reads. Verify POD5 input path and model compatibility."
    fi
    log_info "Dorado basecalling validated: ${read_count} reads."
fi

log_info "5b. Aligning to draft assembly..."
run_cmd bash -c 'dorado aligner "$1" "$2" --threads "$3" | samtools sort -@ "$3" -o aligned.sorted.bam' \
    _ "${assembly}" all_reads_w_moves.bam "${threads}"
run_cmd samtools index -@ "${threads}" aligned.sorted.bam

log_info "5c. Checking for move tables & polishing..."
if [[ -z "${dry_run:-}" ]]; then
    # Temporarily disable pipefail: grep -m 1 -q exits early, causing
    # samtools view to receive SIGPIPE (exit 141) under pipefail.
    set +o pipefail
    if ! samtools view aligned.sorted.bam | grep -m 1 -q "mv:B:c"; then
        log_warn "WARNING: No move tables (mv:B:c tags) detected in BAM."
        log_warn "         Dorado polishing requires move tables. Check if your dorado_model supports them."
        log_warn "         Polishing may fail or be ineffective."
    fi
    set -o pipefail
fi

run_cmd bash -c 'dorado polish "$1" "$2" \
    --bacteria --threads "$3" --infer-threads "$3" \
    > "$4"' \
    _ aligned.sorted.bam "${assembly}" "${threads}" "${polished_assembly}"

log_info ">>> CHECK: Compare ${polished_assembly} to ${assembly}"

# --- Optional BAM cleanup ---
if [[ -n "${cleanup_bam}" && -z "${dry_run:-}" ]]; then
    log_info "Cleaning up intermediate BAM files..."
    run_cmd rm -f all_reads_w_moves.bam aligned.sorted.bam aligned.sorted.bam.bai
fi

# ==============================================================================
# STEP 6 — Dnaapler Reorientation
# ==============================================================================

log_step "Step 6: Dnaapler reorientation"

log_info "6a. Reorienting contigs..."
run_cmd dnaapler all \
    -i "${polished_assembly}" \
    -o dnaapler_out \
    -p "${sample_name}" \
    -t "${threads}"

run_cmd ln -sf "$(pwd)/dnaapler_out/${sample_name}_reoriented.fasta" "${reoriented_assembly}"

log_info ">>> CHECK: Review dnaapler_out/ for start gene identification."

# --- Summary -----------------------------------------------------------------
log_step "Polishing & reorientation complete."
log_info "Polished:    ${polished_assembly}"
log_info "Reoriented:  ${reoriented_assembly}"
