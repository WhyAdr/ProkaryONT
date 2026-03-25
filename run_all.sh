#!/usr/bin/env bash
# ==============================================================================
# run_all.sh — Master orchestrator (config file driven)
# ==============================================================================
# Usage:
#   bash run_all.sh --config pipeline.conf
#
# Reads all variables from the config file, then calls each pipeline script
# in sequence with the appropriate flags.
# ==============================================================================

source "$(dirname "$0")/00_setup.sh"

script_dir="$(dirname "$0")"

# --- Usage -------------------------------------------------------------------
usage() {
    echo "Usage: $(basename "$0") --config FILE [--skip-curation] [--dry-run]"
    echo ""
    echo "Required:"
    echo "  --config FILE     Path to pipeline configuration file"
    echo ""
    echo "Optional:"
    echo "  --skip-curation   Skip all manual curation pauses"
    echo "  --dry-run         Print commands without executing"
    echo "  --help            Show this help"
    exit 0
}

# --- Parse arguments ---------------------------------------------------------
config_file=""
skip_curation=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)         config_file="$2"; shift 2 ;;
        --skip-curation)  skip_curation=true; shift ;;
        --dry-run)        dry_run=true; shift ;;
        --help|-h)        usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

require_arg "--config" "${config_file}"
load_config "${config_file}"

# --- Validate required config keys ---
require_arg "input_fastq (in config)" "${input_fastq:-}"
require_arg "sequencing_summary (in config)" "${sequencing_summary:-}"
require_arg "pod5_dir (in config)" "${pod5_dir:-}"
require_arg "sample_name (in config)" "${sample_name:-}"

# --- Build common flags ------------------------------------------------------
common_flags=()
[[ -n "${dry_run:-}" ]] && common_flags+=(--dry-run)

pipeline_start=$(date +%s)

log_step "========================================="
log_step "Starting full pipeline: ${sample_name}"
log_step "========================================="
log_info "Config:  ${config_file}"
log_info "Threads: ${threads:-128}"
echo ""

# ==============================================================================
# Stage 1 — QC & Estimation
# ==============================================================================

stage_start=$(date +%s)
log_step ">>> Stage 1: QC & genome size estimation"

if ! bash "${script_dir}/01_qc_estimate.sh" \
    --input-fastq "${input_fastq}" \
    --sequencing-summary "${sequencing_summary}" \
    --threads "${threads:-128}" \
    --kmer-size "${meryl_kmer_size:-21}" \
    --nanoplot-color "${nanoplot_color:-green}" \
    "${common_flags[@]+"${common_flags[@]}"}"; then
    log_error "01_qc_estimate.sh failed. Pipeline aborted."
fi

log_info "Stage 1 completed in $(( $(date +%s) - stage_start ))s"
echo ""

# ==============================================================================
# Stage 2 — Filter & Assemble
# ==============================================================================

stage_start=$(date +%s)
log_step ">>> Stage 2: Filtering & assembly"

filter_flags=(
    --input-fastq "${input_fastq}"
    --threads "${threads:-128}"
    --read-type "${read_type:-ont_r10}"
    --sample-name "${sample_name}"
    --min-length "${filtlong_min_length:-200}"
    --min-qscore "${min_qscore:-7}"
    --keep-percent "${filtlong_keep_percent:-90}"
    --subsample-count "${subsample_count:-4}"
    --parallel-jobs "${parallel_jobs:-4}"
    --canu-parallel-jobs "${canu_parallel_jobs:-2}"
)
[[ -n "${genome_size_override:-}" ]] && filter_flags+=(--genome-size "${genome_size_override}")
[[ -n "${min_read_depth:-}" ]]      && filter_flags+=(--min-read-depth "${min_read_depth}")
[[ -n "${skip_curation}" ]]          && filter_flags+=(--skip-curation)

if ! bash "${script_dir}/02_filter_assemble.sh" \
    "${filter_flags[@]}" \
    "${common_flags[@]+"${common_flags[@]}"}"; then
    log_error "02_filter_assemble.sh failed. Pipeline aborted."
fi

log_info "Stage 2 completed in $(( $(date +%s) - stage_start ))s"
echo ""

# ==============================================================================
# Stage 3 — Polish & Orient
# ==============================================================================

stage_start=$(date +%s)
log_step ">>> Stage 3: Polishing & reorientation"

if ! bash "${script_dir}/03_polish_orient.sh" \
    --assembly autocycler_consensus.fasta \
    --pod5-dir "${pod5_dir}" \
    --threads "${threads:-128}" \
    --sample-name "${sample_name}" \
    --dorado-model "${dorado_model:-sup}" \
    --min-qscore "${dorado_min_qscore:-7}" \
    "${common_flags[@]+"${common_flags[@]}"}"; then
    log_error "03_polish_orient.sh failed. Pipeline aborted."
fi

log_info "Stage 3 completed in $(( $(date +%s) - stage_start ))s"
echo ""

# ==============================================================================
# Final summary
# ==============================================================================

pipeline_end=$(date +%s)
elapsed=$(( pipeline_end - pipeline_start ))
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))

log_step "========================================="
log_step "Assembly pipeline complete!"
log_step "========================================="
log_info "Sample:         ${sample_name}"
log_info "Total runtime:  ${hours}h ${minutes}m ${seconds}s"
log_info ""
log_info "Key outputs:"
log_info "  Consensus:    autocycler_consensus.fasta"
log_info "  Polished:     polished_assembly.fasta"
log_info "  Final:        dnaapler_reoriented.fasta"
log_info "  Metrics:      metrics.tsv"
log_info "  Log:          ${log_file}"
log_info ""
log_info "Next steps (run manually):"
log_info "  1. bash scripts/04_taxonomy.sh --config scripts/pipeline.conf"
log_info "  2. Review taxonomy, update genus/species/gram in pipeline.conf"
log_info "  3. bash scripts/05_annotate_assess.sh --assembly dnaapler_reoriented.fasta --bakta-db /path/to/bakta_db"
