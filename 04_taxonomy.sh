#!/usr/bin/env bash
# ==============================================================================
# 04_taxonomy.sh — Taxonomic classification, MLST typing & 16S extraction
# ==============================================================================
# Usage:
#   bash 04_taxonomy.sh --assembly reoriented.fasta --gtdbtk-db /db/gtdbtk
#
# Optional flags:
#   --config FILE            Path to pipeline.conf (values override defaults)
#   --threads N              Number of threads (default: 128)
#   --barrnap-kingdom K      Barrnap kingdom (bac, arc, mito, euk) (default: bac)
#   --dry-run                Print commands without executing
#   --help                   Show this help
# ==============================================================================

source "$(dirname "$0")/00_setup.sh"

# --- Defaults ----------------------------------------------------------------
threads="${threads:-128}"
assembly="${assembly:-}"
gtdbtk_data_path="${gtdbtk_data_path:-}"
barrnap_kingdom="${barrnap_kingdom:-bac}"
config_file="${config_file:-}"

# --- Usage -------------------------------------------------------------------
usage() {
    echo "Usage: $(basename "$0") --assembly FILE --gtdbtk-db DIR [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --assembly FILE              Final assembly FASTA (e.g. dnaapler_reoriented.fasta)"
    echo "  --gtdbtk-db DIR              GTDB-Tk database path"
    echo ""
    echo "Optional:"
    echo "  --config FILE                Path to pipeline.conf"
    echo "  --threads N                  Number of threads (default: 128)"
    echo "  --barrnap-kingdom K          Barrnap kingdom: bac, arc, mito, euk (default: bac)"
    echo "  --dry-run                    Print commands without executing"
    echo "  --help                       Show this help"
    exit 0
}

# --- Parse arguments ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --assembly)          assembly="$2"; shift 2 ;;
        --gtdbtk-db)         gtdbtk_data_path="$2"; shift 2 ;;
        --config)            config_file="$2"; shift 2 ;;
        --threads)           threads="$2"; shift 2 ;;
        --barrnap-kingdom)   barrnap_kingdom="$2"; shift 2 ;;
        --dry-run)           dry_run=true; shift ;;
        --help|-h)       usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

# --- Load config (CLI flags parsed above take priority) ----------------------
[[ -n "${config_file}" ]] && load_config "${config_file}"

# --- Validate ----------------------------------------------------------------
require_arg "--assembly" "${assembly}"
require_arg "--gtdbtk-db" "${gtdbtk_data_path}"
require_file "${assembly}" "03_polish_orient.sh"

require_tool gtdbtk

for tool in mlst barrnap; do
    command -v "${tool}" &>/dev/null || log_warn "'${tool}' not found — will skip."
done

# --- Derived paths -----------------------------------------------------------
taxonomy_dir="$(pwd)/07_taxonomy"

# ==============================================================================
# STEP 7a — GTDB-Tk Classification
# ==============================================================================

log_step "Step 7a: GTDB-Tk classification"

[[ -n "${gtdbtk_data_path}" ]] && export GTDBTK_DATA_PATH="${gtdbtk_data_path}"

run_cmd mkdir -p "${taxonomy_dir}/genomes"
run_cmd cp "${assembly}" "${taxonomy_dir}/genomes/"

run_cmd gtdbtk classify_wf \
    --genome_dir "${taxonomy_dir}/genomes" \
    --extension fasta \
    --out_dir "${taxonomy_dir}" \
    --cpus "${threads}"

log_info ">>> CHECK: Review ${taxonomy_dir}/ for taxonomic assignment."

if [[ -z "${dry_run:-}" ]] && ls "${taxonomy_dir}"/gtdbtk.*.summary.tsv 1> /dev/null 2>&1; then
    gtdbtk_class=$(awk -F'\t' 'NR==2{print $2}' "${taxonomy_dir}"/gtdbtk.*.summary.tsv 2>/dev/null)
    if [[ -n "${gtdbtk_class}" ]]; then
        log_info "    Classification: ${gtdbtk_class}"
    fi
fi

log_info "    Update genus/species/gram in pipeline.conf before running 05_annotate_assess.sh."

# ==============================================================================
# STEP 7b — MLST Typing
# ==============================================================================

if command -v mlst &>/dev/null; then
    log_step "Step 7b: MLST typing"

    run_cmd mlst "${assembly}" \
        | tee "${taxonomy_dir}/mlst_result.tsv"

    log_info "MLST result:"
    log_info "  $(cat "${taxonomy_dir}/mlst_result.tsv")"
else
    log_warn "Skipping MLST (mlst not found)."
fi

# ==============================================================================
# STEP 7c — Barrnap 16S rRNA Extraction
# ==============================================================================

if command -v barrnap &>/dev/null; then
    log_step "Step 7c: Barrnap 16S rRNA extraction"

    log_info "Predicting rRNA genes..."
    run_cmd barrnap --threads "${threads}" --kingdom "${barrnap_kingdom}" \
        "${assembly}" > "${taxonomy_dir}/rrna_predictions.gff3"

    # Extract 16S sequences from the assembly using the GFF coordinates
    if grep -q "16S" "${taxonomy_dir}/rrna_predictions.gff3" 2>/dev/null; then
        run_cmd barrnap --threads "${threads}" --kingdom "${barrnap_kingdom}" \
            --outseq "${taxonomy_dir}/16S_sequences.fasta" \
            "${assembly}" > /dev/null

        # Keep only 16S sequences
        awk '/^>.*16S_rRNA/{p=1} /^>/ && !/16S_rRNA/{p=0} p' "${taxonomy_dir}/16S_sequences.fasta" \
            > "${taxonomy_dir}/16S_only.fasta"

        count_16s=$(grep -c '^>' "${taxonomy_dir}/16S_only.fasta" 2>/dev/null || echo 0)
        log_info "Extracted ${count_16s} 16S rRNA sequence(s) → ${taxonomy_dir}/16S_only.fasta"
        log_info ">>> TIP: BLAST these at https://blast.ncbi.nlm.nih.gov/ for independent taxonomy check."
    else
        log_warn "No 16S rRNA genes found in assembly."
    fi
else
    log_warn "Skipping Barrnap (barrnap not found)."
fi

# --- Summary -----------------------------------------------------------------
log_step "Taxonomy & typing complete."
log_info "GTDB-Tk:   ${taxonomy_dir}/"
log_info "MLST:      ${taxonomy_dir}/mlst_result.tsv"
log_info "16S rRNA:  ${taxonomy_dir}/16S_only.fasta"
log_info ""
log_info ">>> NEXT: Review results, update genus/species/gram in config, then run 05_annotate_assess.sh"
