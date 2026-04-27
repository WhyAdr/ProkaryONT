#!/usr/bin/env bash
# ==============================================================================
# 06e_reconcile_merge.sh — Reconcile & merge all annotation evidence
# ==============================================================================
# Usage:
#   bash 06e_reconcile_merge.sh --assembly polished.fasta \
#       --consensus-gff consensus_genes.gff3 --bakta-dir 14_trackA/bakta
#
# Merges gene predictions from Track B consensus with functional annotations
# from Tracks A, B, and C into unified output: GFF3 + TSV + GenBank + JSON.
#
# Optional flags:
#   --rast-gto FILE          RASTtk GTO for SEED subsystem import
#   --dram-dir DIR           DRAM output directory
#   --emapper-out DIR        eggNOG-mapper output directory
#   --interproscan-out DIR   InterProScan output directory
#   --kofamscan-out DIR      KofamScan output directory
#   --diamond-out DIR        DIAMOND/SwissProt output directory
#   --phanotate-gff FILE     Phanotate GFF from Track C
#   --overlap-threshold N    Coordinate match threshold (default: 0.80)
#   --sample-name STR        Sample name for outputs (default: MyBacteria)
#   --config FILE            Path to pipeline.conf
#   --dry-run                Print commands without executing
#   --help                   Show this help
# ==============================================================================

source "$(dirname "$0")/00_setup.sh"

# --- Defaults ----------------------------------------------------------------
assembly="${assembly:-}"
consensus_gff="${consensus_gff:-}"
bakta_dir="${bakta_dir:-}"
rast_gto="${rast_gto:-}"
dram_dir="${dram_dir:-}"
emapper_out="${emapper_out:-}"
interproscan_out="${interproscan_out:-}"
kofamscan_out="${kofamscan_out:-}"
diamond_out="${diamond_out:-}"
phanotate_gff="${phanotate_gff:-}"
overlap_threshold="${overlap_threshold:-0.80}"
sample_name="${sample_name:-MyBacteria}"

# --- Usage -------------------------------------------------------------------
usage() {
    echo "Usage: $(basename "$0") --assembly FILE --consensus-gff FILE --bakta-dir DIR [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --assembly FILE              Assembly FASTA"
    echo "  --consensus-gff FILE         Consensus gene GFF3 from 06a"
    echo "  --bakta-dir DIR              Bakta output directory from 06b"
    echo ""
    echo "Optional:"
    echo "  --rast-gto FILE              RASTtk GTO for SEED subsystem import"
    echo "  --dram-dir DIR               DRAM output directory from 06b"
    echo "  --emapper-out DIR            eggNOG-mapper output directory from 06c"
    echo "  --interproscan-out DIR       InterProScan output directory from 06c"
    echo "  --kofamscan-out DIR          KofamScan output directory from 06c"
    echo "  --diamond-out DIR            DIAMOND/SwissProt output directory from 06c"
    echo "  --phanotate-gff FILE         Phanotate GFF3 from 06d"
    echo "  --overlap-threshold N        Coordinate match threshold (default: 0.80)"
    echo "  --sample-name STR            Sample name for outputs (default: MyBacteria)"
    echo "  --config FILE                Path to pipeline.conf"
    echo "  --dry-run                    Print commands without executing"
    echo "  --help                       Show this help"
    exit 0
}

# --- Parse arguments ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --assembly)           assembly="$2"; shift 2 ;;
        --consensus-gff)      consensus_gff="$2"; shift 2 ;;
        --bakta-dir)          bakta_dir="$2"; shift 2 ;;
        --rast-gto)           rast_gto="$2"; shift 2 ;;
        --dram-dir)           dram_dir="$2"; shift 2 ;;
        --emapper-out)        emapper_out="$2"; shift 2 ;;
        --interproscan-out)   interproscan_out="$2"; shift 2 ;;
        --kofamscan-out)      kofamscan_out="$2"; shift 2 ;;
        --diamond-out)        diamond_out="$2"; shift 2 ;;
        --phanotate-gff)      phanotate_gff="$2"; shift 2 ;;
        --overlap-threshold)  overlap_threshold="$2"; shift 2 ;;
        --sample-name)        sample_name="$2"; shift 2 ;;
        --config)             load_config "$2"; shift 2 ;;
        --dry-run)            dry_run=true; shift ;;
        --help|-h)            usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

# --- Validate ----------------------------------------------------------------
require_arg "--assembly" "${assembly}"
require_arg "--consensus-gff" "${consensus_gff}"
require_arg "--bakta-dir" "${bakta_dir}"
require_file "${assembly}"
require_file "${consensus_gff}" "06a_predict_genes.sh"
require_dir "${bakta_dir}" "06b_annotate_trackA.sh"

require_tool python3

# Check for the reconciliation script
script_dir="$(cd "$(dirname "$0")" && pwd)"
reconcile_script="${script_dir}/utils/reconcile_annotations.py"
require_file "${reconcile_script}" "Ensure scripts/utils/reconcile_annotations.py exists"

# --- Derived paths -----------------------------------------------------------
output_dir="$(pwd)/17_final_annotation"
run_cmd mkdir -p "${output_dir}"

# ==============================================================================
# STEP E1 — Reconcile & Merge All Evidence
# ==============================================================================

log_step "Step E1: Reconciling and merging all annotation evidence"

# Build argument list for the reconciliation script
reconcile_args=(
    --assembly "${assembly}"
    --consensus-gff "${consensus_gff}"
    --bakta-dir "${bakta_dir}"
    --overlap-threshold "${overlap_threshold}"
    --sample-name "${sample_name}"
    --output-dir "${output_dir}"
    --log-file "${output_dir}/reconciliation.log"
)

# Optional Track A evidence
if [[ -n "${rast_gto}" && -f "${rast_gto}" ]]; then
    reconcile_args+=(--rast-gto "${rast_gto}")
    log_info "Including RASTtk SEED subsystem annotations."
fi

if [[ -n "${dram_dir}" && -d "${dram_dir}" ]]; then
    reconcile_args+=(--dram-dir "${dram_dir}")
    log_info "Including DRAM annotations (MEROPS, dbCAN, VOGDB)."
fi

# Optional Track B evidence
if [[ -n "${emapper_out}" && -d "${emapper_out}" ]]; then
    reconcile_args+=(--emapper-dir "${emapper_out}")
    log_info "Including eggNOG-mapper annotations."
fi

if [[ -n "${interproscan_out}" && -d "${interproscan_out}" ]]; then
    reconcile_args+=(--interproscan-dir "${interproscan_out}")
    log_info "Including InterProScan domain annotations."
fi

if [[ -n "${kofamscan_out}" && -d "${kofamscan_out}" ]]; then
    reconcile_args+=(--kofamscan-dir "${kofamscan_out}")
    log_info "Including KofamScan KEGG KO assignments."
fi

if [[ -n "${diamond_out}" && -d "${diamond_out}" ]]; then
    reconcile_args+=(--diamond-dir "${diamond_out}")
    log_info "Including DIAMOND/SwissProt hits."
fi

# Optional Track C evidence
if [[ -n "${phanotate_gff}" && -f "${phanotate_gff}" ]]; then
    reconcile_args+=(--phanotate-gff "${phanotate_gff}")
    log_info "Including Phanotate prophage gene predictions."
fi

run_cmd python3 "${reconcile_script}" "${reconcile_args[@]}"

# ==============================================================================
# STEP E2 — Report Results
# ==============================================================================

if [[ -z "${dry_run:-}" ]]; then
    log_step "Step E2: Reporting results"

    if [[ -f "${output_dir}/final_annotation.gff3" ]]; then
        final_cds=$(grep -c $'\tCDS\t' "${output_dir}/final_annotation.gff3" 2>/dev/null || echo 0)
        log_info "Final annotation: ${final_cds} CDS features."
    fi

    if [[ -f "${output_dir}/annotation_matrix.tsv" ]]; then
        tsv_cols=$(head -n1 "${output_dir}/annotation_matrix.tsv" | awk -F'\t' '{print NF}' 2>/dev/null || echo 0)
        tsv_rows=$(tail -n +2 "${output_dir}/annotation_matrix.tsv" | wc -l 2>/dev/null || echo 0)
        log_info "Annotation matrix: ${tsv_rows} genes × ${tsv_cols} columns."
    fi

    if [[ -f "${output_dir}/summary.json" ]]; then
        log_info "Per-isolate summary: ${output_dir}/summary.json"
    fi
fi

# --- Summary -----------------------------------------------------------------
log_step "Ensemble annotation pipeline complete!"
log_info ""
log_info "=== OUTPUT FILES ==="
log_info "GFF3:       ${output_dir}/final_annotation.gff3"
log_info "TSV:        ${output_dir}/annotation_matrix.tsv"
log_info "GenBank:    ${output_dir}/final_annotation.gbk"
log_info "JSON:       ${output_dir}/summary.json"
log_info "Merge log:  ${output_dir}/reconciliation.log"
log_info ""
log_info "=== EVIDENCE SOURCES USED ==="
log_info "Consensus genes:  ${consensus_gff}"
log_info "Bakta:            ${bakta_dir}/"
[[ -n "${rast_gto}" ]]        && log_info "RASTtk:           ${rast_gto}"
[[ -n "${dram_dir}" ]]        && log_info "DRAM:             ${dram_dir}/"
[[ -n "${emapper_out}" ]]     && log_info "eggNOG-mapper:    ${emapper_out}/"
[[ -n "${interproscan_out}" ]]&& log_info "InterProScan:     ${interproscan_out}/"
[[ -n "${kofamscan_out}" ]]   && log_info "KofamScan:        ${kofamscan_out}/"
[[ -n "${diamond_out}" ]]     && log_info "DIAMOND/SwissProt:${diamond_out}/"
[[ -n "${phanotate_gff}" ]]   && log_info "Phanotate:        ${phanotate_gff}"
