#!/usr/bin/env bash
# ==============================================================================
# 06d_annotate_trackC.sh — Track C: Prophage detection & gene calling
# ==============================================================================
# Usage:
#   bash 06d_annotate_trackC.sh --assembly polished.fasta --genomad-db /db/genomad
#
# Runs geNomad to detect prophage regions, extracts them, then runs Phanotate
# for phage-optimised gene prediction on those regions only.
#
# Optional flags:
#   --threads N              Number of threads (default: 128)
#   --config FILE            Path to pipeline.conf
#   --dry-run                Print commands without executing
#   --help                   Show this help
# ==============================================================================

source "$(dirname "$0")/00_setup.sh"

# --- Defaults ----------------------------------------------------------------
threads="${threads:-128}"
assembly="${assembly:-}"
genomad_db="${genomad_db:-}"

# --- Usage -------------------------------------------------------------------
usage() {
    echo "Usage: $(basename "$0") --assembly FILE --genomad-db DIR [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --assembly FILE              Polished assembly FASTA"
    echo "  --genomad-db DIR             geNomad database path"
    echo ""
    echo "Optional:"
    echo "  --threads N                  Number of threads (default: 128)"
    echo "  --config FILE                Path to pipeline.conf"
    echo "  --dry-run                    Print commands without executing"
    echo "  --help                       Show this help"
    exit 0
}

# --- Parse arguments ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --assembly)     assembly="$2"; shift 2 ;;
        --genomad-db)   genomad_db="$2"; shift 2 ;;
        --threads)      threads="$2"; shift 2 ;;
        --config)       load_config "$2"; shift 2 ;;
        --dry-run)      dry_run=true; shift ;;
        --help|-h)      usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

# --- Validate ----------------------------------------------------------------
require_arg "--assembly" "${assembly}"
require_arg "--genomad-db" "${genomad_db}"
require_file "${assembly}"
require_dir "${genomad_db}"

# Hard requirement
require_tool genomad

# Soft requirements
for tool in phanotate.py seqkit; do
    command -v "${tool}" &>/dev/null || log_warn "'${tool}' not found — will skip dependent steps."
done

# --- Derived paths -----------------------------------------------------------
trackC_dir="$(pwd)/16_trackC_prophage"
genomad_dir="${trackC_dir}/genomad"
phanotate_dir="${trackC_dir}/phanotate"

run_cmd mkdir -p "${genomad_dir}"
run_cmd mkdir -p "${phanotate_dir}"

# ==============================================================================
# STEP C1 — geNomad Prophage Detection
# ==============================================================================

log_step "Step C1: geNomad prophage detection"

run_cmd genomad end-to-end \
    "${assembly}" \
    "${genomad_dir}" \
    "${genomad_db}" \
    --threads "${threads}"

# --- Check for prophage regions -----------------------------------------------
assembly_basename=$(basename "${assembly}" .fasta)
assembly_basename=$(basename "${assembly_basename}" .fa)
assembly_basename=$(basename "${assembly_basename}" .fna)
provirus_fasta="${genomad_dir}/${assembly_basename}_find_proviruses/${assembly_basename}_provirus.fna"
provirus_summary="${genomad_dir}/${assembly_basename}_summary/${assembly_basename}_virus_summary.tsv"

prophage_count=0
if [[ -z "${dry_run:-}" ]]; then
    if [[ -f "${provirus_fasta}" ]]; then
        prophage_count=$(grep -c '^>' "${provirus_fasta}" 2>/dev/null || echo 0)
    fi
fi
log_info "geNomad detected ${prophage_count} prophage region(s)."

if [[ "${prophage_count}" -eq 0 && -z "${dry_run:-}" ]]; then
    log_info "No prophage regions found. Skipping Phanotate."
    log_step "Track C complete (no prophages detected)."
    log_info "geNomad output: ${genomad_dir}/"
    exit 0
fi

# ==============================================================================
# STEP C2 — Extract Prophage Regions
# ==============================================================================

log_step "Step C2: Extracting prophage regions for Phanotate"

# geNomad's provirus FASTA already contains the extracted prophage sequences
# with coordinates encoded in the header. We use this directly.
if [[ -z "${dry_run:-}" && -f "${provirus_fasta}" ]]; then
    run_cmd cp "${provirus_fasta}" "${phanotate_dir}/prophage_regions.fasta"
    log_info "Prophage sequences: ${phanotate_dir}/prophage_regions.fasta"
else
    run_cmd cp "${provirus_fasta}" "${phanotate_dir}/prophage_regions.fasta"
fi

# ==============================================================================
# STEP C3 — Phanotate Gene Prediction (prophage regions only)
# ==============================================================================

if command -v phanotate.py &>/dev/null; then
    log_step "Step C3: Phanotate gene prediction on prophage regions"

    run_cmd phanotate.py \
        "${phanotate_dir}/prophage_regions.fasta" \
        -o "${phanotate_dir}/phanotate_genes.gff3" \
        -f gff3

    phanotate_count=0
    if [[ -z "${dry_run:-}" && -f "${phanotate_dir}/phanotate_genes.gff3" ]]; then
        phanotate_count=$(grep -c $'\tCDS\t' "${phanotate_dir}/phanotate_genes.gff3" 2>/dev/null || echo 0)
    fi
    log_info "Phanotate predicted ${phanotate_count} phage CDS features."
else
    log_warn "Step C3: Skipping Phanotate (phanotate.py not found)."
    log_warn "Prophage regions detected by geNomad but gene prediction skipped."
fi

# --- Summary -----------------------------------------------------------------
log_step "Track C (prophage) complete."
log_info "geNomad:    ${genomad_dir}/"
log_info "Phanotate:  ${phanotate_dir}/"
if [[ -f "${provirus_summary}" ]]; then
    log_info "Summary:    ${provirus_summary}"
fi
log_info ""
log_info ">>> NEXT: Run 06e_reconcile_merge.sh to integrate prophage genes into the final annotation."
