#!/usr/bin/env bash
# ==============================================================================
# 06b_annotate_trackA.sh — Track A: Full-pipeline annotators
# ==============================================================================
# Usage:
#   bash 06b_annotate_trackA.sh --assembly polished.fasta --bakta-db /db/bakta
#
# Runs Bakta (hard requirement) and DRAM (soft) on the assembly.
# Imports RASTtk results from a manual webserver run if provided.
#
# Optional flags:
#   --rast-gto FILE          Import RASTtk GTO from BV-BRC webserver
#   --dram-db DIR            DRAM database path
#   --genus STR              Bakta genus (default: from config)
#   --species STR            Bakta species (default: from config)
#   --gram STR               Bakta gram stain (default: ?)
#   --translation-table N    Bakta translation table (default: 11)
#   --threads N              Number of threads (default: 128)
#   --sample-name STR        Sample name for output prefixes (default: MyBacteria)
#   --config FILE            Path to pipeline.conf
#   --dry-run                Print commands without executing
#   --help                   Show this help
# ==============================================================================

source "$(dirname "$0")/00_setup.sh"

# --- Defaults ----------------------------------------------------------------
threads="${threads:-128}"
assembly="${assembly:-}"
bakta_db="${bakta_db:-}"
rast_gto="${rast_gto:-}"
dram_db="${dram_db:-}"
bakta_genus="${bakta_genus:-Escherichia}"
bakta_species="${bakta_species:-coli}"
bakta_gram="${bakta_gram:-?}"
bakta_translation_table="${bakta_translation_table:-11}"
sample_name="${sample_name:-MyBacteria}"

# --- Usage -------------------------------------------------------------------
usage() {
    echo "Usage: $(basename "$0") --assembly FILE --bakta-db DIR [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --assembly FILE              Polished assembly FASTA"
    echo "  --bakta-db DIR               Bakta database path"
    echo ""
    echo "Optional:"
    echo "  --rast-gto FILE              RASTtk GTO file (from BV-BRC webserver)"
    echo "  --dram-db DIR                DRAM database path"
    echo "  --genus STR                  Bakta genus (default: Escherichia)"
    echo "  --species STR                Bakta species (default: coli)"
    echo "  --gram STR                   Bakta gram stain (default: ?)"
    echo "  --translation-table N        Translation table (default: 11)"
    echo "  --threads N                  Number of threads (default: 128)"
    echo "  --sample-name STR            Sample name prefix (default: MyBacteria)"
    echo "  --config FILE                Path to pipeline.conf"
    echo "  --dry-run                    Print commands without executing"
    echo "  --help                       Show this help"
    exit 0
}

# --- Parse arguments ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --assembly)           assembly="$2"; shift 2 ;;
        --bakta-db)           bakta_db="$2"; shift 2 ;;
        --rast-gto)           rast_gto="$2"; shift 2 ;;
        --dram-db)            dram_db="$2"; shift 2 ;;
        --genus)              bakta_genus="$2"; shift 2 ;;
        --species)            bakta_species="$2"; shift 2 ;;
        --gram)               bakta_gram="$2"; shift 2 ;;
        --translation-table)  bakta_translation_table="$2"; shift 2 ;;
        --threads)            threads="$2"; shift 2 ;;
        --sample-name)        sample_name="$2"; shift 2 ;;
        --config)             load_config "$2"; shift 2 ;;
        --dry-run)            dry_run=true; shift ;;
        --help|-h)            usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

# --- Validate ----------------------------------------------------------------
require_arg "--assembly" "${assembly}"
require_arg "--bakta-db" "${bakta_db}"
require_file "${assembly}"
require_dir "${bakta_db}"

# Hard requirement
require_tool bakta

# Soft requirements
if command -v DRAM.py &>/dev/null; then
    has_dram=true
else
    log_warn "'DRAM.py' not found — DRAM annotation will be skipped."
    has_dram=false
fi

# Optional RASTtk import
if [[ -n "${rast_gto}" ]]; then
    require_file "${rast_gto}" "Upload assembly to BV-BRC and download the GTO output"
    log_info "RASTtk GTO provided: ${rast_gto}"
else
    log_info "No RASTtk GTO provided (--rast-gto). SEED subsystem annotation skipped."
fi

# --- Sanitize locus tag (same logic as 05_annotate_assess.sh) ----------------
bakta_locus_tag=$(echo "${sample_name}" | tr -cd '[:alnum:]' | cut -c 1-12)
[[ ${#bakta_locus_tag} -lt 3 ]] && bakta_locus_tag="LOCUS"

# --- Derived paths -----------------------------------------------------------
trackA_dir="$(pwd)/14_trackA"
bakta_dir="${trackA_dir}/bakta"
rast_dir="${trackA_dir}/rasttk"
dram_dir="${trackA_dir}/dram"

run_cmd mkdir -p "${bakta_dir}"
run_cmd mkdir -p "${rast_dir}"
run_cmd mkdir -p "${dram_dir}"

# ==============================================================================
# STEP A1 — Bakta Annotation
# ==============================================================================

log_step "Step A1: Bakta annotation (${bakta_genus} ${bakta_species}, gram=${bakta_gram}, locus=${bakta_locus_tag})"

run_cmd bakta \
    --db "${bakta_db}" \
    --output "${bakta_dir}" \
    --prefix "${sample_name}" \
    --genus "${bakta_genus}" \
    --species "${bakta_species}" \
    --complete \
    --translation-table "${bakta_translation_table}" \
    --gram "${bakta_gram}" \
    --compliant \
    --locus-tag "${bakta_locus_tag}" \
    --threads "${threads}" \
    "${assembly}"

if [[ -z "${dry_run:-}" && -f "${bakta_dir}/${sample_name}.gff3" ]]; then
    bakta_cds=$(grep -c $'\tCDS\t' "${bakta_dir}/${sample_name}.gff3" 2>/dev/null || echo 0)
    log_info "Bakta annotated ${bakta_cds} CDS features."
fi

# ==============================================================================
# STEP A2 — RASTtk Import (from BV-BRC webserver)
# ==============================================================================

if [[ -n "${rast_gto}" ]]; then
    log_step "Step A2: Importing RASTtk results"

    run_cmd cp "${rast_gto}" "${rast_dir}/rast_output.gto"
    log_info "RASTtk GTO copied to ${rast_dir}/rast_output.gto"

    # If rast-export-genome is available, extract the subsystem table
    if command -v rast-export-genome &>/dev/null; then
        log_info "Extracting SEED subsystem assignments..."
        run_cmd bash -c 'rast-export-genome -i "$1" --format seed_subsystem_tab > "$2"' \
            _ "${rast_dir}/rast_output.gto" "${rast_dir}/seed_subsystems.tsv"
        run_cmd bash -c 'rast-export-genome -i "$1" --format gff > "$2"' \
            _ "${rast_dir}/rast_output.gto" "${rast_dir}/rast_annotation.gff"
    else
        log_warn "'rast-export-genome' not found. GTO file saved but SEED table not extracted."
        log_warn "Install BV-BRC CLI to extract subsystem assignments, or export from the web portal."
    fi
else
    log_info "Step A2: RASTtk skipped (no --rast-gto provided)."
fi

# ==============================================================================
# STEP A3 — DRAM Annotation
# ==============================================================================

if [[ "${has_dram}" == "true" ]]; then
    log_step "Step A3: DRAM annotation (KEGG, Pfam, dbCAN, MEROPS, VOGDB)"

    dram_args=(
        -i "${assembly}"
        -o "${dram_dir}"
        --threads "${threads}"
    )

    # DRAM may use a pre-configured database or require explicit path
    if [[ -n "${dram_db}" ]]; then
        log_info "DRAM database: ${dram_db}"
    fi

    run_cmd DRAM.py annotate "${dram_args[@]}"

    # Generate the distillation summary if annotation succeeded
    if [[ -z "${dry_run:-}" && -f "${dram_dir}/annotations.tsv" ]]; then
        run_cmd DRAM.py distill \
            -i "${dram_dir}/annotations.tsv" \
            -o "${dram_dir}/distillation"
        log_info "DRAM annotation and distillation complete."
    fi
else
    log_info "Step A3: DRAM skipped (DRAM.py not found)."
fi

# --- Summary -----------------------------------------------------------------
log_step "Track A annotation complete."
log_info "Bakta:   ${bakta_dir}/"
if [[ -n "${rast_gto}" ]]; then
    log_info "RASTtk:  ${rast_dir}/"
fi
if [[ "${has_dram}" == "true" ]]; then
    log_info "DRAM:    ${dram_dir}/"
fi
log_info ""
log_info ">>> NEXT: Also run 06c_annotate_trackB.sh and 06d_annotate_trackC.sh"
log_info "         Then run 06e_reconcile_merge.sh to merge all results."
