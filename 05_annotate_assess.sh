#!/usr/bin/env bash
# ==============================================================================
# 05_annotate_assess.sh — Annotation & assembly quality assessment
# ==============================================================================
# Usage:
#   bash 05_annotate_assess.sh --assembly reoriented.fasta --bakta-db /db/bakta
#
# Optional flags:
#   --threads N              Number of threads (default: 128)
#   --sample-name NAME       Sample name (default: MyBacteria)
#   --genus STR              Bakta genus (default: Escherichia)
#   --species STR            Bakta species (default: coli)
#   --gram STR               Bakta gram stain: -|+|? (default: ?)
#   --translation-table N    Bakta translation table (default: 11)
#   --busco-lineage STR      BUSCO lineage dataset (default: bacteria_odb10)
#   --merqury-path DIR       Path to Merqury scripts
#   --meryl-db DIR           Path to Meryl k-mer database
#   --dry-run                Print commands without executing
#   --help                   Show this help
# ==============================================================================

source "$(dirname "$0")/00_setup.sh"

# --- Defaults ----------------------------------------------------------------
threads="${threads:-128}"
assembly="${assembly:-}"
sample_name="${sample_name:-MyBacteria}"
bakta_db="${bakta_db:-}"
merqury_path="${merqury_path:-}"
meryl_db="${meryl_db:-$(pwd)/02_genome_size/genome.meryl}"
bakta_genus="${bakta_genus:-Escherichia}"
bakta_species="${bakta_species:-coli}"
bakta_gram="${bakta_gram:-?}"
bakta_translation_table="${bakta_translation_table:-11}"
busco_lineage="${busco_lineage:-bacteria_odb10}"

# --- Usage -------------------------------------------------------------------
usage() {
    echo "Usage: $(basename "$0") --assembly FILE --bakta-db DIR [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --assembly FILE              Final assembly FASTA (e.g. dnaapler_reoriented.fasta)"
    echo "  --bakta-db DIR               Bakta database path"
    echo ""
    echo "Optional:"
    echo "  --threads N                  Number of threads (default: 128)"
    echo "  --sample-name NAME           Sample name (default: MyBacteria)"
    echo "  --genus STR                  Bakta genus (default: Escherichia)"
    echo "  --species STR                Bakta species (default: coli)"
    echo "  --gram STR                   Gram stain: -|+|? (default: ?)"
    echo "  --translation-table N        Translation table (default: 11)"
    echo "  --busco-lineage STR          BUSCO lineage dataset (default: bacteria_odb10)"
    echo "  --merqury-path DIR           Path to Merqury scripts (optional)"
    echo "  --meryl-db DIR               Meryl k-mer database (default: 02_genome_size/genome.meryl)"
    echo "  --dry-run                    Print commands without executing"
    echo "  --help                       Show this help"
    exit 0
}

# --- Parse arguments ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --assembly)          assembly="$2"; shift 2 ;;
        --bakta-db)          bakta_db="$2"; shift 2 ;;
        --threads)           threads="$2"; shift 2 ;;
        --sample-name)       sample_name="$2"; shift 2 ;;
        --genus)             bakta_genus="$2"; shift 2 ;;
        --species)           bakta_species="$2"; shift 2 ;;
        --gram)              bakta_gram="$2"; shift 2 ;;
        --translation-table) bakta_translation_table="$2"; shift 2 ;;
        --busco-lineage)     busco_lineage="$2"; shift 2 ;;
        --merqury-path)      merqury_path="$2"; shift 2 ;;
        --meryl-db)          meryl_db="$2"; shift 2 ;;
        --dry-run)           dry_run=true; shift ;;
        --help|-h)           usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

# --- Validate ----------------------------------------------------------------
require_arg "--assembly" "${assembly}"
require_arg "--bakta-db" "${bakta_db}"
require_file "${assembly}" "03_polish_orient.sh"

require_tool bakta

for tool in quast checkm2 busco; do
    command -v "${tool}" &>/dev/null || log_warn "'${tool}' not found — will skip."
done

# --- Derived paths -----------------------------------------------------------
merqury_dir="$(pwd)/09_merqury"
quast_dir="$(pwd)/10_quast"
checkm2_dir="$(pwd)/11_checkm2"
busco_dir="$(pwd)/12_busco"

# --- Sanitize locus tag for Bakta (must be 3-12 alphanumeric chars) ---
bakta_locus_tag=$(echo "${sample_name}" | tr -cd '[:alnum:]' | cut -c 1-12)
[[ ${#bakta_locus_tag} -lt 3 ]] && bakta_locus_tag="LOCUS"

# ==============================================================================
# STEP 8 — Bakta Annotation
# ==============================================================================

log_step "Step 8: Bakta annotation (${bakta_genus} ${bakta_species}, gram=${bakta_gram}, locus=${bakta_locus_tag})"

run_cmd bakta \
    --db "${bakta_db}" \
    --output bakta_result \
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

log_info ">>> CHECK: Review bakta_result/ for GFF3, GenBank, and summary files."

# ==============================================================================
# STEP 9 — QUAST Assembly Statistics
# ==============================================================================

if command -v quast &>/dev/null; then
    log_step "Step 9: QUAST assembly statistics"

    mkdir -p "${quast_dir}"
    run_cmd quast "${assembly}" \
        --output-dir "${quast_dir}" \
        --threads "${threads}" \
        --min-contig 0 \
        --labels "${sample_name}"

    log_info "QUAST summary:"
    if [[ -f "${quast_dir}/report.txt" ]]; then
        while IFS= read -r line; do
            log_info "  ${line}"
        done < "${quast_dir}/report.txt"
    fi
else
    log_warn "Skipping QUAST (quast not found)."
fi

# ==============================================================================
# STEP 10 — CheckM2 Completeness & Contamination
# ==============================================================================

if command -v checkm2 &>/dev/null; then
    log_step "Step 10: CheckM2 completeness & contamination"

    mkdir -p "${checkm2_dir}/genomes"
    cp "${assembly}" "${checkm2_dir}/genomes/"

    run_cmd checkm2 predict \
        --input "${checkm2_dir}/genomes" \
        --output-directory "${checkm2_dir}" \
        --threads "${threads}" \
        --force

    if [[ -f "${checkm2_dir}/quality_report.tsv" ]]; then
        log_info "CheckM2 results:"
        while IFS=$'\t' read -r name completeness contamination rest; do
            [[ "${name}" == "Name" ]] && continue
            log_info "  ${name}: completeness=${completeness}%, contamination=${contamination}%"
        done < "${checkm2_dir}/quality_report.tsv"
    fi
else
    log_warn "Skipping CheckM2 (checkm2 not found)."
fi

# ==============================================================================
# STEP 11 — BUSCO
# ==============================================================================

if command -v busco &>/dev/null; then
    log_step "Step 11: BUSCO assessment (lineage=${busco_lineage})"

    run_cmd busco \
        -i "${assembly}" \
        -o "${busco_dir}" \
        -m genome \
        -l "${busco_lineage}" \
        -c "${threads}" \
        -f

    busco_out_name="$(basename "${busco_dir}")"
    if [[ -f "${busco_dir}/short_summary.specific.${busco_lineage}.${busco_out_name}.txt" ]]; then
        log_info "BUSCO summary:"
        grep -E "^\s+C:" "${busco_dir}/short_summary.specific.${busco_lineage}.${busco_out_name}.txt" \
            | while IFS= read -r line; do log_info "  ${line}"; done
    else
        # Try to find the summary file with a glob
        shopt -s nullglob
        for f in "${busco_dir}"/short_summary*.txt; do
            log_info "BUSCO summary:"
            grep -E "^\s+C:" "$f" | while IFS= read -r line; do log_info "  ${line}"; done
            break
        done
        shopt -u nullglob
    fi
else
    log_warn "Skipping BUSCO (busco not found)."
fi

# ==============================================================================
# STEP 12 — Merqury QV Assessment
# ==============================================================================

log_step "Step 12: Merqury assessment"

if [[ -z "${merqury_path}" ]]; then
    log_warn "No --merqury-path provided. Skipping Merqury."
elif [[ ! -d "${meryl_db}" ]]; then
    log_warn "Meryl k-mer database not found at ${meryl_db}. Skipping Merqury."
    log_warn "Run 01_qc_estimate.sh first, or specify --meryl-db."
else
    require_file "${merqury_path}/merqury.sh" ""

    mkdir -p "${merqury_dir}"
    (
        cd "${merqury_dir}"
        run_cmd "${merqury_path}/merqury.sh" \
            "${meryl_db}" \
            "${assembly}" \
            "${sample_name}"
    )

    log_info ">>> CHECK: Review ${merqury_dir}/ for QV and completeness."
fi

# --- Summary -----------------------------------------------------------------
log_step "Annotation & assessment complete."
log_info "Bakta:     bakta_result/"
log_info "QUAST:     ${quast_dir}/"
log_info "CheckM2:   ${checkm2_dir}/"
log_info "BUSCO:     ${busco_dir}/"
log_info "Merqury:   ${merqury_dir}/"
