#!/usr/bin/env bash
# ==============================================================================
# 06a_predict_genes.sh — Multi-caller gene prediction & consensus merging
# ==============================================================================
# Usage:
#   bash 06a_predict_genes.sh --assembly polished.fasta
#
# Runs Pyrodigal and Glimmer3, optionally imports GeneMarkS-2 predictions,
# then merges gene calls into a consensus set using majority vote (default)
# or union mode.
#
# Optional flags:
#   --genemark-gff FILE      Import GeneMarkS-2 GFF from webserver run
#   --merge-mode MODE        Gene merge strategy: majority (default) or union
#   --overlap-threshold N    Reciprocal overlap for "same gene" (default: 0.80)
#   --threads N              Number of threads (default: 128)
#   --config FILE            Path to pipeline.conf
#   --dry-run                Print commands without executing
#   --help                   Show this help
# ==============================================================================

source "$(dirname "$0")/00_setup.sh"

# --- Defaults ----------------------------------------------------------------
threads="${threads:-128}"
assembly="${assembly:-}"
genemark_gff="${genemark_gff:-}"
gene_merge_mode="${gene_merge_mode:-majority}"
overlap_threshold="${overlap_threshold:-0.80}"

# --- Usage -------------------------------------------------------------------
usage() {
    echo "Usage: $(basename "$0") --assembly FILE [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --assembly FILE              Polished/oriented assembly FASTA"
    echo ""
    echo "Optional:"
    echo "  --genemark-gff FILE          GeneMarkS-2 GFF (from webserver, manual)"
    echo "  --merge-mode MODE            Gene merge mode: majority (default) or union"
    echo "  --overlap-threshold N        Reciprocal overlap threshold (default: 0.80)"
    echo "  --threads N                  Number of threads (default: 128)"
    echo "  --config FILE                Path to pipeline.conf"
    echo "  --dry-run                    Print commands without executing"
    echo "  --help                       Show this help"
    exit 0
}

# --- Parse arguments ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --assembly)           assembly="$2"; shift 2 ;;
        --genemark-gff)       genemark_gff="$2"; shift 2 ;;
        --merge-mode)         gene_merge_mode="$2"; shift 2 ;;
        --overlap-threshold)  overlap_threshold="$2"; shift 2 ;;
        --threads)            threads="$2"; shift 2 ;;
        --config)             load_config "$2"; shift 2 ;;
        --dry-run)            dry_run=true; shift ;;
        --help|-h)            usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

# --- Validate ----------------------------------------------------------------
require_arg "--assembly" "${assembly}"
require_file "${assembly}" "03_polish_orient.sh or 05_annotate_assess.sh"

# Validate merge mode
if [[ "${gene_merge_mode}" != "majority" && "${gene_merge_mode}" != "union" ]]; then
    log_error "--merge-mode must be 'majority' or 'union', got: ${gene_merge_mode}"
fi

# Hard requirements
require_tool pyrodigal
require_tool python3

# Glimmer3 comes as a suite — check for the main binary
if ! command -v glimmer3 &>/dev/null; then
    log_warn "'glimmer3' not found in \$PATH. Glimmer3 predictions will be skipped."
    has_glimmer=false
else
    has_glimmer=true
    for glimmer_dep in long-orfs build-icm extract; do
        require_tool "${glimmer_dep}"
    done
fi

# Check for the merge script
script_dir="$(cd "$(dirname "$0")" && pwd)"
merge_script="${script_dir}/utils/merge_gene_predictions.py"
require_file "${merge_script}" "Ensure scripts/utils/merge_gene_predictions.py exists"

# Optional: GeneMarkS-2 GFF
if [[ -n "${genemark_gff}" ]]; then
    require_file "${genemark_gff}" "Upload assembly to GeneMark webserver and download GFF"
    log_info "GeneMarkS-2 GFF provided: ${genemark_gff}"
else
    log_info "No GeneMarkS-2 GFF provided (--genemark-gff). Proceeding without."
fi

# --- Derived paths -----------------------------------------------------------
predict_dir="$(pwd)/13_gene_prediction"
pyrodigal_dir="${predict_dir}/pyrodigal"
glimmer_dir="${predict_dir}/glimmer3"
consensus_dir="${predict_dir}/consensus"

run_cmd mkdir -p "${pyrodigal_dir}"
run_cmd mkdir -p "${glimmer_dir}"
run_cmd mkdir -p "${consensus_dir}"

# ==============================================================================
# STEP 1 — Pyrodigal Gene Prediction
# ==============================================================================

log_step "Step 1: Pyrodigal gene prediction"

run_cmd pyrodigal \
    -i "${assembly}" \
    -f gff \
    -o "${pyrodigal_dir}/pyrodigal.gff" \
    -a "${pyrodigal_dir}/pyrodigal_proteins.faa" \
    -d "${pyrodigal_dir}/pyrodigal_genes.fna" \
    -p single

pyrodigal_count=0
if [[ -z "${dry_run:-}" && -f "${pyrodigal_dir}/pyrodigal.gff" ]]; then
    pyrodigal_count=$(grep -c $'\tCDS\t' "${pyrodigal_dir}/pyrodigal.gff" 2>/dev/null || echo 0)
fi
log_info "Pyrodigal predicted ${pyrodigal_count} CDS features."

# ==============================================================================
# STEP 2 — Glimmer3 Gene Prediction
# ==============================================================================

if [[ "${has_glimmer}" == "true" ]]; then
    log_step "Step 2: Glimmer3 gene prediction"

    glimmer_tag="${glimmer_dir}/glimmer"

    # Step 2a: Find long ORFs for training
    log_info "2a. Finding long ORFs for ICM training..."
    run_cmd long-orfs -n -t 1.15 "${assembly}" "${glimmer_dir}/long_orfs.txt"

    # Step 2b: Build Interpolated Context Model (ICM)
    log_info "2b. Building ICM from long ORFs..."
    if [[ -z "${dry_run:-}" ]]; then
        run_cmd bash -c 'extract "$1" "$2" | build-icm -r "$3"' \
            _ "${assembly}" "${glimmer_dir}/long_orfs.txt" "${glimmer_dir}/glimmer.icm"
    else
        run_cmd extract "${assembly}" "${glimmer_dir}/long_orfs.txt" \| build-icm -r "${glimmer_dir}/glimmer.icm"
    fi

    # Step 2c: Predict genes
    log_info "2c. Predicting genes with Glimmer3..."
    run_cmd glimmer3 -o50 -g110 -t30 \
        "${assembly}" "${glimmer_dir}/glimmer.icm" "${glimmer_tag}"

    glimmer_count=0
    if [[ -z "${dry_run:-}" && -f "${glimmer_tag}.predict" ]]; then
        glimmer_count=$(grep -c '^[^>]' "${glimmer_tag}.predict" 2>/dev/null || echo 0)
    fi
    log_info "Glimmer3 predicted ${glimmer_count} ORFs."
else
    log_warn "Step 2: Skipping Glimmer3 (not installed)."
fi

# ==============================================================================
# STEP 3 — (Optional) GeneMarkS-2 Import
# ==============================================================================

if [[ -n "${genemark_gff}" ]]; then
    log_step "Step 3: Importing GeneMarkS-2 predictions"

    genemark_target="${predict_dir}/genemark/genemark.gff"
    run_cmd mkdir -p "${predict_dir}/genemark"
    run_cmd cp "${genemark_gff}" "${genemark_target}"

    genemark_count=0
    if [[ -z "${dry_run:-}" && -f "${genemark_target}" ]]; then
        genemark_count=$(grep -c $'\tCDS\t' "${genemark_target}" 2>/dev/null || echo 0)
    fi
    log_info "GeneMarkS-2 GFF contains ${genemark_count} CDS features."
else
    log_info "Step 3: No GeneMarkS-2 GFF — skipping import."
fi

# ==============================================================================
# STEP 4 — Merge Gene Predictions into Consensus
# ==============================================================================

log_step "Step 4: Merging gene predictions (mode: ${gene_merge_mode}, overlap: ${overlap_threshold})"

# Build the argument list for the merge script
merge_args=(
    --assembly "${assembly}"
    --pyrodigal-gff "${pyrodigal_dir}/pyrodigal.gff"
    --merge-mode "${gene_merge_mode}"
    --overlap-threshold "${overlap_threshold}"
    --output-dir "${consensus_dir}"
    --log-file "${consensus_dir}/gene_merge.log"
)

if [[ "${has_glimmer}" == "true" ]]; then
    merge_args+=(--glimmer-predict "${glimmer_dir}/glimmer.predict")
    merge_args+=(--glimmer-detail "${glimmer_dir}/glimmer.detail")
fi

if [[ -n "${genemark_gff}" ]]; then
    merge_args+=(--genemark-gff "${predict_dir}/genemark/genemark.gff")
fi

run_cmd python3 "${merge_script}" "${merge_args[@]}"

# --- Report results ----------------------------------------------------------
consensus_count=0
if [[ -z "${dry_run:-}" && -f "${consensus_dir}/consensus_genes.gff3" ]]; then
    consensus_count=$(grep -c $'\tCDS\t' "${consensus_dir}/consensus_genes.gff3" 2>/dev/null || echo 0)
fi

log_info "Consensus gene set: ${consensus_count} CDS features."

# --- Summary -----------------------------------------------------------------
log_step "Gene prediction complete."
log_info "Pyrodigal GFF: ${pyrodigal_dir}/pyrodigal.gff"
if [[ "${has_glimmer}" == "true" ]]; then
    log_info "Glimmer3:      ${glimmer_dir}/glimmer.predict"
fi
if [[ -n "${genemark_gff}" ]]; then
    log_info "GeneMarkS-2:   ${predict_dir}/genemark/genemark.gff"
fi
log_info "Consensus GFF: ${consensus_dir}/consensus_genes.gff3"
log_info "Consensus FAA: ${consensus_dir}/consensus_proteins.faa"
log_info "Consensus FNA: ${consensus_dir}/consensus_cds.fna"
log_info "Merge log:     ${consensus_dir}/gene_merge.log"
log_info ""
log_info ">>> NEXT: Run 06b_annotate_trackA.sh (assembly-level annotation)"
log_info "         Run 06c_annotate_trackB.sh --consensus-proteins ${consensus_dir}/consensus_proteins.faa"
log_info "         Run 06d_annotate_trackC.sh (prophage detection)"
