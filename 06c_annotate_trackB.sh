#!/usr/bin/env bash
# ==============================================================================
# 06c_annotate_trackB.sh — Track B: Protein-level sequential annotation
# ==============================================================================
# Usage:
#   bash 06c_annotate_trackB.sh --consensus-proteins consensus_proteins.faa
#
# Runs eggNOG-mapper, InterProScan, KofamScan, and optionally DIAMOND
# against SwissProt on the consensus protein set, sequentially.
#
# Optional flags:
#   --interproscan-db DIR    InterProScan data directory
#   --kofam-profiles DIR     KofamScan HMM profiles directory
#   --kofam-ko-list FILE     KofamScan KO list file
#   --swissprot-dmnd FILE    DIAMOND-formatted SwissProt database (optional)
#   --threads N              Number of threads (default: 128)
#   --config FILE            Path to pipeline.conf
#   --dry-run                Print commands without executing
#   --help                   Show this help
# ==============================================================================

source "$(dirname "$0")/00_setup.sh"

# --- Defaults ----------------------------------------------------------------
threads="${threads:-128}"
consensus_proteins="${consensus_proteins:-}"
interproscan_db="${interproscan_db:-}"
kofam_profiles="${kofam_profiles:-}"
kofam_ko_list="${kofam_ko_list:-}"
swissprot_dmnd="${swissprot_dmnd:-}"

# --- Usage -------------------------------------------------------------------
usage() {
    echo "Usage: $(basename "$0") --consensus-proteins FILE [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --consensus-proteins FILE    Consensus protein FASTA from 06a"
    echo ""
    echo "Optional:"
    echo "  --interproscan-db DIR        InterProScan data directory"
    echo "  --kofam-profiles DIR         KofamScan HMM profiles directory"
    echo "  --kofam-ko-list FILE         KofamScan KO list file"
    echo "  --swissprot-dmnd FILE        DIAMOND SwissProt DB (optional; skipped if omitted)"
    echo "  --threads N                  Number of threads (default: 128)"
    echo "  --config FILE                Path to pipeline.conf"
    echo "  --dry-run                    Print commands without executing"
    echo "  --help                       Show this help"
    exit 0
}

# --- Parse arguments ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --consensus-proteins) consensus_proteins="$2"; shift 2 ;;
        --interproscan-db)    interproscan_db="$2"; shift 2 ;;
        --kofam-profiles)     kofam_profiles="$2"; shift 2 ;;
        --kofam-ko-list)      kofam_ko_list="$2"; shift 2 ;;
        --swissprot-dmnd)     swissprot_dmnd="$2"; shift 2 ;;
        --threads)            threads="$2"; shift 2 ;;
        --config)             load_config "$2"; shift 2 ;;
        --dry-run)            dry_run=true; shift ;;
        --help|-h)            usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

# --- Validate ----------------------------------------------------------------
require_arg "--consensus-proteins" "${consensus_proteins}"
require_file "${consensus_proteins}" "06a_predict_genes.sh"

# Hard requirement: eggNOG-mapper
require_tool emapper.py

# Soft requirements
for tool_info in "interproscan.sh:InterProScan" "exec_annotation:KofamScan" "diamond:DIAMOND"; do
    tool_bin="${tool_info%%:*}"
    tool_name="${tool_info##*:}"
    command -v "${tool_bin}" &>/dev/null || log_warn "'${tool_bin}' not found — ${tool_name} will be skipped."
done

# --- Derived paths -----------------------------------------------------------
trackB_dir="$(pwd)/15_trackB"
emapper_dir="${trackB_dir}/eggnog"
ips_dir="${trackB_dir}/interproscan"
kofam_dir="${trackB_dir}/kofamscan"
diamond_dir="${trackB_dir}/diamond_swissprot"

run_cmd mkdir -p "${emapper_dir}"
run_cmd mkdir -p "${ips_dir}"
run_cmd mkdir -p "${kofam_dir}"
run_cmd mkdir -p "${diamond_dir}"

# ==============================================================================
# STEP B1 — eggNOG-mapper
# ==============================================================================

log_step "Step B1: eggNOG-mapper (orthology-based annotation)"

run_cmd emapper.py \
    -i "${consensus_proteins}" \
    --output "${emapper_dir}/emapper_results" \
    --cpu "${threads}" \
    -m diamond \
    --itype proteins \
    --tax_scope auto \
    --go_evidence non-electronic \
    --decorate_gff yes

if [[ -z "${dry_run:-}" && -f "${emapper_dir}/emapper_results.emapper.annotations" ]]; then
    emapper_hits=$(tail -n +5 "${emapper_dir}/emapper_results.emapper.annotations" | grep -cv '^#' 2>/dev/null || echo 0)
    log_info "eggNOG-mapper annotated ${emapper_hits} proteins."
fi

# ==============================================================================
# STEP B2 — InterProScan
# ==============================================================================

if command -v interproscan.sh &>/dev/null; then
    log_step "Step B2: InterProScan (multi-domain annotation)"

    ips_args=(
        -i "${consensus_proteins}"
        -o "${ips_dir}/interproscan_results.tsv"
        -f TSV,GFF3
        -goterms
        -pa
        -cpu "${threads}"
    )

    # Add database path if provided
    if [[ -n "${interproscan_db}" ]]; then
        ips_args+=(-dp "${interproscan_db}")
    fi

    run_cmd interproscan.sh "${ips_args[@]}"

    if [[ -z "${dry_run:-}" && -f "${ips_dir}/interproscan_results.tsv" ]]; then
        ips_hits=$(wc -l < "${ips_dir}/interproscan_results.tsv" 2>/dev/null || echo 0)
        log_info "InterProScan produced ${ips_hits} domain/signature hits."
    fi
else
    log_warn "Step B2: Skipping InterProScan (interproscan.sh not found)."
fi

# ==============================================================================
# STEP B3 — KofamScan
# ==============================================================================

if command -v exec_annotation &>/dev/null; then
    log_step "Step B3: KofamScan (KEGG KO with adaptive thresholds)"

    kofam_args=(
        -o "${kofam_dir}/kofamscan_results.txt"
        --cpu "${threads}"
        -f detail-tsv
    )

    if [[ -n "${kofam_profiles}" ]]; then
        kofam_args+=(--profile "${kofam_profiles}")
    fi
    if [[ -n "${kofam_ko_list}" ]]; then
        kofam_args+=(--ko-list "${kofam_ko_list}")
    fi

    run_cmd exec_annotation "${kofam_args[@]}" "${consensus_proteins}"

    if [[ -z "${dry_run:-}" && -f "${kofam_dir}/kofamscan_results.txt" ]]; then
        # Count significant hits (lines starting with *)
        kofam_sig=$(grep -c '^\*' "${kofam_dir}/kofamscan_results.txt" 2>/dev/null || echo 0)
        log_info "KofamScan found ${kofam_sig} significant KO assignments."
    fi
else
    log_warn "Step B3: Skipping KofamScan (exec_annotation not found)."
fi

# ==============================================================================
# STEP B4 — DIAMOND vs SwissProt (Optional)
# ==============================================================================

if [[ -n "${swissprot_dmnd}" ]] && command -v diamond &>/dev/null; then
    log_step "Step B4: DIAMOND vs SwissProt (curated functional hits)"

    require_file "${swissprot_dmnd}" "Download SwissProt and run: diamond makedb --in uniprot_sprot.fasta --db swissprot"

    run_cmd diamond blastp \
        --query "${consensus_proteins}" \
        --db "${swissprot_dmnd}" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
        --max-target-seqs 5 \
        --evalue 1e-10 \
        --threads "${threads}" \
        --out "${diamond_dir}/diamond_swissprot.tsv"

    if [[ -z "${dry_run:-}" && -f "${diamond_dir}/diamond_swissprot.tsv" ]]; then
        diamond_hits=$(cut -f1 "${diamond_dir}/diamond_swissprot.tsv" | sort -u | wc -l 2>/dev/null || echo 0)
        log_info "DIAMOND found SwissProt hits for ${diamond_hits} proteins."
    fi
elif [[ -z "${swissprot_dmnd}" ]]; then
    log_info "Step B4: DIAMOND/SwissProt skipped (no --swissprot-dmnd provided)."
else
    log_warn "Step B4: Skipping DIAMOND (diamond not found)."
fi

# --- Summary -----------------------------------------------------------------
log_step "Track B annotation complete."
log_info "eggNOG-mapper: ${emapper_dir}/"
log_info "InterProScan:  ${ips_dir}/"
log_info "KofamScan:     ${kofam_dir}/"
log_info "DIAMOND:       ${diamond_dir}/"
log_info ""
log_info ">>> NEXT: Also run 06b_annotate_trackA.sh and 06d_annotate_trackC.sh"
log_info "         Then run 06e_reconcile_merge.sh to merge all results."
