#!/usr/bin/env bash
# ==============================================================================
# 02_preprocess_filter.sh — Adapter trimming, end-cropping & quality/length filtering
# ==============================================================================
# Usage:
#   bash 02_preprocess_filter.sh --input-fastq reads.fq.gz
#
# Optional flags:
#   --config FILE                 Path to pipeline.conf (values override defaults)
#   --threads N                   Number of threads (default: 128)
#   --min-length N                Filtlong minimum read length (default: 200)
#   --min-qscore N                Seqkit minimum average read quality (default: 7)
#   --keep-percent N              Filtlong keep percent (default: 90)
#   --enable-porechopabi-trim     Enable Porechop_ABI adapter trimming (default: off)
#   --enable-chopper-trim         Enable Chopper end-trimming (default: off — gate on Stage 1/2 SNIKT+Porechop_ABI triage)
#   --chopper-trim-approach NAME  fixed-crop | trim-by-quality (default: fixed-crop; verify full list via `chopper --help`)
#   --headcrop N                  Chopper fixed-crop: bases to crop from read start (default: 50)
#   --tailcrop N                  Chopper fixed-crop: bases to crop from read end (default: 30)
#   --chopper-trim-cutoff N       Chopper trim-by-quality cutoff (verify exact chopper flag before use)
#   --dry-run                     Print commands without executing
#   --help                        Show this help
# ==============================================================================

source "$(dirname "$0")/00_setup.sh"

# --- Defaults ----------------------------------------------------------------
threads="${threads:-128}"
input_fastq="${input_fastq:-}"
filtlong_min_length="${filtlong_min_length:-200}"
min_qscore="${min_qscore:-7}"
filtlong_keep_percent="${filtlong_keep_percent:-90}"
config_file="${config_file:-}"
enable_porechopabi_trim="${enable_porechopabi_trim:-false}"
enable_chopper_trim="${enable_chopper_trim:-false}"
chopper_trim_approach="${chopper_trim_approach:-fixed-crop}"
chopper_headcrop="${chopper_headcrop:-50}"
chopper_tailcrop="${chopper_tailcrop:-30}"
chopper_trim_cutoff="${chopper_trim_cutoff:-}"
chopper_split_window="${chopper_split_window:-1}"
enable_fastcat_lint="${enable_fastcat_lint:-false}"
lint_threshold="${lint_threshold:-20}"
lint_window="${lint_window:-64}"
lint_max_proportion="${lint_max_proportion:-0.95}"

# --- Usage -------------------------------------------------------------------
usage() {
    echo "Usage: $(basename "$0") --input-fastq FILE [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --input-fastq FILE          Raw FASTQ file (gzipped ok)"
    echo ""
    echo "Optional:"
    echo "  --config FILE               Path to pipeline.conf"
    echo "  --threads N                 Number of threads (default: 128)"
    echo "  --min-length N              Filtlong min read length in bp (default: 200)"
    echo "  --min-qscore N              Seqkit min average read quality (default: 7)"
    echo "  --keep-percent N            Filtlong keep percent (default: 90)"
    echo "  --enable-porechopabi-trim   Enable Porechop_ABI adapter trimming"
    echo "  --enable-chopper-trim       Enable Chopper end-trimming"
    echo "  --chopper-trim-approach N   fixed-crop | trim-by-quality | best-read-segment | split-by-low-quality (default: fixed-crop)"
    echo "  --headcrop N                Chopper headcrop bases (default: 50)"
    echo "  --tailcrop N                Chopper tailcrop bases (default: 30)"
    echo "  --chopper-trim-cutoff N     Chopper quality cutoff/threshold (required for quality-aware modes)"
    echo "  --split-window N            Min consecutive low-quality bases to split (default: 1)"
    echo "  --enable-fastcat-lint       Enable Fastcat lint low-complexity filtering"
    echo "  --lint-threshold N          DUST threshold T (default: 20)"
    echo "  --lint-window N             DUST window size W (default: 64)"
    echo "  --lint-max-proportion N     Max tolerated masked fraction to keep read (default: 0.95)"
    echo "  --dry-run                   Print commands without executing"
    echo "  --help                      Show this help"
    exit 0
}

# --- Parse arguments ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-fastq)            input_fastq="$2"; shift 2 ;;
        --config)                 config_file="$2"; shift 2 ;;
        --threads)                threads="$2"; shift 2 ;;
        --min-length)             filtlong_min_length="$2"; shift 2 ;;
        --min-qscore)             min_qscore="$2"; shift 2 ;;
        --keep-percent)           filtlong_keep_percent="$2"; shift 2 ;;
        --enable-porechopabi-trim) enable_porechopabi_trim=true; shift ;;
        --enable-chopper-trim)    enable_chopper_trim=true; shift ;;
        --chopper-trim-approach)  chopper_trim_approach="$2"; shift 2 ;;
        --headcrop)               chopper_headcrop="$2"; shift 2 ;;
        --tailcrop)               chopper_tailcrop="$2"; shift 2 ;;
        --chopper-trim-cutoff)    chopper_trim_cutoff="$2"; shift 2 ;;
        --split-window)           chopper_split_window="$2"; shift 2 ;;
        --enable-fastcat-lint)    enable_fastcat_lint=true; shift ;;
        --lint-threshold)         lint_threshold="$2"; shift 2 ;;
        --lint-window)            lint_window="$2"; shift 2 ;;
        --lint-max-proportion)    lint_max_proportion="$2"; shift 2 ;;
        --dry-run)                dry_run=true; shift ;;
        --help|-h)                usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

# --- Load config (CLI flags parsed above take priority) ----------------------
[[ -n "${config_file}" ]] && load_config "${config_file}"

# Set GZIP_BIN now that threads are parsed
export GZIP_BIN="$(get_gzip_cmd)"

# --- Validate ----------------------------------------------------------------
require_arg "--input-fastq" "${input_fastq}"
require_file "${input_fastq}" ""

require_tool filtlong
require_tool seqkit
require_tool NanoPlot
require_tool snikt.R
require_tool fastplong
[[ "${enable_porechopabi_trim}" == "true" ]] && require_tool porechop_abi
[[ "${enable_chopper_trim}" == "true" ]] && require_tool chopper
[[ "${enable_fastcat_lint}" == "true" ]] && require_tool fastcat

# --- Derived paths -----------------------------------------------------------
qc_dir="$(pwd)/01_qc"
filtered_reads="$(pwd)/filtered_input.fastq.gz"

# --- Temporary files for cleanup ---------------------------------------------
_trim_tmp=""
_crop_tmp=""
_lint_tmp=""
_qfilt_tmp=""

cleanup() {
    [[ -n "${_trim_tmp:-}" && -f "${_trim_tmp}" ]] && rm -f "${_trim_tmp}"
    [[ -n "${_crop_tmp:-}" && -f "${_crop_tmp}" ]] && rm -f "${_crop_tmp}"
    [[ -n "${_lint_tmp:-}" && -f "${_lint_tmp}" ]] && rm -f "${_lint_tmp}"
    [[ -n "${_qfilt_tmp:-}" && -f "${_qfilt_tmp}" ]] && rm -f "${_qfilt_tmp}"
    true
}
trap cleanup EXIT

# ==============================================================================
# STEP 4 — Porechop_ABI Adapter Trimming (optional)
# ==============================================================================

_preproc_input="${input_fastq}"

if [[ "${enable_porechopabi_trim}" == "true" ]]; then
    log_step "Step 4: Trimming adapters with Porechop_ABI..."
    # See 01_qc_estimate.sh STEP 2 for the barcode side-effect caveat.
    _trim_tmp="$(mktemp --suffix=.fastq)"
    run_cmd porechop_abi -abi -i "${_preproc_input}" -o "${_trim_tmp}" --threads "${threads}"
    _preproc_input="${_trim_tmp}"

    # ==========================================================================
    # STEP 5 — SNIKT Re-Assessment (post-trim)
    # ==========================================================================
    # Only runs when Porechop_ABI actually trimmed something — re-running
    # against unchanged reads would just reproduce the Stage 1 baseline.
    log_step "Step 5: Re-assessing systemic end-bias with SNIKT (post Porechop_ABI trim)..."
    mkdir -p "${qc_dir}/snikt_post_porechop"
    run_cmd bash -c 'cd "$1" && snikt.R --notrim "$2"' \
        _ "${qc_dir}/snikt_post_porechop" "${_preproc_input}"
    log_info ">>> CHECK: Compare ${qc_dir}/snikt_post_porechop/ to ${qc_dir}/snikt_baseline/ to decide whether Chopper end-trimming is still warranted."
fi

# ==============================================================================
# STEP 6a — Chopper End-Trimming (optional, triage-gated)
# ==============================================================================

if [[ "${enable_chopper_trim}" == "true" ]]; then
    log_step "Step 6a: End-trimming with Chopper (approach=${chopper_trim_approach})..."
    _crop_tmp="$(mktemp --suffix=.fastq)"
    case "${chopper_trim_approach}" in
        fixed-crop)
            run_cmd bash -c 'zcat -f "$1" | chopper --trim-approach fixed-crop --headcrop "$2" --tailcrop "$3" --threads "$4" > "$5"' \
                _ "${_preproc_input}" "${chopper_headcrop}" "${chopper_tailcrop}" "${threads}" "${_crop_tmp}"
            ;;
        trim-by-quality|best-read-segment)
            require_arg "--chopper-trim-cutoff" "${chopper_trim_cutoff}"
            run_cmd bash -c 'zcat -f "$1" | chopper --trim-approach "$2" --cutoff "$3" --threads "$4" > "$5"' \
                _ "${_preproc_input}" "${chopper_trim_approach}" "${chopper_trim_cutoff}" "${threads}" "${_crop_tmp}"
            ;;
        split-by-low-quality)
            require_arg "--chopper-trim-cutoff" "${chopper_trim_cutoff}"
            run_cmd bash -c 'zcat -f "$1" | chopper --trim-approach split-by-low-quality --cutoff "$2" --split-window "$3" --threads "$4" > "$5"' \
                _ "${_preproc_input}" "${chopper_trim_cutoff}" "${chopper_split_window}" "${threads}" "${_crop_tmp}"
            ;;
        *)
            log_error "Unknown chopper_trim_approach: ${chopper_trim_approach}. Verify the full option list via 'chopper --help'."
            ;;
    esac
    [[ "${_preproc_input}" != "${input_fastq}" ]] && rm -f "${_preproc_input}"
    _preproc_input="${_crop_tmp}"
fi

# ==============================================================================
# STEP 6b — Fastcat Lint Low-Complexity Filtering (optional)
# ==============================================================================

if [[ "${enable_fastcat_lint}" == "true" ]]; then
    log_step "Step 6b: Low-complexity filtering with Fastcat lint..."
    _lint_tmp="$(mktemp --suffix=.fastq)"
    mkdir -p "${qc_dir}"
    
    run_cmd bash -c 'fastcat lint --threshold "$1" --window "$2" --max-proportion "$3" "$4" 2>>"$5" > "$6"' \
        _ "${lint_threshold}" "${lint_window}" "${lint_max_proportion}" "${_preproc_input}" "${qc_dir}/02_lint_rejected.log" "${_lint_tmp}"
    
    [[ "${_preproc_input}" != "${input_fastq}" ]] && rm -f "${_preproc_input}"
    _preproc_input="${_lint_tmp}"
    
    rejected_count=$(grep -c "skipping" "${qc_dir}/02_lint_rejected.log" 2>/dev/null || echo 0)
    log_info "Fastcat lint completed. Rejected ${rejected_count} low-complexity reads (logged to 01_qc/02_lint_rejected.log)."
fi

# ==============================================================================
# STEP 7 — Filtering with Seqkit + Filtlong (logic unchanged from the original
# 02_filter_assemble.sh STEP 2 — preserve the rationale comment below verbatim)
# ==============================================================================

log_step "Step 7: Filtering reads (min_qscore=${min_qscore}, min_length=${filtlong_min_length}, keep=${filtlong_keep_percent}%)"

if [[ -s "${filtered_reads}" ]] && gzip -t "${filtered_reads}" 2>/dev/null; then
    log_info "Found valid existing ${filtered_reads}. Skipping filtering..."
else
    if [[ -f "${filtered_reads}" ]]; then
        log_warn "Existing filtered reads are incomplete or corrupt. Overwriting..."
        rm -f "${filtered_reads}"
    fi

    # Step 7a: Quality-filter with seqkit.
    # NOTE: Seqkit is used here because its `--min-qual` flag enforces a
    # probability-space mean Phred threshold filter (sum error probabilities,
    # average, convert back via -10*log10) — the same approach Chopper's
    # ave_qual() uses internally. Fastplong's passFilter() instead sums raw
    # integer Phred values with integer division (flat arithmetic mean,
    # truncated) — a different, less precise calculation, which is why it is
    # NOT used here despite being available (see STEP 8 for its QC-only role).
    # This is distinct from Filtlong's `--min_mean_q`, which uses a custom
    # Z-score distribution metric that can let reads with small patches of
    # high-quality bases mask overall poor quality.
    # Writing to a temporary file is intentional: depending on compile-time
    # flags, certain builds/versions of Filtlong can crash or silently hang
    # when reading from standard input (stdin). Using an intermediate file
    # has been validated as the most reliable cross-platform solution.
    _qfilt_tmp="$(mktemp --suffix=.fastq)"
    run_cmd seqkit seq --min-qual "${min_qscore}" "${_preproc_input}" -o "${_qfilt_tmp}"

    # Step 7b: Length / keep-percent filter with filtlong
    run_cmd bash -c 'set -o pipefail; filtlong --min_length "$1" --keep_percent "$2" "$3" | ${GZIP_BIN} > "$4"' \
        _ "${filtlong_min_length}" "${filtlong_keep_percent}" "${_qfilt_tmp}" "${filtered_reads}"

    rm -f "${_qfilt_tmp}"
fi

# Clean up any remaining preprocessing temp files
[[ "${_preproc_input}" != "${input_fastq}" ]] && rm -f "${_preproc_input}"

# ==============================================================================
# STEP 8 — Post-Filter QC: NanoPlot + Fastplong Report
# ==============================================================================

log_info "Running NanoPlot on filtered reads..."
mkdir -p "${qc_dir}"
run_cmd NanoPlot --threads "${threads}" -c green \
    --fastq "${filtered_reads}" \
    -o "${qc_dir}/NanoPlot_FiltPol" \
    --loglength --plots dot --N50

log_info "Running Fastplong for QC reporting only — its Q-score formula is a flat arithmetic mean of raw Phred integers with integer truncation, which diverges from Seqkit's probability-space mean (see STEP 7 comment), so it is never used for filtering decisions here."
run_cmd fastplong -i "${filtered_reads}" -o /dev/null \
    --html "${qc_dir}/fastplong_report.html" --json "${qc_dir}/fastplong_report.json" \
    --thread "${threads}"

log_info ">>> CHECK: Compare ${qc_dir}/NanoPlot_FiltPol/ to ${qc_dir}/NanoPlot_sample/; review ${qc_dir}/fastplong_report.html"

# --- Summary -----------------------------------------------------------------
log_step "Preprocessing & filtering complete."
log_info "Filtered reads: ${filtered_reads}"
log_info "NanoPlot:       ${qc_dir}/NanoPlot_FiltPol/"
log_info "Fastplong:      ${qc_dir}/fastplong_report.html"
