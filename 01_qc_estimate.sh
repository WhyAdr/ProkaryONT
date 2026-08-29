#!/usr/bin/env bash
# ==============================================================================
# 01_qc_estimate.sh — Read QC, genome size estimation, k-mer profiling
# ==============================================================================
# Usage:
#   bash 01_qc_estimate.sh --input-fastq reads.fq.gz --sequencing-summary summary.txt
#
# Optional flags:
#   --config FILE           Path to pipeline.conf (values override defaults)
#   --threads N             Number of threads (default: 128)
#   --memory N              Memory in GB for Meryl (default: auto from system)
#   --kmer-size N           K-mer size for Meryl (default: 21)
#   --nanoplot-color C      Color for NanoPlot graphs (default: green)
#   --lint-threshold N      SDUST repetition threshold for profiling (default: 20)
#   --lint-window N         SDUST window size for profiling (default: 64)
#   --dry-run               Print commands without executing
#   -h, --help              Show this help
# ==============================================================================

source "$(dirname "$0")/00_setup.sh"

# --- Defaults ----------------------------------------------------------------
threads="${threads:-}"
input_fastq="${input_fastq:-}"
sequencing_summary="${sequencing_summary:-}"
meryl_memory="${meryl_memory:-}"
meryl_kmer_size="${meryl_kmer_size:-}"
nanoplot_color="${nanoplot_color:-}"
lrge_seed="${lrge_seed:-}"
# Keep these empty until after config loading so Stage 1 has the intended
# precedence: CLI > pipeline.conf > hard-coded default.
lint_threshold="${lint_threshold:-}"
lint_window="${lint_window:-}"
config_file="${config_file:-}"

# --- Usage -------------------------------------------------------------------
usage() {
    echo "Usage: $(basename "$0") --input-fastq FILE --sequencing-summary FILE [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --input-fastq FILE          Raw FASTQ file (gzipped ok)"
    echo "  --sequencing-summary FILE   Sequencing summary from MinKNOW"
    echo ""
    echo "Optional:"
    echo "  --config FILE               Path to pipeline.conf"
    echo "  --threads N                 Number of threads (default: 128)"
    echo "  --memory N                  Memory in GB for Meryl (default: auto)"
    echo "  --kmer-size N               K-mer size for Meryl (default: 21)"
    echo "  --nanoplot-color C          Color for NanoPlot graphs (default: green)"
    echo "  --lint-threshold N          SDUST repetition threshold for profiling (default: 20)"
    echo "  --lint-window N             SDUST window size for profiling (default: 64)"
    echo "  --dry-run                   Print commands without executing"
    echo "  -h, --help                  Show this help"
    exit 0
}

# --- Parse arguments ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-fastq)        input_fastq="$2"; shift 2 ;;
        --sequencing-summary) sequencing_summary="$2"; shift 2 ;;
        --config)             config_file="$2"; shift 2 ;;
        --threads)            threads="$2"; shift 2 ;;
        --memory)             meryl_memory="$2"; shift 2 ;;
        --kmer-size)          meryl_kmer_size="$2"; shift 2 ;;
        --nanoplot-color)     nanoplot_color="$2"; shift 2 ;;
        --lint-threshold)     lint_threshold="$2"; shift 2 ;;
        --lint-window)        lint_window="$2"; shift 2 ;;
        --dry-run)            dry_run=true; shift ;;
        --help|-h)            usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

# --- Load config (CLI flags parsed above take priority) ----------------------
[[ -n "${config_file}" ]] && load_config "${config_file}"

threads="${threads:-128}"
meryl_kmer_size="${meryl_kmer_size:-21}"
nanoplot_color="${nanoplot_color:-green}"
lrge_seed="${lrge_seed:-123}"
lint_threshold="${lint_threshold:-20}"
lint_window="${lint_window:-64}"

# Set GZIP_BIN now that threads are parsed
export GZIP_BIN="$(get_gzip_cmd)"

# --- Validate ----------------------------------------------------------------
require_arg "--input-fastq" "${input_fastq}"
require_arg "--sequencing-summary" "${sequencing_summary}"
require_file "${input_fastq}" ""
require_file "${sequencing_summary}" ""

require_tool NanoPlot
require_tool lrge
require_tool raven
require_tool meryl
require_tool fastcat
require_tool porechop_abi
require_tool snikt.R

# --- Derived paths -----------------------------------------------------------
qc_dir="$(pwd)/01_qc"
genome_size_dir="$(pwd)/02_genome_size"

# ==============================================================================
# STEP 1 — NanoPlot QC
# ==============================================================================

log_step "Step 1: NanoPlot QC"
mkdir -p "${qc_dir}"

if [[ -f "${qc_dir}/run_summary_profile/NanoPlot-report.html" ]]; then
    log_info "NanoPlot summary report already exists. Skipping..."
else
    log_info "Running NanoPlot on sequencing summary..."
    run_cmd NanoPlot --threads "${threads}" \
        --summary "${sequencing_summary}" \
        --loglength \
        -o "${qc_dir}/run_summary_profile"
fi

if [[ -f "${qc_dir}/NanoPlot_sample/NanoPlot-report.html" ]]; then
    log_info "NanoPlot sample report already exists. Skipping..."
else
    log_info "Running NanoPlot on raw FASTQ..."
    run_cmd NanoPlot --threads "${threads}" \
        -c "${nanoplot_color}" --fastq "${input_fastq}" \
        -o "${qc_dir}/NanoPlot_sample" \
        --loglength --plots dot --N50
fi

log_info ">>> CHECK: Open ${qc_dir}/NanoPlot_sample/NanoPlot-report.html"

# ==============================================================================
# STEP 2 — Preprocessing Diagnostics: Fastcat, SDUST, Porechop_ABI Scan & SNIKT
# ==============================================================================

log_step "Step 2: Preprocessing diagnostics"

# --- Fastcat: fast per-file stats + length/quality histograms ---------------
log_info "Running Fastcat per-file summary statistics..."
mkdir -p "${qc_dir}/fastcat_histograms"
run_cmd bash -c 'fastcat --histograms "$1" -f "$2" "$3" > /dev/null' \
    _ "${qc_dir}/fastcat_histograms" "${qc_dir}/fastcat_per_file_stats.tsv" "${input_fastq}"

# --- Fastcat lint: exhaustive SDUST low-complexity profile ------------------
# This diagnostic pass deliberately uses max-proportion=0.0. Fastcat sends
# every SDUST-positive read to stderr and writes only zero-SDUST reads to
# stdout, which gives an exhaustive read-level partition without retaining a
# second FASTQ. Fastcat reports fractions rounded to two decimal places, so do
# not infer exact masked-base counts from this output.
sdust_dir="${qc_dir}/fastcat_sdust"
mkdir -p "${sdust_dir}"
log_info "Profiling read low-complexity burden with Fastcat SDUST (T=${lint_threshold}, W=${lint_window}; diagnostic only)..."

if [[ -n "${dry_run:-}" ]]; then
    run_cmd bash -c 'set -o pipefail; fastcat lint --threshold "$1" --window "$2" --max-proportion 0.0 "$3" 2> "$4" | awk "END { print int(NR / 4) }" > "$5"' \
        _ "${lint_threshold}" "${lint_window}" "${input_fastq}" "${sdust_dir}/lint_p0.stderr.log" "${sdust_dir}/zero_sdust_reads.txt"
    log_info "[DRY-RUN] SDUST profile artifacts would be written to ${sdust_dir}/"
else
    run_cmd bash -c 'set -o pipefail; fastcat lint --threshold "$1" --window "$2" --max-proportion 0.0 "$3" 2> "$4" | awk "END { print int(NR / 4) }" > "$5"' \
        _ "${lint_threshold}" "${lint_window}" "${input_fastq}" "${sdust_dir}/lint_p0.stderr.log" "${sdust_dir}/zero_sdust_reads.txt"

    {
        printf 'read_id\treported_masked_fraction\n'
        sed -nE 's/.*Read ([^ ]+) masked fraction ([0-9.]+) exceeds threshold [0-9.]+, skipping\..*/\1\t\2/p' \
            "${sdust_dir}/lint_p0.stderr.log"
    } > "${sdust_dir}/sdust_positive_reads.tsv"

    zero_sdust_reads="$(awk 'NR == 1 { print $1; exit }' "${sdust_dir}/zero_sdust_reads.txt")"
    zero_sdust_reads="${zero_sdust_reads:-0}"
    sdust_positive_reads="$(awk 'NR > 1 { count++ } END { print count + 0 }' "${sdust_dir}/sdust_positive_reads.tsv")"
    skipped_log_lines="$(awk '/skipping\./ { count++ } END { print count + 0 }' "${sdust_dir}/lint_p0.stderr.log")"

    if [[ "${skipped_log_lines}" -ne "${sdust_positive_reads}" ]]; then
        log_error "Could not parse every Fastcat lint rejection from ${sdust_dir}/lint_p0.stderr.log; no SDUST summary was produced."
    fi

    total_sdust_reads=$((zero_sdust_reads + sdust_positive_reads))
    zero_sdust_fraction="$(awk -v zero="${zero_sdust_reads}" -v total="${total_sdust_reads}" 'BEGIN { if (total) printf "%.6f", zero / total; else printf "0.000000" }')"
    sdust_positive_fraction="$(awk -v positive="${sdust_positive_reads}" -v total="${total_sdust_reads}" 'BEGIN { if (total) printf "%.6f", positive / total; else printf "0.000000" }')"
    reported_gt_050="$(awk -F '\t' 'NR > 1 && $2 > 0.50 { count++ } END { print count + 0 }' "${sdust_dir}/sdust_positive_reads.tsv")"
    reported_gt_080="$(awk -F '\t' 'NR > 1 && $2 > 0.80 { count++ } END { print count + 0 }' "${sdust_dir}/sdust_positive_reads.tsv")"
    reported_gt_090="$(awk -F '\t' 'NR > 1 && $2 > 0.90 { count++ } END { print count + 0 }' "${sdust_dir}/sdust_positive_reads.tsv")"
    reported_gt_095="$(awk -F '\t' 'NR > 1 && $2 > 0.95 { count++ } END { print count + 0 }' "${sdust_dir}/sdust_positive_reads.tsv")"
    reported_fraction_max="$(awk -F '\t' 'NR > 1 && $2 > max { max = $2 } END { printf "%.2f", max + 0 }' "${sdust_dir}/sdust_positive_reads.tsv")"

    {
        printf 'masked_fraction\tread_count\n'
        printf '0.00\t%s\n' "${zero_sdust_reads}"
        awk -F '\t' 'NR > 1 { count[$2]++ } END { for (fraction in count) print fraction, count[fraction] }' OFS='\t' \
            "${sdust_dir}/sdust_positive_reads.tsv" | sort -n
    } > "${sdust_dir}/sdust_fraction_hist.tsv"

    {
        head -n 1 "${sdust_dir}/sdust_positive_reads.tsv"
        awk -F '\t' 'NR > 1 && $2 >= 0.80' "${sdust_dir}/sdust_positive_reads.tsv" | sort -t $'\t' -k2,2nr
    } > "${sdust_dir}/sdust_high_burden_reads.tsv"

    {
        printf 'metric\tvalue\n'
        printf 'sdust_threshold\t%s\n' "${lint_threshold}"
        printf 'sdust_window\t%s\n' "${lint_window}"
        printf 'total_reads\t%s\n' "${total_sdust_reads}"
        printf 'zero_sdust_reads\t%s\n' "${zero_sdust_reads}"
        printf 'sdust_positive_reads\t%s\n' "${sdust_positive_reads}"
        printf 'zero_sdust_fraction\t%s\n' "${zero_sdust_fraction}"
        printf 'sdust_positive_fraction\t%s\n' "${sdust_positive_fraction}"
        printf 'reported_fraction_gt_0.50\t%s\n' "${reported_gt_050}"
        printf 'reported_fraction_gt_0.80\t%s\n' "${reported_gt_080}"
        printf 'reported_fraction_gt_0.90\t%s\n' "${reported_gt_090}"
        printf 'reported_fraction_gt_0.95\t%s\n' "${reported_gt_095}"
        printf 'reported_fraction_max\t%s\n' "${reported_fraction_max}"
    } > "${sdust_dir}/sdust_summary.tsv"

    log_info "SDUST complexity profile: total=${total_sdust_reads}; zero-SDUST=${zero_sdust_reads}; SDUST-positive=${sdust_positive_reads}; reported >0.80=${reported_gt_080}; >0.90=${reported_gt_090}; >0.95=${reported_gt_095}; max=${reported_fraction_max}"
    log_info ">>> CHECK: ${sdust_dir}/sdust_summary.tsv, sdust_fraction_hist.tsv, sdust_high_burden_reads.tsv"
fi

# --- Porechop_ABI: ab-initio adapter discovery scan (report only) -----------
log_info "Running Porechop_ABI ab-initio adapter discovery (report only, no trimming)..."
# NOTE: Porechop_ABI is not barcode-aware by design (its own docs say
# demultiplexing should be handled upstream). Its ab-initio module detects
# over-represented k-mers at read termini agnostically, though, so it may
# incidentally flag residual barcode sequence sitting just inside the
# adapter boundary as a side effect — treat any such hits as a bonus, not
# a guaranteed capability.
run_cmd bash -c 'porechop_abi -abi --guess_adapter_only -i "$1" --threads "$2" > "$3" 2>&1' \
    _ "${input_fastq}" "${threads}" "${qc_dir}/porechop_abi_scan.log"

# --- SNIKT: baseline systemic end-bias report (pre-trim) --------------------
log_info "Running SNIKT baseline contamination check (--notrim, diagnostic only)..."
mkdir -p "${qc_dir}/snikt_baseline"
run_cmd bash -c 'cd "$1" && snikt.R --notrim "$2"' \
    _ "${qc_dir}/snikt_baseline" "${input_fastq}"

log_info ">>> CHECK: Review ${qc_dir}/fastcat_per_file_stats.tsv, fastcat_sdust/, porechop_abi_scan.log, snikt_baseline/ — these inform whether Stage 2 trimming is warranted"
log_info ">>> CHECK SDUST tail before enabling Stage 2 Fastcat filtering."

# ==============================================================================
# STEP 3 — Genome Size Estimation
# ==============================================================================

log_step "Step 3: Genome size estimation"
mkdir -p "${genome_size_dir}"

log_info "Running LRGE..."
if [[ -s "${genome_size_dir}/lrge_output.txt" ]]; then
    lrge_size=$(cat "${genome_size_dir}/lrge_output.txt")
    log_info "Found existing LRGE estimate: ${lrge_size}. Skipping..."
else
    run_cmd bash -c 'lrge -P ont -t "$1" -s "$2" "$3" > "$4"' \
        _ "${threads}" "${lrge_seed}" "${input_fastq}" "${genome_size_dir}/lrge_output.txt"
    if [[ -z "${dry_run:-}" && -f "${genome_size_dir}/lrge_output.txt" ]]; then
        lrge_size=$(cat "${genome_size_dir}/lrge_output.txt")
        log_info "LRGE estimated genome size: ${lrge_size}"
    fi
fi

log_info "Running Raven quick assembly (0 polish)..."
if [[ -s "${genome_size_dir}/assembly.fasta" ]]; then
    log_info "Found existing Raven assembly. Skipping..."
else
    run_cmd bash -c 'raven --threads "$1" -p 0 \
        --graphical-fragment-assembly "$2" \
        "$3" > "$4"' \
        _ "${threads}" "${genome_size_dir}/assemblyGraph.gfa" "${input_fastq}" "${genome_size_dir}/assembly.fasta"
fi

if [[ -z "${dry_run:-}" && -f "${genome_size_dir}/assembly.fasta" ]]; then
    raven_size=$(grep -v '^>' "${genome_size_dir}/assembly.fasta" | tr -d '\r\n' | wc -c)
    log_info "Raven assembly size: ${raven_size} bp"
    run_cmd rm -f "${genome_size_dir}/assemblyGraph.gfa"
fi

if [[ -n "${lrge_size:-}" && -n "${raven_size:-}" ]]; then
    # Heuristic: Raven assembly length is generally more accurate than LRGE k-mer estimation
    mean_genome_size=$(( (2 * raven_size + lrge_size) / 3 ))
    echo "${mean_genome_size}" > "${genome_size_dir}/mean_genome_size.txt"
    log_info "Weighted Mean Genome Size ((2*Raven + LRGE) / 3): ${mean_genome_size} bp"
fi

log_info "Running Meryl k-mer counting..."
if [[ -z "${meryl_memory}" ]]; then
    meryl_memory=$( { free -g 2>/dev/null || echo "Mem: 0 0 0 0 0 0 16"; } | awk '/^Mem:/{print ($7 != "" ? $7 : $4)}')
    meryl_memory="${meryl_memory:-16}"
    [[ "${meryl_memory}" -eq 0 ]] && meryl_memory=16
fi

if [[ -d "${genome_size_dir}/genome.meryl" ]]; then
    log_info "Found existing Meryl database. Skipping k-mer counting..."
else
    run_cmd meryl count compress k="${meryl_kmer_size}" memory="${meryl_memory}" threads="${threads}" \
        output "${genome_size_dir}/genome.meryl" "${input_fastq}"
fi

if [[ -d "${genome_size_dir}/genome.meryl" && ! -f "${genome_size_dir}/genome.hist" ]]; then
    run_cmd bash -c 'meryl histogram "$1" > "$2"' \
        _ "${genome_size_dir}/genome.meryl" "${genome_size_dir}/genome.hist"
else
    log_info "Meryl histogram already exists. Skipping..."
fi

# --- Summary -----------------------------------------------------------------
log_step "QC & estimation complete."
log_info "LRGE estimate:  ${lrge_size:-N/A}"
log_info "Raven size:     ${raven_size:-N/A} bp"
log_info "Mean (Weight):  ${mean_genome_size:-N/A} bp"
log_info "NanoPlot:       ${qc_dir}/"
log_info "Meryl db:       ${genome_size_dir}/genome.meryl"
log_info "Diagnostics:    ${qc_dir}/ (fastcat_per_file_stats.tsv / fastcat_sdust/ / porechop_abi_scan.log / snikt_baseline/)"
