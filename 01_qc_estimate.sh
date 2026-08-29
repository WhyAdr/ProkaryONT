#!/usr/bin/env bash
# Stage 1: non-destructive read QC and independent genome-size diagnostics.

source "$(dirname "$0")/00_setup.sh"
# Do not create the legacy pipeline.log while parsing/validating this stage.
log_file=""

threads="${threads:-}"
input_fastq="${input_fastq:-}"
sequencing_summary="${sequencing_summary:-}"
sample_name="${sample_name:-}"
output_dir="${output_dir:-}"
tmp_dir="${tmp_dir:-}"
config_file="${config_file:-}"
meryl_memory="${meryl_memory:-}"
meryl_kmer_size="${meryl_kmer_size:-}"
nanoplot_color="${nanoplot_color:-}"
lrge_seed="${lrge_seed:-}"
lint_threshold="${lint_threshold:-}"
lint_window="${lint_window:-}"
lint_max_proportion="${lint_max_proportion:-}"
genome_size_disagreement_threshold="${genome_size_disagreement_threshold:-}"
input_sha256="${input_sha256:-}"
force="${force:-}"
skip_nanoplot="${skip_nanoplot:-}"
skip_fastcat_stats="${skip_fastcat_stats:-}"
skip_sdust="${skip_sdust:-}"
skip_porechop_scan="${skip_porechop_scan:-}"
skip_snikt="${skip_snikt:-}"
skip_lrge="${skip_lrge:-}"
skip_raven="${skip_raven:-}"
skip_meryl="${skip_meryl:-}"

usage() {
    cat <<'EOF'
Usage: 01_qc_estimate.sh --input-fastq FILE [OPTIONS]

Required:
  --input-fastq FILE             Raw FASTQ (plain or gzip-compressed)

Inputs and layout:
  --sequencing-summary FILE      Optional MinKNOW/Guppy sequencing summary
  --config FILE                  Path to pipeline.conf
  --output-dir DIR               Output root (default: current directory)
  --sample-name NAME             Manifest/sample label (default: input basename)
  --tmp-dir DIR                  Temporary directory (default: OUTPUT/.prokaryont_tmp)

Resources and diagnostics:
  --threads N                    Threads (default: SLURM allocation, then nproc)
  --memory N                     Meryl memory in GB (default: 70% of visible allocation)
  --kmer-size N                  Meryl k-mer size (default: 21)
  --nanoplot-color COLOR         NanoPlot color (default: green)
  --lint-threshold N             SDUST repetition threshold T (default: 20)
  --lint-window N                SDUST window W (default: 64)
  --lint-max-proportion N        Exact count-only SDUST preview threshold (default: 0.95)
  --genome-disagreement N        LRGE/Raven relative warning threshold (default: 0.20)

Default-on diagnostic skip flags:
  --skip-nanoplot                Skip FASTQ and sequencing-summary NanoPlot reports
  --skip-fastcat-stats           Skip Fastcat 1.x FASTQ statistics
  --skip-sdust                   Skip diagnostic SDUST census and preview
  --skip-porechop-scan           Skip Porechop_ABI adapter discovery scan
  --skip-snikt                   Skip SNIKT baseline report
  --skip-lrge                    Skip LRGE overlap estimate
  --skip-raven                   Skip unpolished Raven assembly-span cross-check
  --skip-meryl                   Skip Meryl k-mer profile

Execution policy:
  --input-sha256                 Include full input SHA-256 hashes in identity
  --force                        Regenerate outputs after a missing/mismatched manifest
  --dry-run                      Validate and print commands without filesystem mutation
  -h, --help                     Show this help
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-fastq|--sequencing-summary|--config|--output-dir|--sample-name|--tmp-dir|--threads|--memory|--kmer-size|--nanoplot-color|--lint-threshold|--lint-window|--lint-max-proportion|--genome-disagreement)
            require_option_value "$1" "${2:-}"
            case "$1" in
                --input-fastq) input_fastq="$2" ;;
                --sequencing-summary) sequencing_summary="$2" ;;
                --config) config_file="$2" ;;
                --output-dir) output_dir="$2" ;;
                --sample-name) sample_name="$2" ;;
                --tmp-dir) tmp_dir="$2" ;;
                --threads) threads="$2" ;;
                --memory) meryl_memory="$2" ;;
                --kmer-size) meryl_kmer_size="$2" ;;
                --nanoplot-color) nanoplot_color="$2" ;;
                --lint-threshold) lint_threshold="$2" ;;
                --lint-window) lint_window="$2" ;;
                --lint-max-proportion) lint_max_proportion="$2" ;;
                --genome-disagreement) genome_size_disagreement_threshold="$2" ;;
            esac
            shift 2
            ;;
        --skip-nanoplot) skip_nanoplot=true; shift ;;
        --skip-fastcat-stats) skip_fastcat_stats=true; shift ;;
        --skip-sdust) skip_sdust=true; shift ;;
        --skip-porechop-scan) skip_porechop_scan=true; shift ;;
        --skip-snikt) skip_snikt=true; shift ;;
        --skip-lrge) skip_lrge=true; shift ;;
        --skip-raven) skip_raven=true; shift ;;
        --skip-meryl) skip_meryl=true; shift ;;
        --input-sha256) input_sha256=true; shift ;;
        --force) force=true; shift ;;
        --dry-run) dry_run=true; shift ;;
        --help|-h) usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

[[ -n "${config_file}" ]] && load_config "${config_file}"

threads="${threads:-$(detect_threads)}"
meryl_kmer_size="${meryl_kmer_size:-21}"
nanoplot_color="${nanoplot_color:-green}"
lrge_seed="${lrge_seed:-123}"
lint_threshold="${lint_threshold:-20}"
lint_window="${lint_window:-64}"
lint_max_proportion="${lint_max_proportion:-0.95}"
genome_size_disagreement_threshold="${genome_size_disagreement_threshold:-0.20}"
input_sha256="${input_sha256:-false}"
force="${force:-false}"
skip_nanoplot="${skip_nanoplot:-false}"
skip_fastcat_stats="${skip_fastcat_stats:-false}"
skip_sdust="${skip_sdust:-false}"
skip_porechop_scan="${skip_porechop_scan:-false}"
skip_snikt="${skip_snikt:-false}"
skip_lrge="${skip_lrge:-false}"
skip_raven="${skip_raven:-false}"
skip_meryl="${skip_meryl:-false}"

# Validate every option before tool discovery or output creation.
require_arg "--input-fastq" "${input_fastq}"
require_file "${input_fastq}"
[[ -z "${sequencing_summary}" ]] || require_file "${sequencing_summary}"
require_positive_int "--threads" "${threads}"
require_positive_int "--kmer-size" "${meryl_kmer_size}"
require_positive_int "--lint-threshold" "${lint_threshold}"
require_positive_int "--lint-window" "${lint_window}"
require_nonnegative_int "lrge_seed" "${lrge_seed}"
require_number_range "--lint-max-proportion" "${lint_max_proportion}" 0 1 true true
require_number_range "--genome-disagreement" "${genome_size_disagreement_threshold}" 0 1 false true
for setting in input_sha256 force skip_nanoplot skip_fastcat_stats skip_sdust \
    skip_porechop_scan skip_snikt skip_lrge skip_raven skip_meryl; do
    require_boolean "${setting}" "${!setting}"
done

python_bin="$(stage_python)"
require_file "${stage_contract_py}"
input_fastq="$(canonical_path "${input_fastq}")"
if [[ -n "${sequencing_summary}" ]]; then
    sequencing_summary="$(canonical_path "${sequencing_summary}")"
fi
output_dir="$(canonical_path "${output_dir:-$(pwd)}")"
if [[ -z "${tmp_dir}" ]]; then
    tmp_dir="${output_dir}/.prokaryont_tmp"
elif [[ "${tmp_dir}" != /* ]]; then
    tmp_dir="${output_dir}/${tmp_dir}"
fi
tmp_dir="$(canonical_path "${tmp_dir}")"
sample_name="${sample_name:-$(infer_sample_name "${input_fastq}")}"
require_safe_sample_name "${sample_name}"

if [[ "${skip_meryl}" == "false" ]]; then
    if [[ -z "${meryl_memory}" ]]; then
        meryl_memory="$(detect_memory_gb || true)"
        [[ -n "${meryl_memory}" ]] || log_error "Could not detect an allocated memory limit for default-on Meryl; supply --memory or use --skip-meryl."
    fi
    require_positive_int "--memory" "${meryl_memory}"
fi

# Tool requirements follow enabled modules exactly.
[[ "${skip_nanoplot}" == "true" ]] || require_tool NanoPlot
if [[ "${skip_fastcat_stats}" == "false" || "${skip_sdust}" == "false" ]]; then
    require_tool fastcat
    fastcat_fastq_help="$(fastcat fastq --help 2>&1 || true)"
    grep -q -- '--force-error' <<< "${fastcat_fastq_help}" || log_error "Fastcat does not expose the required 1.x 'fastq' interface."
fi
if [[ "${skip_sdust}" == "false" ]]; then
    fastcat_lint_help="$(fastcat lint --help 2>&1 || true)"
    grep -q -- '--max-proportion' <<< "${fastcat_lint_help}" || log_error "Fastcat does not expose the required 1.x 'lint' interface."
fi
[[ "${skip_porechop_scan}" == "true" ]] || require_tool porechop_abi
[[ "${skip_snikt}" == "true" ]] || require_tool snikt.R
[[ "${skip_lrge}" == "true" ]] || require_tool lrge
[[ "${skip_raven}" == "true" ]] || require_tool raven
[[ "${skip_meryl}" == "true" ]] || require_tool meryl

qc_dir="${output_dir}/01_qc"
genome_size_dir="${output_dir}/02_genome_size"
manifest="${output_dir}/.prokaryont_stage1.complete.tsv"

IFS=$'\t' read -r _ input_size input_mtime_ns input_digest <<< "$(input_identity "${input_fastq}" "${input_sha256}")"
summary_size="not-supplied"
summary_mtime_ns="not-supplied"
summary_digest="not-supplied"
if [[ -n "${sequencing_summary}" ]]; then
    IFS=$'\t' read -r _ summary_size summary_mtime_ns summary_digest <<< "$(input_identity "${sequencing_summary}" "${input_sha256}")"
fi

fastcat_version="skipped"
nanoplot_version="skipped"
porechop_version="skipped"
snikt_version="skipped"
lrge_version="skipped"
raven_version="skipped"
meryl_version="skipped"
[[ "${skip_nanoplot}" == "true" ]] || nanoplot_version="$(tool_version NanoPlot)"
[[ "${skip_fastcat_stats}" == "true" && "${skip_sdust}" == "true" ]] || fastcat_version="$(tool_version fastcat)"
[[ "${skip_porechop_scan}" == "true" ]] || porechop_version="$(tool_version porechop_abi)"
[[ "${skip_snikt}" == "true" ]] || snikt_version="$(tool_version snikt.R)"
[[ "${skip_lrge}" == "true" ]] || lrge_version="$(tool_version lrge)"
[[ "${skip_raven}" == "true" ]] || raven_version="$(tool_version raven)"
[[ "${skip_meryl}" == "true" ]] || meryl_version="$(tool_version meryl)"
repository_revision="$(git -C "${_prokaryont_setup_dir}" rev-parse HEAD 2>/dev/null || printf 'unavailable')"

stage_signature="$({
    printf 'schema=prokaryont.stage1.signature.v1\n'
    printf 'input_path=%s\ninput_size=%s\ninput_mtime_ns=%s\ninput_sha256=%s\n' \
        "${input_fastq}" "${input_size}" "${input_mtime_ns}" "${input_digest}"
    printf 'summary_path=%s\nsummary_size=%s\nsummary_mtime_ns=%s\nsummary_sha256=%s\n' \
        "${sequencing_summary:-not-supplied}" "${summary_size}" "${summary_mtime_ns}" "${summary_digest}"
    printf 'sample_name=%s\nthreads=%s\nmeryl_memory=%s\nmeryl_kmer_size=%s\n' \
        "${sample_name}" "${threads}" "${meryl_memory:-skipped}" "${meryl_kmer_size}"
    printf 'nanoplot_color=%s\nlrge_seed=%s\nlint_threshold=%s\nlint_window=%s\n' \
        "${nanoplot_color}" "${lrge_seed}" "${lint_threshold}" "${lint_window}"
    printf 'lint_max_proportion=%s\ngenome_disagreement=%s\n' \
        "${lint_max_proportion}" "${genome_size_disagreement_threshold}"
    for setting in skip_nanoplot skip_fastcat_stats skip_sdust skip_porechop_scan \
        skip_snikt skip_lrge skip_raven skip_meryl; do
        printf '%s=%s\n' "${setting}" "${!setting}"
    done
    printf 'nanoplot_version=%s\nfastcat_version=%s\nporechop_version=%s\nsnikt_version=%s\n' \
        "${nanoplot_version}" "${fastcat_version}" "${porechop_version}" "${snikt_version}"
    printf 'lrge_version=%s\nraven_version=%s\nmeryl_version=%s\nrevision=%s\n' \
        "${lrge_version}" "${raven_version}" "${meryl_version}" "${repository_revision}"
} | stable_signature)"

_stage1_artifacts_present() {
    local path
    for path in \
        "${qc_dir}/run_summary_profile" "${qc_dir}/NanoPlot_sample" \
        "${qc_dir}/fastcat_output" "${qc_dir}/fastcat_histograms" \
        "${qc_dir}/fastcat_per_file_stats.tsv" "${qc_dir}/fastcat_sdust" \
        "${qc_dir}/porechop_abi_scan.log" "${qc_dir}/snikt_baseline" \
        "${genome_size_dir}/lrge_output.txt" "${genome_size_dir}/assembly.fasta" \
        "${genome_size_dir}/genome.meryl" "${genome_size_dir}/genome.hist" \
        "${genome_size_dir}/genome_size_estimates.tsv" "${genome_size_dir}/mean_genome_size.txt"; do
        [[ -e "${path}" ]] && return 0
    done
    return 1
}

if [[ -z "${dry_run:-}" ]]; then
    if [[ -f "${manifest}" && "$(manifest_signature "${manifest}")" == "${stage_signature}" ]] && \
        validate_completion_manifest "${manifest}"; then
        log_info "Stage 1 completion manifest and expected outputs match; nothing to do."
        exit 0
    fi
    if [[ "${force}" != "true" ]]; then
        if [[ -f "${manifest}" ]]; then
            log_error "Stage 1 completion manifest or outputs do not match this run. Re-run with --force to regenerate managed Stage 1 artifacts."
        elif _stage1_artifacts_present; then
            log_error "Stage 1 artifacts exist without a matching completion manifest. Re-run with --force to regenerate them safely."
        fi
    fi
fi

validate_fastq_nonempty "${input_fastq}" || log_error "Input FASTQ is empty, corrupt, truncated, or structurally invalid."

run_id="$(date -u +%Y%m%dT%H%M%SZ).$$"
stage_tmp=""
provenance_tsv=""
cleanup() {
    local exit_code=$?
    if [[ -n "${provenance_tsv:-}" && -f "${provenance_tsv}" ]]; then
        printf 'ended_utc\t%s\nexit_status\t%s\n' \
            "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${exit_code}" >> "${provenance_tsv}" || true
    fi
    [[ -n "${stage_tmp:-}" && -d "${stage_tmp}" ]] && rm -rf -- "${stage_tmp}" || true
    return "${exit_code}"
}
trap cleanup EXIT

if [[ -n "${dry_run:-}" ]]; then
    stage_tmp="${tmp_dir}/stage1.DRY_RUN"
else
    mkdir -p "${output_dir}/provenance" "${qc_dir}" "${genome_size_dir}" "${tmp_dir}"
    stage_tmp="$(mktemp -d "${tmp_dir}/stage1.XXXXXX")"
    init_stage_logging "${output_dir}/provenance/${sample_name}.stage1.${run_id}.log"
    provenance_tsv="${output_dir}/provenance/${sample_name}.stage1.${run_id}.tsv"
    {
        printf 'key\tvalue\n'
        printf 'schema\tprokaryont.stage-provenance.v1\nrun_id\t%s\nstage\t1\nsample_name\t%s\n' "${run_id}" "${sample_name}"
        printf 'started_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'input_path\t%s\ninput_size\t%s\ninput_mtime_ns\t%s\ninput_sha256\t%s\n' "${input_fastq}" "${input_size}" "${input_mtime_ns}" "${input_digest}"
        printf 'sequencing_summary\t%s\noutput_dir\t%s\ntmp_dir\t%s\n' "${sequencing_summary:-not-supplied}" "${output_dir}" "${tmp_dir}"
        printf 'signature\t%s\nrepository_revision\t%s\nthreads\t%s\nmeryl_memory_gb\t%s\n' "${stage_signature}" "${repository_revision}" "${threads}" "${meryl_memory:-skipped}"
        printf 'meryl_kmer_size\t%s\nnanoplot_color\t%s\nlrge_seed\t%s\n' "${meryl_kmer_size}" "${nanoplot_color}" "${lrge_seed}"
        printf 'lint_threshold\t%s\nlint_window\t%s\nlint_max_proportion\t%s\ngenome_disagreement_threshold\t%s\n' \
            "${lint_threshold}" "${lint_window}" "${lint_max_proportion}" "${genome_size_disagreement_threshold}"
        for setting in skip_nanoplot skip_fastcat_stats skip_sdust skip_porechop_scan \
            skip_snikt skip_lrge skip_raven skip_meryl; do
            printf 'option.%s\t%s\n' "${setting}" "${!setting}"
        done
        printf 'slurm_job_id\t%s\nslurm_cpus_per_task\t%s\n' "${SLURM_JOB_ID:-not-set}" "${SLURM_CPUS_PER_TASK:-not-set}"
        for tool in NanoPlot fastcat porechop_abi snikt.R lrge raven meryl; do
            printf 'tool_path.%s\t%s\n' "${tool}" "$(tool_path "${tool}")"
        done
        printf 'tool_version.NanoPlot\t%s\ntool_version.fastcat\t%s\ntool_version.porechop_abi\t%s\n' "${nanoplot_version}" "${fastcat_version}" "${porechop_version}"
        printf 'tool_version.snikt\t%s\ntool_version.lrge\t%s\ntool_version.raven\t%s\ntool_version.meryl\t%s\n' "${snikt_version}" "${lrge_version}" "${raven_version}" "${meryl_version}"
    } > "${provenance_tsv}"
fi

stage_qc="${stage_tmp}/01_qc"
stage_genome="${stage_tmp}/02_genome_size"
run_cmd mkdir -p "${stage_qc}" "${stage_genome}"

log_step "Stage 1: non-destructive read QC"

if [[ "${skip_nanoplot}" == "false" ]]; then
    if [[ -n "${sequencing_summary}" ]]; then
        run_cmd NanoPlot --threads "${threads}" --summary "${sequencing_summary}" \
            --loglength -o "${stage_qc}/run_summary_profile"
    fi
    run_cmd NanoPlot --threads "${threads}" -c "${nanoplot_color}" \
        --fastq "${input_fastq}" -o "${stage_qc}/NanoPlot_sample" \
        --loglength --plots dot --N50
fi

if [[ "${skip_fastcat_stats}" == "false" ]]; then
    fastcat_tmp_output="${stage_tmp}/fastcat-output"
    run_cmd bash -c 'fastcat fastq --force-error --output "$1" "$2" > /dev/null' \
        _ "${fastcat_tmp_output}" "${input_fastq}"
    if [[ -z "${dry_run:-}" ]]; then
        [[ -s "${fastcat_tmp_output}/per_file.txt" ]] || log_error "Fastcat did not produce a non-empty per_file.txt."
        [[ -s "${fastcat_tmp_output}/histograms/length.hist" ]] || log_error "Fastcat did not produce length.hist."
        [[ -s "${fastcat_tmp_output}/histograms/quality.hist" ]] || log_error "Fastcat did not produce quality.hist."
        mv -- "${fastcat_tmp_output}" "${stage_qc}/fastcat_output"
        cp -- "${stage_qc}/fastcat_output/per_file.txt" "${stage_qc}/fastcat_per_file_stats.tsv"
        cp -R -- "${stage_qc}/fastcat_output/histograms" "${stage_qc}/fastcat_histograms"
    fi
fi

if [[ "${skip_sdust}" == "false" ]]; then
    sdust_dir="${stage_qc}/fastcat_sdust"
    run_cmd mkdir -p "${sdust_dir}"
    run_cmd bash -c 'set -o pipefail; fastcat lint --threshold "$1" --window "$2" --max-proportion 0.0 "$3" 2> "$4" | awk "END { print int(NR / 4) }" > "$5"' \
        _ "${lint_threshold}" "${lint_window}" "${input_fastq}" \
        "${sdust_dir}/lint_p0.stderr.log" "${sdust_dir}/zero_sdust_reads.txt"
    run_cmd bash -c 'fastcat lint --threshold "$1" --window "$2" --max-proportion "$3" "$4" > /dev/null 2> "$5"' \
        _ "${lint_threshold}" "${lint_window}" "${lint_max_proportion}" \
        "${input_fastq}" "${sdust_dir}/lint_preview.stderr.log"

    if [[ -z "${dry_run:-}" ]]; then
        {
            printf 'read_id\treported_masked_fraction_rounded\n'
            sed -nE 's/.*Read ([^ ]+) masked fraction ([0-9.]+) exceeds threshold [0-9.]+, skipping\..*/\1\t\2/p' \
                "${sdust_dir}/lint_p0.stderr.log"
        } > "${sdust_dir}/sdust_positive_reads.tsv"
        zero_sdust_reads="$(awk 'NR == 1 { print $1; exit }' "${sdust_dir}/zero_sdust_reads.txt")"
        zero_sdust_reads="${zero_sdust_reads:-0}"
        sdust_positive_reads="$(awk 'NR > 1 { count++ } END { print count + 0 }' "${sdust_dir}/sdust_positive_reads.tsv")"
        skipped_log_lines="$(awk '/skipping\./ { count++ } END { print count + 0 }' "${sdust_dir}/lint_p0.stderr.log")"
        [[ "${skipped_log_lines}" -eq "${sdust_positive_reads}" ]] || \
            log_error "Could not parse every Fastcat lint rejection; no SDUST summary was published."
        preview_rejected="$(awk '/skipping\./ { count++ } END { print count + 0 }' "${sdust_dir}/lint_preview.stderr.log")"
        total_sdust_reads=$((zero_sdust_reads + sdust_positive_reads))
        zero_fraction="$(awk -v z="${zero_sdust_reads}" -v t="${total_sdust_reads}" 'BEGIN { printf "%.6f", t ? z/t : 0 }')"
        positive_fraction="$(awk -v p="${sdust_positive_reads}" -v t="${total_sdust_reads}" 'BEGIN { printf "%.6f", t ? p/t : 0 }')"
        reported_gt_050="$(awk -F '\t' 'NR > 1 && $2 > 0.50 { n++ } END { print n + 0 }' "${sdust_dir}/sdust_positive_reads.tsv")"
        reported_gt_080="$(awk -F '\t' 'NR > 1 && $2 > 0.80 { n++ } END { print n + 0 }' "${sdust_dir}/sdust_positive_reads.tsv")"
        reported_gt_090="$(awk -F '\t' 'NR > 1 && $2 > 0.90 { n++ } END { print n + 0 }' "${sdust_dir}/sdust_positive_reads.tsv")"
        reported_gt_095="$(awk -F '\t' 'NR > 1 && $2 > 0.95 { n++ } END { print n + 0 }' "${sdust_dir}/sdust_positive_reads.tsv")"
        reported_max="$(awk -F '\t' 'NR > 1 && $2 > max { max=$2 } END { printf "%.2f", max + 0 }' "${sdust_dir}/sdust_positive_reads.tsv")"
        {
            printf 'class\treported_masked_fraction_rounded\tread_count\n'
            printf 'exact_zero\t0.00\t%s\n' "${zero_sdust_reads}"
            awk -F '\t' 'NR > 1 { count[$2]++ } END { for (fraction in count) print "positive_reported", fraction, count[fraction] }' OFS='\t' \
                "${sdust_dir}/sdust_positive_reads.tsv" | sort -t $'\t' -k2,2n
        } > "${sdust_dir}/sdust_fraction_hist.tsv"
        {
            head -n 1 "${sdust_dir}/sdust_positive_reads.tsv"
            awk -F '\t' 'NR > 1 && $2 >= 0.80' "${sdust_dir}/sdust_positive_reads.tsv" | sort -t $'\t' -k2,2nr
        } > "${sdust_dir}/sdust_high_burden_reads.tsv"
        {
            printf 'metric\tvalue\n'
            printf 'sdust_threshold\t%s\nsdust_window\t%s\ntotal_reads\t%s\n' "${lint_threshold}" "${lint_window}" "${total_sdust_reads}"
            printf 'zero_sdust_reads\t%s\nsdust_positive_reads\t%s\n' "${zero_sdust_reads}" "${sdust_positive_reads}"
            printf 'zero_sdust_fraction\t%s\nsdust_positive_fraction\t%s\n' "${zero_fraction}" "${positive_fraction}"
            printf 'reported_fraction_gt_0.50\t%s\nreported_fraction_gt_0.80\t%s\n' "${reported_gt_050}" "${reported_gt_080}"
            printf 'reported_fraction_gt_0.90\t%s\nreported_fraction_gt_0.95\t%s\n' "${reported_gt_090}" "${reported_gt_095}"
            printf 'reported_fraction_max\t%s\npreview_max_proportion\t%s\npreview_rejected_reads_exact_threshold\t%s\n' \
                "${reported_max}" "${lint_max_proportion}" "${preview_rejected}"
        } > "${sdust_dir}/sdust_summary.tsv"
    fi
fi

if [[ "${skip_porechop_scan}" == "false" ]]; then
    run_cmd bash -c 'porechop_abi -abi --guess_adapter_only -i "$1" --threads "$2" > "$3" 2>&1' \
        _ "${input_fastq}" "${threads}" "${stage_qc}/porechop_abi_scan.log"
fi

if [[ "${skip_snikt}" == "false" ]]; then
    run_cmd mkdir -p "${stage_qc}/snikt_baseline"
    run_cmd bash -c 'cd "$1" && snikt.R --notrim "$2"' \
        _ "${stage_qc}/snikt_baseline" "${input_fastq}"
fi

log_step "Stage 1: independent genome-size diagnostics"
lrge_size=""
raven_size=""
if [[ "${skip_lrge}" == "false" ]]; then
    run_cmd bash -c 'lrge -P ont -t "$1" -s "$2" "$3" > "$4"' \
        _ "${threads}" "${lrge_seed}" "${input_fastq}" "${stage_genome}/lrge_output.txt"
    if [[ -z "${dry_run:-}" ]]; then
        lrge_size="$(awk 'NF { print $1; exit }' "${stage_genome}/lrge_output.txt")"
        [[ "${lrge_size}" =~ ^[1-9][0-9]*$ ]] || log_error "LRGE did not produce a positive integer estimate."
    fi
fi

if [[ "${skip_raven}" == "false" ]]; then
    run_cmd bash -c 'raven --threads "$1" -p 0 --graphical-fragment-assembly "$2" "$3" > "$4"' \
        _ "${threads}" "${stage_genome}/assemblyGraph.gfa" "${input_fastq}" "${stage_genome}/assembly.fasta"
    if [[ -z "${dry_run:-}" ]]; then
        raven_size="$(awk '!/^>/ { gsub(/[[:space:]]/, ""); n += length($0) } END { print n + 0 }' "${stage_genome}/assembly.fasta")"
        [[ "${raven_size}" =~ ^[1-9][0-9]*$ ]] || log_error "Raven did not produce a non-empty assembly span."
        rm -f -- "${stage_genome}/assemblyGraph.gfa"
    fi
fi

if [[ "${skip_meryl}" == "false" ]]; then
    run_cmd meryl count compress "k=${meryl_kmer_size}" "memory=${meryl_memory}" "threads=${threads}" \
        output "${stage_genome}/genome.meryl" "${input_fastq}"
    run_cmd bash -c 'meryl histogram "$1" > "$2"' \
        _ "${stage_genome}/genome.meryl" "${stage_genome}/genome.hist"
fi

if [[ -z "${dry_run:-}" ]]; then
    absolute_difference="NA"
    relative_difference="NA"
    disagreement_status="not_available"
    if [[ -n "${lrge_size}" && -n "${raven_size}" ]]; then
        absolute_difference=$((lrge_size > raven_size ? lrge_size - raven_size : raven_size - lrge_size))
        relative_difference="$(awk -v a="${lrge_size}" -v b="${raven_size}" 'BEGIN { printf "%.6f", (a == 0 && b == 0) ? 0 : (a > b ? a-b : b-a) / ((a+b)/2) }')"
        if awk -v value="${relative_difference}" -v threshold="${genome_size_disagreement_threshold}" \
            'BEGIN { exit value >= threshold ? 0 : 1 }'; then
            disagreement_status="warning"
            log_warn "LRGE and Raven differ by ${relative_difference} (threshold ${genome_size_disagreement_threshold}); treat both as diagnostics."
        else
            disagreement_status="within_threshold"
        fi
    fi
    {
        printf 'metric\tvalue\tunit\ttool_version\n'
        printf 'lrge_estimate\t%s\tbp\t%s\n' "${lrge_size:-NA}" "${lrge_version}"
        printf 'raven_assembly_span\t%s\tbp\t%s\n' "${raven_size:-NA}" "${raven_version}"
        printf 'absolute_difference\t%s\tbp\tNA\n' "${absolute_difference}"
        printf 'relative_difference\t%s\tfraction\tNA\n' "${relative_difference}"
        printf 'disagreement_threshold\t%s\tfraction\tNA\n' "${genome_size_disagreement_threshold}"
        printf 'disagreement_status\t%s\tcategory\tNA\n' "${disagreement_status}"
    } > "${stage_genome}/genome_size_estimates.tsv"
fi

if [[ -n "${dry_run:-}" ]]; then
    log_info "[DRY-RUN] No directories, logs, temporary files, manifests, or results were created."
    exit 0
fi

# All producers succeeded. Replace only Stage 1-managed artifacts, then write
# the completion marker last. Stage 2 QC content under 01_qc is untouched.
managed_paths=(
    "${qc_dir}/run_summary_profile" "${qc_dir}/NanoPlot_sample"
    "${qc_dir}/fastcat_output" "${qc_dir}/fastcat_histograms"
    "${qc_dir}/fastcat_per_file_stats.tsv" "${qc_dir}/fastcat_sdust"
    "${qc_dir}/porechop_abi_scan.log" "${qc_dir}/snikt_baseline"
    "${genome_size_dir}/lrge_output.txt" "${genome_size_dir}/assembly.fasta"
    "${genome_size_dir}/assemblyGraph.gfa" "${genome_size_dir}/genome.meryl"
    "${genome_size_dir}/genome.hist" "${genome_size_dir}/genome_size_estimates.tsv"
    "${genome_size_dir}/mean_genome_size.txt"
)
rm -rf -- "${managed_paths[@]}"

_publish_if_present() {
    local source_path="$1"
    local destination_path="$2"
    [[ -e "${source_path}" ]] || return 0
    mkdir -p "$(dirname "${destination_path}")"
    mv -- "${source_path}" "${destination_path}"
}

for name in run_summary_profile NanoPlot_sample fastcat_output fastcat_histograms \
    fastcat_per_file_stats.tsv fastcat_sdust porechop_abi_scan.log snikt_baseline; do
    _publish_if_present "${stage_qc}/${name}" "${qc_dir}/${name}"
done
for name in lrge_output.txt assembly.fasta genome.meryl genome.hist genome_size_estimates.tsv; do
    _publish_if_present "${stage_genome}/${name}" "${genome_size_dir}/${name}"
done

expected=("nonempty:${genome_size_dir}/genome_size_estimates.tsv" "nonempty:${provenance_tsv}" "nonempty:${log_file}")
[[ "${skip_nanoplot}" == "true" ]] || expected+=("nonempty:${qc_dir}/NanoPlot_sample/NanoPlot-report.html")
if [[ "${skip_nanoplot}" == "false" && -n "${sequencing_summary}" ]]; then
    expected+=("nonempty:${qc_dir}/run_summary_profile/NanoPlot-report.html")
fi
if [[ "${skip_fastcat_stats}" == "false" ]]; then
    expected+=("nonempty:${qc_dir}/fastcat_per_file_stats.tsv" "nonempty:${qc_dir}/fastcat_histograms/length.hist" "nonempty:${qc_dir}/fastcat_histograms/quality.hist")
fi
if [[ "${skip_sdust}" == "false" ]]; then
    expected+=("nonempty:${qc_dir}/fastcat_sdust/sdust_summary.tsv" "nonempty:${qc_dir}/fastcat_sdust/sdust_fraction_hist.tsv" "nonempty:${qc_dir}/fastcat_sdust/sdust_positive_reads.tsv" "nonempty:${qc_dir}/fastcat_sdust/sdust_high_burden_reads.tsv")
fi
[[ "${skip_porechop_scan}" == "true" ]] || expected+=("file:${qc_dir}/porechop_abi_scan.log")
[[ "${skip_snikt}" == "true" ]] || expected+=("dir:${qc_dir}/snikt_baseline")
[[ "${skip_lrge}" == "true" ]] || expected+=("nonempty:${genome_size_dir}/lrge_output.txt")
[[ "${skip_raven}" == "true" ]] || expected+=("nonempty:${genome_size_dir}/assembly.fasta")
[[ "${skip_meryl}" == "true" ]] || expected+=("dir:${genome_size_dir}/genome.meryl" "nonempty:${genome_size_dir}/genome.hist")

write_completion_manifest "${manifest}" 1 "${stage_signature}" "${run_id}" \
    "field:sample_name=${sample_name}" "field:input_path=${input_fastq}" \
    "field:input_size=${input_size}" "field:input_mtime_ns=${input_mtime_ns}" \
    "field:input_sha256=${input_digest}" "field:output_dir=${output_dir}" \
    "field:provenance=${provenance_tsv}" "${expected[@]}"
validate_completion_manifest "${manifest}" || log_error "Stage 1 completion manifest validation failed after publication."

log_step "Stage 1 complete"
log_info "QC outputs: ${qc_dir}"
log_info "Genome-size diagnostics: ${genome_size_dir}/genome_size_estimates.tsv"
log_info "Completion manifest: ${manifest}"
