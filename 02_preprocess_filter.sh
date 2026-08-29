#!/usr/bin/env bash
# Stage 2: optional trimming plus quality/length filtering and final QC.

source "$(dirname "$0")/00_setup.sh"
log_file=""

threads="${threads:-}"
input_fastq="${input_fastq:-}"
output_dir="${output_dir:-}"
output_fastq="${output_fastq:-}"
sample_name="${sample_name:-}"
tmp_dir="${tmp_dir:-}"
config_file="${config_file:-}"
filtlong_min_length="${filtlong_min_length:-}"
min_qscore="${min_qscore:-}"
filtlong_keep_percent="${filtlong_keep_percent:-}"
enable_porechopabi_trim="${enable_porechopabi_trim:-}"
enable_chopper_trim="${enable_chopper_trim:-}"
chopper_trim_approach="${chopper_trim_approach:-}"
chopper_headcrop="${chopper_headcrop:-}"
chopper_tailcrop="${chopper_tailcrop:-}"
chopper_trim_cutoff="${chopper_trim_cutoff:-}"
chopper_split_window="${chopper_split_window:-}"
enable_fastcat_lint="${enable_fastcat_lint:-}"
lint_threshold="${lint_threshold:-}"
lint_window="${lint_window:-}"
lint_max_proportion="${lint_max_proportion:-}"
genome_size_override="${genome_size_override:-}"
input_sha256="${input_sha256:-}"
force="${force:-}"
skip_snikt="${skip_snikt:-}"
skip_nanoplot="${skip_nanoplot:-}"
skip_fastplong="${skip_fastplong:-}"
skip_intermediate_metrics="${skip_intermediate_metrics:-}"

usage() {
    cat <<'EOF'
Usage: 02_preprocess_filter.sh --input-fastq FILE [OPTIONS]

Required:
  --input-fastq FILE             Raw FASTQ (plain or gzip-compressed)

Layout and identity:
  --config FILE                  Path to pipeline.conf
  --output-dir DIR               Output root (default: current directory)
  --output-fastq FILE            Final gzip FASTQ (default: OUTPUT/filtered_input.fastq.gz)
                                Relative paths are resolved under --output-dir
  --sample-name NAME             Manifest/sample label (default: input basename)
  --tmp-dir DIR                  Large-intermediate directory (default: OUTPUT/.prokaryont_tmp)
  --input-sha256                 Include full input SHA-256 in the completion signature
  --force                        Regenerate after a missing/mismatched completion manifest

Filtering:
  --threads N                    Threads (default: SLURM allocation, then nproc)
  --min-length N                 Minimum read length (default: 200)
  --min-qscore N                 Probability-space mean read Q threshold (default: 7)
  --keep-percent N               Optional Filtlong ranking; must be >0 and <100
  --enable-porechopabi-trim      Enable Porechop_ABI adapter trimming
  --enable-chopper-trim          Enable Chopper end trimming
  --chopper-trim-approach MODE   fixed-crop | trim-by-quality | best-read-segment | split-by-low-quality
  --headcrop N                   Fixed-crop bases from the read start (default: 50)
  --tailcrop N                   Fixed-crop bases from the read end (default: 30)
  --chopper-trim-cutoff N        Required Q cutoff for quality-aware Chopper modes
  --split-window N               Consecutive low-Q bases required to split (default: 1)
  --enable-fastcat-lint          Enable explicit SDUST read filtering
  --lint-threshold N             SDUST threshold T (default: 20)
  --lint-window N                SDUST window W (default: 64)
  --lint-max-proportion N        Maximum masked fraction to retain (default: 0.95)

Reporting:
  --genome-size N                Authoritative genome size for coverage reporting only
  --skip-snikt                   Skip post-Porechop SNIKT diagnostic
  --skip-nanoplot                Skip final NanoPlot diagnostic
  --skip-fastplong               Skip final Fastplong diagnostic
  --skip-intermediate-metrics    Report input/final metrics only (reduces scans)
  --dry-run                      Validate and print commands without filesystem mutation
  -h, --help                     Show this help
EOF
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-fastq|--config|--output-dir|--output-fastq|--sample-name|--tmp-dir|--threads|--min-length|--min-qscore|--keep-percent|--chopper-trim-approach|--headcrop|--tailcrop|--chopper-trim-cutoff|--split-window|--lint-threshold|--lint-window|--lint-max-proportion|--genome-size)
            require_option_value "$1" "${2:-}"
            case "$1" in
                --input-fastq) input_fastq="$2" ;;
                --config) config_file="$2" ;;
                --output-dir) output_dir="$2" ;;
                --output-fastq) output_fastq="$2" ;;
                --sample-name) sample_name="$2" ;;
                --tmp-dir) tmp_dir="$2" ;;
                --threads) threads="$2" ;;
                --min-length) filtlong_min_length="$2" ;;
                --min-qscore) min_qscore="$2" ;;
                --keep-percent) filtlong_keep_percent="$2" ;;
                --chopper-trim-approach) chopper_trim_approach="$2" ;;
                --headcrop) chopper_headcrop="$2" ;;
                --tailcrop) chopper_tailcrop="$2" ;;
                --chopper-trim-cutoff) chopper_trim_cutoff="$2" ;;
                --split-window) chopper_split_window="$2" ;;
                --lint-threshold) lint_threshold="$2" ;;
                --lint-window) lint_window="$2" ;;
                --lint-max-proportion) lint_max_proportion="$2" ;;
                --genome-size) genome_size_override="$2" ;;
            esac
            shift 2
            ;;
        --enable-porechopabi-trim) enable_porechopabi_trim=true; shift ;;
        --enable-chopper-trim) enable_chopper_trim=true; shift ;;
        --enable-fastcat-lint) enable_fastcat_lint=true; shift ;;
        --skip-snikt) skip_snikt=true; shift ;;
        --skip-nanoplot) skip_nanoplot=true; shift ;;
        --skip-fastplong) skip_fastplong=true; shift ;;
        --skip-intermediate-metrics) skip_intermediate_metrics=true; shift ;;
        --input-sha256) input_sha256=true; shift ;;
        --force) force=true; shift ;;
        --dry-run) dry_run=true; shift ;;
        --help|-h) usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

[[ -n "${config_file}" ]] && load_config "${config_file}"

threads="${threads:-$(detect_threads)}"
filtlong_min_length="${filtlong_min_length:-200}"
min_qscore="${min_qscore:-7}"
enable_porechopabi_trim="${enable_porechopabi_trim:-false}"
enable_chopper_trim="${enable_chopper_trim:-false}"
chopper_trim_approach="${chopper_trim_approach:-fixed-crop}"
chopper_headcrop="${chopper_headcrop:-50}"
chopper_tailcrop="${chopper_tailcrop:-30}"
chopper_split_window="${chopper_split_window:-1}"
enable_fastcat_lint="${enable_fastcat_lint:-false}"
lint_threshold="${lint_threshold:-20}"
lint_window="${lint_window:-64}"
lint_max_proportion="${lint_max_proportion:-0.95}"
input_sha256="${input_sha256:-false}"
force="${force:-false}"
skip_snikt="${skip_snikt:-false}"
skip_nanoplot="${skip_nanoplot:-false}"
skip_fastplong="${skip_fastplong:-false}"
skip_intermediate_metrics="${skip_intermediate_metrics:-false}"

# Complete validation occurs before any output is created or any producer runs.
require_arg "--input-fastq" "${input_fastq}"
require_file "${input_fastq}"
require_positive_int "--threads" "${threads}"
require_positive_int "--min-length" "${filtlong_min_length}"
require_number_range "--min-qscore" "${min_qscore}" 0 93 true true
require_positive_int "--lint-threshold" "${lint_threshold}"
require_positive_int "--lint-window" "${lint_window}"
require_number_range "--lint-max-proportion" "${lint_max_proportion}" 0 1 true true
require_nonnegative_int "--headcrop" "${chopper_headcrop}"
require_nonnegative_int "--tailcrop" "${chopper_tailcrop}"
require_positive_int "--split-window" "${chopper_split_window}"
if [[ -n "${filtlong_keep_percent}" ]]; then
    require_number_range "--keep-percent" "${filtlong_keep_percent}" 0 100 false false
fi
if [[ -n "${genome_size_override}" ]]; then
    require_positive_int "--genome-size" "${genome_size_override}"
fi
for setting in enable_porechopabi_trim enable_chopper_trim enable_fastcat_lint \
    input_sha256 force skip_snikt skip_nanoplot skip_fastplong skip_intermediate_metrics; do
    require_boolean "${setting}" "${!setting}"
done
case "${chopper_trim_approach}" in
    fixed-crop)
        [[ -z "${chopper_trim_cutoff}" ]] || log_error "--chopper-trim-cutoff is not used by fixed-crop."
        ;;
    trim-by-quality|best-read-segment|split-by-low-quality)
        if [[ "${enable_chopper_trim}" == "true" ]]; then
            require_arg "--chopper-trim-cutoff" "${chopper_trim_cutoff}"
        fi
        [[ -z "${chopper_trim_cutoff}" ]] || require_number_range "--chopper-trim-cutoff" "${chopper_trim_cutoff}" 0 93 true true
        if [[ "${chopper_trim_approach}" != "split-by-low-quality" && "${chopper_split_window}" != "1" ]]; then
            log_error "--split-window is only meaningful with split-by-low-quality."
        fi
        ;;
    *) log_error "Unknown Chopper trimming mode: ${chopper_trim_approach}" ;;
esac

python_bin="$(stage_python)"
require_file "${stage_contract_py}"
input_fastq="$(canonical_path "${input_fastq}")"
output_dir="$(canonical_path "${output_dir:-$(pwd)}")"
if [[ -z "${output_fastq}" ]]; then
    output_fastq="${output_dir}/filtered_input.fastq.gz"
elif [[ "${output_fastq}" != /* ]]; then
    output_fastq="${output_dir}/${output_fastq}"
fi
output_fastq="$(canonical_path "${output_fastq}")"
if [[ -z "${tmp_dir}" ]]; then
    tmp_dir="${output_dir}/.prokaryont_tmp"
elif [[ "${tmp_dir}" != /* ]]; then
    tmp_dir="${output_dir}/${tmp_dir}"
fi
tmp_dir="$(canonical_path "${tmp_dir}")"
sample_name="${sample_name:-$(infer_sample_name "${input_fastq}")}"
require_safe_sample_name "${sample_name}"

if [[ "${input_fastq}" == "${output_fastq}" ]] || { [[ -e "${output_fastq}" ]] && [[ "${input_fastq}" -ef "${output_fastq}" ]]; }; then
    log_error "Input and output FASTQ resolve to the same file; refusing to overwrite the input."
fi

require_tool seqkit
require_tool gzip
[[ -z "${filtlong_keep_percent}" ]] || require_tool filtlong
[[ "${enable_porechopabi_trim}" == "false" ]] || require_tool porechop_abi
if [[ "${enable_porechopabi_trim}" == "true" && "${skip_snikt}" == "false" ]]; then
    require_tool snikt.R
fi
if [[ "${enable_chopper_trim}" == "true" ]]; then
    require_tool chopper
    chopper_help="$(chopper --help 2>&1 || true)"
    grep -q -- '--trim-approach' <<< "${chopper_help}" || log_error "Chopper lacks --trim-approach; install Chopper 0.13.x."
    grep -q -- '--input' <<< "${chopper_help}" || log_error "Chopper lacks direct --input support."
    if [[ "${chopper_trim_approach}" == "split-by-low-quality" ]]; then
        grep -q -- '--split-window' <<< "${chopper_help}" || log_error "Chopper lacks --split-window; install Chopper 0.13.x."
    fi
fi
if [[ "${enable_fastcat_lint}" == "true" ]]; then
    require_tool fastcat
    fastcat_lint_help="$(fastcat lint --help 2>&1 || true)"
    grep -q -- '--max-proportion' <<< "${fastcat_lint_help}" || log_error "Fastcat does not expose the required 1.x 'lint' interface."
fi
[[ "${skip_nanoplot}" == "true" ]] || require_tool NanoPlot
[[ "${skip_fastplong}" == "true" ]] || require_tool fastplong

qc_dir="${output_dir}/01_qc"
metrics_path="${output_dir}/02_filtering_metrics.tsv"
manifest="${output_dir}/.prokaryont_stage2.complete.tsv"

IFS=$'\t' read -r _ input_size input_mtime_ns input_digest <<< "$(input_identity "${input_fastq}" "${input_sha256}")"
seqkit_version="$(tool_version seqkit)"
gzip_version="$(tool_version gzip)"
compressor_name="gzip"
compressor_version="${gzip_version}"
if command -v pigz &>/dev/null; then
    compressor_name="pigz"
    compressor_version="$(tool_version pigz)"
fi
filtlong_version="skipped"
porechop_version="skipped"
snikt_version="skipped"
chopper_version="skipped"
fastcat_version="skipped"
nanoplot_version="skipped"
fastplong_version="skipped"
[[ -z "${filtlong_keep_percent}" ]] || filtlong_version="$(tool_version filtlong)"
[[ "${enable_porechopabi_trim}" == "false" ]] || porechop_version="$(tool_version porechop_abi)"
[[ "${enable_porechopabi_trim}" == "false" || "${skip_snikt}" == "true" ]] || snikt_version="$(tool_version snikt.R)"
[[ "${enable_chopper_trim}" == "false" ]] || chopper_version="$(tool_version chopper)"
[[ "${enable_fastcat_lint}" == "false" ]] || fastcat_version="$(tool_version fastcat)"
[[ "${skip_nanoplot}" == "true" ]] || nanoplot_version="$(tool_version NanoPlot)"
[[ "${skip_fastplong}" == "true" ]] || fastplong_version="$(tool_version fastplong)"
repository_revision="$(git -C "${_prokaryont_setup_dir}" rev-parse HEAD 2>/dev/null || printf 'unavailable')"

stage_signature="$({
    printf 'schema=prokaryont.stage2.signature.v1\n'
    printf 'input_path=%s\ninput_size=%s\ninput_mtime_ns=%s\ninput_sha256=%s\n' "${input_fastq}" "${input_size}" "${input_mtime_ns}" "${input_digest}"
    printf 'output_fastq=%s\nsample_name=%s\nthreads=%s\nmin_length=%s\nmin_qscore=%s\n' "${output_fastq}" "${sample_name}" "${threads}" "${filtlong_min_length}" "${min_qscore}"
    printf 'keep_percent=%s\nporechop=%s\nchopper=%s\nchopper_mode=%s\nheadcrop=%s\ntailcrop=%s\ncutoff=%s\nsplit_window=%s\n' \
        "${filtlong_keep_percent:-disabled}" "${enable_porechopabi_trim}" "${enable_chopper_trim}" "${chopper_trim_approach}" \
        "${chopper_headcrop}" "${chopper_tailcrop}" "${chopper_trim_cutoff:-not-set}" "${chopper_split_window}"
    printf 'fastcat_lint=%s\nlint_threshold=%s\nlint_window=%s\nlint_max=%s\ngenome_size=%s\n' \
        "${enable_fastcat_lint}" "${lint_threshold}" "${lint_window}" "${lint_max_proportion}" "${genome_size_override:-not-set}"
    printf 'skip_snikt=%s\nskip_nanoplot=%s\nskip_fastplong=%s\nskip_intermediate_metrics=%s\n' \
        "${skip_snikt}" "${skip_nanoplot}" "${skip_fastplong}" "${skip_intermediate_metrics}"
    printf 'seqkit_version=%s\nfiltlong_version=%s\nporechop_version=%s\nsnikt_version=%s\n' "${seqkit_version}" "${filtlong_version}" "${porechop_version}" "${snikt_version}"
    printf 'chopper_version=%s\nfastcat_version=%s\nnanoplot_version=%s\nfastplong_version=%s\ncompressor=%s\ncompressor_version=%s\nrevision=%s\n' \
        "${chopper_version}" "${fastcat_version}" "${nanoplot_version}" "${fastplong_version}" "${compressor_name}" "${compressor_version}" "${repository_revision}"
} | stable_signature)"

_stage2_artifacts_present() {
    local path
    for path in "${output_fastq}" "${metrics_path}" "${qc_dir}/NanoPlot_FiltPol" \
        "${qc_dir}/fastplong_report.html" "${qc_dir}/fastplong_report.json" \
        "${qc_dir}/snikt_post_porechop" "${qc_dir}/02_lint_rejected.log"; do
        [[ -e "${path}" ]] && return 0
    done
    return 1
}

if [[ -z "${dry_run:-}" ]]; then
    if [[ -f "${manifest}" && "$(manifest_signature "${manifest}")" == "${stage_signature}" ]] && \
        validate_completion_manifest "${manifest}"; then
        log_info "Stage 2 completion manifest and expected outputs match; nothing to do."
        exit 0
    fi
    if [[ "${force}" != "true" ]]; then
        if [[ -f "${manifest}" ]]; then
            log_error "Stage 2 completion manifest or outputs do not match this run. Re-run with --force to regenerate safely."
        elif _stage2_artifacts_present; then
            log_error "Stage 2 artifacts exist without a matching completion manifest. Re-run with --force to replace the managed outputs atomically."
        fi
    fi
fi

validate_fastq_nonempty "${input_fastq}" || log_error "Input FASTQ is empty, corrupt, truncated, or structurally invalid."

run_id="$(date -u +%Y%m%dT%H%M%SZ).$$"
stage_tmp=""
final_tmp=""
provenance_tsv=""
cleanup() {
    local exit_code=$?
    if [[ -n "${provenance_tsv:-}" && -f "${provenance_tsv}" ]]; then
        printf 'ended_utc\t%s\nexit_status\t%s\n' \
            "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${exit_code}" >> "${provenance_tsv}" || true
    fi
    [[ -n "${stage_tmp:-}" && -d "${stage_tmp}" ]] && rm -rf -- "${stage_tmp}" || true
    [[ -n "${final_tmp:-}" && -f "${final_tmp}" ]] && rm -f -- "${final_tmp}" || true
    return "${exit_code}"
}
trap cleanup EXIT

if [[ -n "${dry_run:-}" ]]; then
    stage_tmp="${tmp_dir}/stage2.DRY_RUN"
    final_tmp="$(dirname "${output_fastq}")/.${output_fastq##*/}.DRY_RUN"
else
    mkdir -p "${output_dir}/provenance" "${qc_dir}" "${tmp_dir}" "$(dirname "${output_fastq}")"
    stage_tmp="$(mktemp -d "${tmp_dir}/stage2.XXXXXX")"
    final_tmp="$(mktemp "$(dirname "${output_fastq}")/.${output_fastq##*/}.tmp.XXXXXX.fastq.gz")"
    init_stage_logging "${output_dir}/provenance/${sample_name}.stage2.${run_id}.log"
    provenance_tsv="${output_dir}/provenance/${sample_name}.stage2.${run_id}.tsv"
    {
        printf 'key\tvalue\n'
        printf 'schema\tprokaryont.stage-provenance.v1\nrun_id\t%s\nstage\t2\nsample_name\t%s\n' "${run_id}" "${sample_name}"
        printf 'started_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'input_path\t%s\ninput_size\t%s\ninput_mtime_ns\t%s\ninput_sha256\t%s\n' "${input_fastq}" "${input_size}" "${input_mtime_ns}" "${input_digest}"
        printf 'output_dir\t%s\noutput_fastq\t%s\ntmp_dir\t%s\nsignature\t%s\n' "${output_dir}" "${output_fastq}" "${tmp_dir}" "${stage_signature}"
        printf 'repository_revision\t%s\nthreads\t%s\nslurm_job_id\t%s\nslurm_cpus_per_task\t%s\n' "${repository_revision}" "${threads}" "${SLURM_JOB_ID:-not-set}" "${SLURM_CPUS_PER_TASK:-not-set}"
        printf 'min_length\t%s\nmin_qscore\t%s\nkeep_percent\t%s\ngenome_size\t%s\n' \
            "${filtlong_min_length}" "${min_qscore}" "${filtlong_keep_percent:-disabled}" "${genome_size_override:-not-set}"
        printf 'enable_porechopabi_trim\t%s\nenable_chopper_trim\t%s\nchopper_trim_approach\t%s\n' \
            "${enable_porechopabi_trim}" "${enable_chopper_trim}" "${chopper_trim_approach}"
        printf 'chopper_headcrop\t%s\nchopper_tailcrop\t%s\nchopper_trim_cutoff\t%s\nchopper_split_window\t%s\n' \
            "${chopper_headcrop}" "${chopper_tailcrop}" "${chopper_trim_cutoff:-not-set}" "${chopper_split_window}"
        printf 'enable_fastcat_lint\t%s\nlint_threshold\t%s\nlint_window\t%s\nlint_max_proportion\t%s\n' \
            "${enable_fastcat_lint}" "${lint_threshold}" "${lint_window}" "${lint_max_proportion}"
        printf 'skip_snikt\t%s\nskip_nanoplot\t%s\nskip_fastplong\t%s\nskip_intermediate_metrics\t%s\n' \
            "${skip_snikt}" "${skip_nanoplot}" "${skip_fastplong}" "${skip_intermediate_metrics}"
        for tool in seqkit filtlong porechop_abi snikt.R chopper fastcat NanoPlot fastplong gzip; do
            printf 'tool_path.%s\t%s\n' "${tool}" "$(tool_path "${tool}")"
        done
        printf 'tool_version.seqkit\t%s\ntool_version.filtlong\t%s\ntool_version.porechop_abi\t%s\n' "${seqkit_version}" "${filtlong_version}" "${porechop_version}"
        printf 'tool_version.snikt\t%s\ntool_version.chopper\t%s\ntool_version.fastcat\t%s\n' "${snikt_version}" "${chopper_version}" "${fastcat_version}"
        printf 'tool_version.NanoPlot\t%s\ntool_version.fastplong\t%s\ntool_version.gzip\t%s\n' "${nanoplot_version}" "${fastplong_version}" "${gzip_version}"
        printf 'compressor\t%s\ncompressor_version\t%s\n' "${compressor_name}" "${compressor_version}"
    } > "${provenance_tsv}"
fi

stage_qc="${stage_tmp}/01_qc"
metrics_tmp="${stage_tmp}/02_filtering_metrics.tsv"
run_cmd mkdir -p "${stage_qc}"

_read_metrics() {
    local reads_file="$1"
    seqkit stats -a -T "${reads_file}" 2>/dev/null | awk -F '\t' '
        NR == 1 {
            for (i=1; i<=NF; i++) column[$i]=i
            next
        }
        NR == 2 {
            reads = column["num_seqs"] ? $(column["num_seqs"]) : "NA"
            bases = column["sum_len"] ? $(column["sum_len"]) : "NA"
            minlen = column["min_len"] ? $(column["min_len"]) : "NA"
            avglen = column["avg_len"] ? $(column["avg_len"]) : "NA"
            median = column["Q2"] ? $(column["Q2"]) : "NA"
            n50 = column["N50"] ? $(column["N50"]) : "NA"
            avgqual = column["AvgQual"] ? $(column["AvgQual"]) : "NA"
            print reads, bases, minlen, avglen, median, n50, avgqual
        }' OFS='\t'
}

initial_reads=""
initial_bases=""
previous_reads=""
previous_bases=""
_append_metrics() {
    local step="$1"
    local reads_file="$2"
    local reads bases minlen avglen median n50 avgqual
    IFS=$'\t' read -r reads bases minlen avglen median n50 avgqual <<< "$(_read_metrics "${reads_file}")"
    [[ "${reads}" =~ ^[0-9]+$ && "${bases}" =~ ^[0-9]+$ ]] || log_error "SeqKit could not report reads/bases for ${step}."
    if [[ -z "${initial_reads}" ]]; then
        initial_reads="${reads}"
        initial_bases="${bases}"
        previous_reads="${reads}"
        previous_bases="${bases}"
    fi
    local read_retention base_retention read_change base_change coverage="NA"
    read_retention="$(awk -v current="${reads}" -v initial="${initial_reads}" 'BEGIN { printf "%.3f", initial ? 100*current/initial : 0 }')"
    base_retention="$(awk -v current="${bases}" -v initial="${initial_bases}" 'BEGIN { printf "%.3f", initial ? 100*current/initial : 0 }')"
    read_change=$((reads - previous_reads))
    base_change=$((bases - previous_bases))
    if [[ -n "${genome_size_override}" ]]; then
        coverage="$(awk -v bases="${bases}" -v genome="${genome_size_override}" 'BEGIN { printf "%.3f", bases/genome }')"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${step}" "${reads}" "${bases}" "${read_change}" "${base_change}" \
        "${read_retention}" "${base_retention}" "${minlen}" "${avglen}" \
        "${median}" "${n50}" "${avgqual}" "${coverage}" >> "${metrics_tmp}"
    previous_reads="${reads}"
    previous_bases="${bases}"
}

if [[ -z "${dry_run:-}" ]]; then
    printf 'step\treads\tbases\tread_change_from_previous\tbase_change_from_previous\tread_retention_pct\tbase_retention_pct\tmin_length\tmean_length\tmedian_length\tn50\tprobability_mean_q\testimated_coverage\n' > "${metrics_tmp}"
    _append_metrics input "${input_fastq}"
fi

log_step "Stage 2: preprocessing and filtering"
preproc_input="${input_fastq}"

if [[ "${enable_porechopabi_trim}" == "true" ]]; then
    porechop_tmp="${stage_tmp}/porechop.fastq"
    run_cmd porechop_abi -abi -i "${preproc_input}" -o "${porechop_tmp}" --threads "${threads}"
    preproc_input="${porechop_tmp}"
    [[ -n "${dry_run:-}" || "${skip_intermediate_metrics}" == "true" ]] || _append_metrics post_porechop "${preproc_input}"
    if [[ "${skip_snikt}" == "false" ]]; then
        run_cmd mkdir -p "${stage_qc}/snikt_post_porechop"
        run_cmd bash -c 'cd "$1" && snikt.R --notrim "$2"' \
            _ "${stage_qc}/snikt_post_porechop" "${preproc_input}"
    fi
fi

if [[ "${enable_chopper_trim}" == "true" ]]; then
    chopper_tmp="${stage_tmp}/chopper.fastq"
    chopper_args=(--trim-approach "${chopper_trim_approach}" --threads "${threads}")
    case "${chopper_trim_approach}" in
        fixed-crop) chopper_args+=(--headcrop "${chopper_headcrop}" --tailcrop "${chopper_tailcrop}") ;;
        trim-by-quality|best-read-segment) chopper_args+=(--cutoff "${chopper_trim_cutoff}") ;;
        split-by-low-quality)
            chopper_args+=(--cutoff "${chopper_trim_cutoff}" --split-window "${chopper_split_window}" --minlength "${filtlong_min_length}")
            ;;
    esac
    chopper_args+=(--input "${preproc_input}")
    run_cmd bash -c 'output=$1; shift; chopper "$@" > "$output"' \
        _ "${chopper_tmp}" "${chopper_args[@]}"
    preproc_input="${chopper_tmp}"
    [[ -n "${dry_run:-}" || "${skip_intermediate_metrics}" == "true" ]] || _append_metrics post_chopper "${preproc_input}"
fi

if [[ "${enable_fastcat_lint}" == "true" ]]; then
    lint_tmp="${stage_tmp}/fastcat_lint.fastq"
    run_cmd bash -c 'fastcat lint --threshold "$1" --window "$2" --max-proportion "$3" "$4" 2> "$5" > "$6"' \
        _ "${lint_threshold}" "${lint_window}" "${lint_max_proportion}" \
        "${preproc_input}" "${stage_qc}/02_lint_rejected.log" "${lint_tmp}"
    preproc_input="${lint_tmp}"
    if [[ -z "${dry_run:-}" ]]; then
        [[ "${skip_intermediate_metrics}" == "true" ]] || _append_metrics post_fastcat_lint "${preproc_input}"
        rejected_count="$(awk '/skipping\./ { n++ } END { print n + 0 }' "${stage_qc}/02_lint_rejected.log")"
        log_info "Fastcat lint rejected ${rejected_count} reads at the exact ${lint_max_proportion} threshold."
    fi
fi

seqkit_tmp="${stage_tmp}/seqkit_filtered.fastq"
run_cmd seqkit seq --min-qual "${min_qscore}" --min-len "${filtlong_min_length}" \
    "${preproc_input}" -o "${seqkit_tmp}"
[[ -n "${dry_run:-}" || "${skip_intermediate_metrics}" == "true" ]] || _append_metrics post_seqkit "${seqkit_tmp}"

if [[ -n "${filtlong_keep_percent}" ]]; then
    filtlong_args=(--keep_percent "${filtlong_keep_percent}" "${seqkit_tmp}")
    if command -v pigz &>/dev/null; then
        run_cmd bash -c 'set -o pipefail; output=$1; threads=$2; shift 2; filtlong "$@" | pigz -p "$threads" > "$output"' \
            _ "${final_tmp}" "${threads}" "${filtlong_args[@]}"
    else
        run_cmd bash -c 'set -o pipefail; output=$1; shift; filtlong "$@" | gzip > "$output"' \
            _ "${final_tmp}" "${filtlong_args[@]}"
    fi
else
    if command -v pigz &>/dev/null; then
        run_cmd bash -c 'pigz -p "$1" -c "$2" > "$3"' _ "${threads}" "${seqkit_tmp}" "${final_tmp}"
    else
        run_cmd bash -c 'gzip -c "$1" > "$2"' _ "${seqkit_tmp}" "${final_tmp}"
    fi
fi

if [[ -n "${dry_run:-}" ]]; then
    if [[ "${skip_nanoplot}" == "false" ]]; then
        run_cmd NanoPlot --threads "${threads}" -c green --fastq "${output_fastq}" \
            -o "${stage_qc}/NanoPlot_FiltPol" --loglength --plots dot --N50
    fi
    if [[ "${skip_fastplong}" == "false" ]]; then
        run_cmd fastplong -i "${output_fastq}" -Q -L -A \
            --html "${stage_qc}/fastplong_report.html" \
            --json "${stage_qc}/fastplong_report.json" --thread "${threads}"
    fi
    log_info "[DRY-RUN] No directories, logs, temporary files, manifests, or results were created."
    exit 0
fi

validate_fastq_nonempty "${final_tmp}" || log_error "Filtered output failed gzip/FASTQ/nonzero validation; the previous final output was not replaced."
_append_metrics final "${final_tmp}"

if [[ "${skip_nanoplot}" == "false" ]]; then
    run_cmd NanoPlot --threads "${threads}" -c green --fastq "${final_tmp}" \
        -o "${stage_qc}/NanoPlot_FiltPol" --loglength --plots dot --N50
fi
if [[ "${skip_fastplong}" == "false" ]]; then
    # Fastplong remains report-only. SeqKit owns the user-facing
    # probability-space mean-Q filter, while Filtlong is used only for optional
    # length/window-quality ranking through --keep_percent.
    run_cmd fastplong -i "${final_tmp}" -Q -L -A \
        --html "${stage_qc}/fastplong_report.html" \
        --json "${stage_qc}/fastplong_report.json" --thread "${threads}"
fi

# Publish all Stage 2 artifacts only after every producer succeeds. The final
# FASTQ temporary is on the destination filesystem, so this replacement is
# atomic and an older result remains untouched until this point.
mv -f -- "${final_tmp}" "${output_fastq}"
final_tmp=""
rm -rf -- "${qc_dir}/NanoPlot_FiltPol" "${qc_dir}/fastplong_report.html" \
    "${qc_dir}/fastplong_report.json" "${qc_dir}/snikt_post_porechop" \
    "${qc_dir}/02_lint_rejected.log" "${metrics_path}"
for name in NanoPlot_FiltPol fastplong_report.html fastplong_report.json \
    snikt_post_porechop 02_lint_rejected.log; do
    [[ -e "${stage_qc}/${name}" ]] && mv -- "${stage_qc}/${name}" "${qc_dir}/${name}"
done
mv -- "${metrics_tmp}" "${metrics_path}"

expected=("fastq:${output_fastq}" "nonempty:${metrics_path}" "nonempty:${provenance_tsv}" "nonempty:${log_file}")
[[ "${skip_nanoplot}" == "true" ]] || expected+=("nonempty:${qc_dir}/NanoPlot_FiltPol/NanoPlot-report.html")
if [[ "${skip_fastplong}" == "false" ]]; then
    expected+=("nonempty:${qc_dir}/fastplong_report.html" "nonempty:${qc_dir}/fastplong_report.json")
fi
if [[ "${enable_porechopabi_trim}" == "true" && "${skip_snikt}" == "false" ]]; then
    expected+=("dir:${qc_dir}/snikt_post_porechop")
fi
[[ "${enable_fastcat_lint}" == "false" ]] || expected+=("file:${qc_dir}/02_lint_rejected.log")

write_completion_manifest "${manifest}" 2 "${stage_signature}" "${run_id}" \
    "field:sample_name=${sample_name}" "field:input_path=${input_fastq}" \
    "field:input_size=${input_size}" "field:input_mtime_ns=${input_mtime_ns}" \
    "field:input_sha256=${input_digest}" "field:output_fastq=${output_fastq}" \
    "field:provenance=${provenance_tsv}" "${expected[@]}"
validate_completion_manifest "${manifest}" || log_error "Stage 2 completion manifest validation failed after publication."

log_step "Stage 2 complete"
log_info "Filtered reads: ${output_fastq}"
log_info "Comprehensive filtering metrics: ${metrics_path}"
log_info "Completion manifest: ${manifest}"
