#!/usr/bin/env bash
# ==============================================================================
# 03_autocycler_assemble.sh — Autocycler multi-assembler subsampling & resolution
# ==============================================================================
# Usage:
#   bash 03_autocycler_assemble.sh --read-type ont_r10 [OPTIONS]
#   (reads filtered_input.fastq.gz by default, or custom reads path via --reads / --input-fastq)
#
# Optional flags:
#   --config FILE           Path to pipeline.conf (values override defaults)
#   --threads N             Number of threads (default: 128)
#   --sample-name NAME      Sample name for metrics (default: MyBacteria)
#   --read-type TYPE        Read type: ont_r9|ont_r10|pacbio_clr|pacbio_hifi (default: ont_r10)
#   --genome-size SIZE      Override genome size (skip re-estimation)
#   --assemblers LIST       Comma-separated Autocycler helper tasks to run
#   --subsample-count N     Number of read subsamples (default: 4)
#   --subsampler NAME       autocycler|rasusa (default: autocycler)
#   --subsample-seed N      Base seed for reproducible subsampling (default: 0)
#   --rasusa-coverage N     Optional Rasusa target coverage per subset
#   --min-read-depth N      Min subset read depth for subsampling (default: 25)
#   --combine-reads PATH    Reads used for final combine depth annotation
#   --combine-depth-kmer N  Combine depth k-mer, odd 11-31 (default: 19)
#   --combine-threads N     Combine threads, 1-100 (default: min(threads, 100))
#   --compress-threads N    Compress threads, 1-100 (default: min(threads, 100))
#   --plassembler-db DIR    Plassembler database used for the PLSDB screen
#   --skip-plsdb-screen     Skip post-consensus PLSDB characterization
#   --parallel-jobs N       Concurrent assembler jobs (default: 4)
#   --canu-parallel-jobs N  Concurrent Canu jobs (default: 2)
#   --trim-mad N            Autocycler trim MAD outlier threshold (default: 5)
#   --trim-min-identity N   Autocycler trim min alignment identity (default: 0.75)
#   --skip-curation         Skip manual curation pauses
#   --dry-run               Print commands without executing
#   --help                  Show this help
#
# To pass extra arguments to individual assemblers, set them in pipeline.conf
# (e.g. flye_extra_args=--meta) or export them before running:
#   export assembler_args=( [flye]="--meta" [raven]="--graphical-fragment-assembly" )
# ==============================================================================

source "$(dirname "$0")/00_setup.sh"

# --- Defaults ----------------------------------------------------------------
threads="${threads:-}"
read_type="${read_type:-}"
sample_name="${sample_name:-}"
genome_size_override="${genome_size_override:-}"
subsample_count="${subsample_count:-}"
subsampler="${subsampler:-}"
subsample_seed="${subsample_seed:-}"
rasusa_coverage="${rasusa_coverage:-}"
min_read_depth="${min_read_depth:-}"
min_depth_rel="${min_depth_rel:-}"
canu_extra_args="${canu_extra_args:-}"
parallel_jobs="${parallel_jobs:-}"
canu_parallel_jobs="${canu_parallel_jobs:-}"
skip_curation="${skip_curation:-}"
trim_mad="${trim_mad:-}"
trim_min_identity="${trim_min_identity:-}"
config_file="${config_file:-}"
reads_path="${reads_path:-}"
autocycler_combine_reads="${autocycler_combine_reads:-}"
autocycler_combine_depth_kmer="${autocycler_combine_depth_kmer:-}"
autocycler_combine_threads="${autocycler_combine_threads:-}"
autocycler_compress_threads="${autocycler_compress_threads:-}"
plassembler_db="${plassembler_db:-}"
skip_plsdb_screen="${skip_plsdb_screen:-}"

# Assembler extra args (associative array — set via config file or before sourcing)
declare -A assembler_args 2>/dev/null || true
: "${assembler_args[flye]:=}"
: "${assembler_args[metamdbg]:=}"
: "${assembler_args[miniasm]:=}"
: "${assembler_args[raven]:=}"
: "${assembler_args[plassembler]:=}"
: "${assembler_args[hifiasm]:=}"
: "${assembler_args[myloasm]:=}"
: "${assembler_args[nextdenovo]:=}"
: "${assembler_args[ilesta]:=}"
: "${assembler_args[lja]:=}"
: "${assembler_args[redbean]:=}"
: "${assembler_args[necat]:=}"

assembler_list="flye,canu,hifiasm,ilesta,lja,raven,miniasm,metamdbg,myloasm,plassembler,nextdenovo,redbean,necat"
assemblers_explicit=false
if [[ -n "${PROKARYONT_ASSEMBLERS:-}" ]]; then
    assembler_list="${PROKARYONT_ASSEMBLERS}"
    assemblers_explicit=true
fi

# --- Usage -------------------------------------------------------------------
usage() {
    echo "Usage: $(basename "$0") --read-type TYPE [OPTIONS]"
    echo ""
    echo "Optional:"
    echo "  --config FILE               Path to pipeline.conf"
    echo "  --reads, --input-fastq PATH Path to reads (default: filtered_input.fastq.gz)"
    echo "  --read-type TYPE            Read type: ont_r9|ont_r10|pacbio_clr|pacbio_hifi (default: ont_r10)"
    echo "  --threads N                 Number of threads (default: 128)"
    echo "  --sample-name NAME          Sample name for metrics (default: MyBacteria)"
    echo "  --genome-size SIZE          Override genome size (skip re-estimation)"
    echo "  --assemblers LIST           Comma-separated Autocycler helper tasks"
    echo "                              Canonical: flye,canu,hifiasm,ilesta,lja,raven,miniasm,metamdbg,myloasm,plassembler,nextdenovo,redbean,necat"
    echo "                              Compatibility alias: wtdbg2 (normalized to redbean)"
    echo "  --subsample-count N         Number of read subsamples (default: 4)"
    echo "  --subsampler NAME           autocycler|rasusa (default: autocycler)"
    echo "  --subsample-seed N          Base seed for reproducible subsampling (default: 0)"
    echo "  --rasusa-coverage FLOAT     Optional Rasusa target coverage per subset"
    echo "  --min-read-depth N          Min subset read depth for subsampling (default: 25)"
    echo "  --min-depth-rel FLOAT       Min relative depth for accepted assemblies (default: 0.1)"
    echo "  --parallel-jobs N           Concurrent assembler jobs (default: 4)"
    echo "  --canu-parallel-jobs N      Concurrent Canu jobs (default: 2)"
    echo "  --trim-mad N                Autocycler trim MAD outlier threshold (default: 5)"
    echo "  --trim-min-identity N       Autocycler trim min alignment identity (default: 0.75)"
    echo "  --combine-reads PATH        Reads for final combine depth annotation (default: Stage 3 input)"
    echo "  --combine-depth-kmer N      Combine depth k-mer, odd 11-31 (default: 19)"
    echo "  --combine-threads N         Combine threads, 1-100 (default: min(global threads, 100))"
    echo "  --compress-threads N        Compress threads, 1-100 (default: min(global threads, 100))"
    echo "  --plassembler-db DIR        Plassembler database for the PLSDB screen"
    echo "  --skip-plsdb-screen         Skip post-consensus PLSDB characterization"
    echo "  --skip-curation             Skip manual curation pauses"
    echo "  --dry-run                   Print commands without executing"
    echo "  --help                      Show this help"
    exit 0
}

# --- Parse arguments ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --reads|--input-fastq) reads_path="$2"; shift 2 ;;
        --read-type)          read_type="$2"; shift 2 ;;
        --config)             config_file="$2"; shift 2 ;;
        --threads)            threads="$2"; shift 2 ;;
        --sample-name)        sample_name="$2"; shift 2 ;;
        --genome-size)        genome_size_override="$2"; shift 2 ;;
        --assemblers)         assembler_list="$2"; assemblers_explicit=true; shift 2 ;;
        --subsample-count)    subsample_count="$2"; shift 2 ;;
        --subsampler)         subsampler="$2"; shift 2 ;;
        --subsample-seed)     subsample_seed="$2"; shift 2 ;;
        --rasusa-coverage)    rasusa_coverage="$2"; shift 2 ;;
        --min-read-depth)     min_read_depth="$2"; shift 2 ;;
        --min-depth-rel)      min_depth_rel="$2"; shift 2 ;;
        --parallel-jobs)      parallel_jobs="$2"; shift 2 ;;
        --canu-parallel-jobs) canu_parallel_jobs="$2"; shift 2 ;;
        --trim-mad)           trim_mad="$2"; shift 2 ;;
        --trim-min-identity)  trim_min_identity="$2"; shift 2 ;;
        --combine-reads)      autocycler_combine_reads="$2"; shift 2 ;;
        --combine-depth-kmer) autocycler_combine_depth_kmer="$2"; shift 2 ;;
        --combine-threads)    autocycler_combine_threads="$2"; shift 2 ;;
        --compress-threads)   autocycler_compress_threads="$2"; shift 2 ;;
        --plassembler-db)     plassembler_db="$2"; shift 2 ;;
        --skip-plsdb-screen)  skip_plsdb_screen=true; shift ;;
        --skip-curation)      skip_curation=true; shift ;;
        --dry-run)            dry_run=true; shift ;;
        --help|-h)            usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

# --- Load config (CLI flags parsed above take priority) ----------------------
[[ -n "${config_file}" ]] && load_config "${config_file}"

if [[ "${assemblers_explicit}" != "true" && -n "${assemblers:-}" ]]; then
    assembler_list="${assemblers}"
fi
threads="${threads:-128}"
read_type="${read_type:-ont_r10}"
sample_name="${sample_name:-MyBacteria}"
subsample_count="${subsample_count:-4}"
min_depth_rel="${min_depth_rel:-0.1}"
parallel_jobs="${parallel_jobs:-4}"
canu_parallel_jobs="${canu_parallel_jobs:-2}"

subsampler="${subsampler:-autocycler}"
subsample_seed="${subsample_seed:-0}"
autocycler_combine_depth_kmer="${autocycler_combine_depth_kmer:-19}"
skip_plsdb_screen="${skip_plsdb_screen:-false}"

# Populate assembler_args from config variables if not already set
for _asm in flye metamdbg miniasm raven plassembler hifiasm myloasm nextdenovo ilesta lja redbean necat; do
    _var="${_asm}_extra_args"
    if [[ -n "${!_var:-}" && -z "${assembler_args[$_asm]:-}" ]]; then
        assembler_args[$_asm]="${!_var}"
    fi
done
if [[ -z "${assembler_args[redbean]:-}" && -n "${wtdbg2_extra_args:-}" ]]; then
    assembler_args[redbean]="${wtdbg2_extra_args}"
fi

normalize_assembler() {
    case "$1" in
        wtdbg2) printf '%s\n' redbean ;;
        *)      printf '%s\n' "$1" ;;
    esac
}

declare -a assemblers=()
declare -A _seen_assemblers=()
IFS=',' read -r -a _requested_assemblers <<< "${assembler_list}"
for _requested in "${_requested_assemblers[@]}"; do
    _requested="${_requested//[[:space:]]/}"
    [[ -z "${_requested}" ]] && continue
    _normalized=$(normalize_assembler "${_requested}")
    if [[ "${_requested}" == "wtdbg2" ]]; then
        log_info "Assembler alias normalized: wtdbg2 -> redbean"
    fi
    if [[ -z "${_seen_assemblers[${_normalized}]:-}" ]]; then
        assemblers+=("${_normalized}")
        _seen_assemblers["${_normalized}"]=1
    fi
done
[[ ${#assemblers[@]} -gt 0 ]] || log_error "No assemblers were selected."

# --- Validate read type early (before any expensive work) --------------------
if [[ ! "${read_type}" =~ ^(ont_r9|ont_r10|pacbio_clr|pacbio_hifi)$ ]]; then
    log_error "Invalid read type '${read_type}'. Valid: ont_r9, ont_r10, pacbio_clr, pacbio_hifi."
fi
for _integer_setting in threads subsample_count parallel_jobs canu_parallel_jobs; do
    _integer_value="${!_integer_setting}"
    if [[ ! "${_integer_value}" =~ ^[1-9][0-9]*$ ]]; then
        log_error "${_integer_setting} must be a positive integer."
    fi
done
if [[ ! "${subsampler}" =~ ^(autocycler|rasusa)$ ]]; then
    log_error "Invalid subsampler '${subsampler}'. Valid: autocycler, rasusa."
fi
if [[ ! "${subsample_seed}" =~ ^[0-9]+$ ]]; then
    log_error "--subsample-seed must be a non-negative integer."
fi
if [[ -n "${rasusa_coverage}" ]] && { [[ ! "${rasusa_coverage}" =~ ^[0-9]+([.][0-9]+)?$ ]] || ! awk -v v="${rasusa_coverage}" 'BEGIN { exit (v > 0) ? 0 : 1 }'; }; then
    log_error "--rasusa-coverage must be a positive number."
fi
if [[ ! "${autocycler_combine_depth_kmer}" =~ ^[0-9]+$ ]] ||
   (( autocycler_combine_depth_kmer < 11 || autocycler_combine_depth_kmer > 31 || autocycler_combine_depth_kmer % 2 == 0 )); then
    log_error "--combine-depth-kmer must be an odd integer from 11 to 31."
fi
if [[ -n "${autocycler_combine_threads}" ]] && { [[ ! "${autocycler_combine_threads}" =~ ^[0-9]+$ ]] || (( autocycler_combine_threads < 1 || autocycler_combine_threads > 100 )); }; then
    log_error "--combine-threads must be an integer from 1 to 100."
fi
if [[ -n "${autocycler_compress_threads}" ]] && { [[ ! "${autocycler_compress_threads}" =~ ^[0-9]+$ ]] || (( autocycler_compress_threads < 1 || autocycler_compress_threads > 100 )); }; then
    log_error "--compress-threads must be an integer from 1 to 100."
fi

# Set GZIP_BIN now that threads are parsed
export GZIP_BIN="$(get_gzip_cmd)"
# --- Validate ----------------------------------------------------------------
require_tool autocycler
require_tool parallel
require_tool seqkit
require_tool minimap2
require_tool samtools

[[ "${subsampler}" == "rasusa" ]] && require_tool rasusa
if [[ "${subsampler}" == "rasusa" ]] && ! rasusa reads --help &>/dev/null; then
    log_error "Installed Rasusa does not expose the modern 'rasusa reads' command. Update Rasusa or use --subsampler autocycler."
fi

for assembler in "${assemblers[@]}"; do
    if ! autocycler helper "${assembler}" --help &>/dev/null; then
        log_error "Installed Autocycler does not expose helper task '${assembler}'. Update Autocycler or remove ${assembler} from the assembler list."
    fi
done

declare -a selected_assembler_tools=()
for assembler in "${assemblers[@]}"; do
    case "${assembler}" in
        metamdbg) selected_assembler_tools+=(metaMDBG) ;;
        redbean)  selected_assembler_tools+=(wtdbg2) ;;
        ilesta)   selected_assembler_tools+=(Ilesta minipolish minimap2 racon) ;;
        lja)      selected_assembler_tools+=(lja) ;;
        miniasm)  selected_assembler_tools+=(miniasm minipolish minimap2 racon) ;;
        *)        selected_assembler_tools+=("${assembler}") ;;
    esac
done
for tool in "${selected_assembler_tools[@]}"; do
    command -v "${tool}" &>/dev/null || log_warn "'${tool}' not found; selected assembler jobs may fail."
done

if [[ "${skip_plsdb_screen}" != "true" ]]; then
    require_tool plassembler
    if [[ -z "${plassembler_db}" && -n "${PLASSEMBLER_DB:-}" && -d "${PLASSEMBLER_DB}" ]]; then
        plassembler_db="${PLASSEMBLER_DB}"
    elif [[ -z "${plassembler_db}" && -n "${CONDA_PREFIX:-}" && -d "${CONDA_PREFIX}/plassembler_db" ]]; then
        plassembler_db="${CONDA_PREFIX}/plassembler_db"
    fi
    if [[ -z "${plassembler_db}" || ! -d "${plassembler_db}" ]]; then
        log_error "PLSDB screening is enabled, but no Plassembler database was found. Set --plassembler-db, plassembler_db, or PLASSEMBLER_DB; or use --skip-plsdb-screen."
    fi
fi

for tool in nucmer mummerplot; do
    command -v "${tool}" &>/dev/null || log_warn "'${tool}' not found; chromosome-scale dotplots will be skipped."
done

# --- Derived values ----------------------------------------------------------
threads_per_job=$(( threads / parallel_jobs ))
canu_threads_per_job=$(( threads / canu_parallel_jobs ))
qc_dir="$(pwd)/01_qc"
genome_size_dir="$(pwd)/02_genome_size"
assemblies_dir="$(pwd)/assemblies"
autocycler_dir="$(pwd)/autocycler_out"
filtered_reads="$(pwd)/filtered_input.fastq.gz"
input_reads="${reads_path:-${filtered_reads}}"
combine_reads="${autocycler_combine_reads:-${input_reads}}"
autocycler_consensus="${autocycler_dir}/consensus_assembly.fasta"

if [[ -n "${reads_path:-}" ]]; then
    require_file "${input_reads}" "(user-supplied --reads path)"
else
    require_file "${input_reads}" "02_preprocess_filter.sh"
fi
require_file "${combine_reads}" "Provide --combine-reads PATH or a valid Stage 3 input FASTQ"

if [[ -n "${autocycler_combine_threads}" ]]; then
    combine_threads="${autocycler_combine_threads}"
else
    combine_threads="${threads}"
    if (( combine_threads > 100 )); then
        combine_threads=100
        log_info "Autocycler combine supports at most 100 threads; using 100 of requested ${threads}."
    fi
fi

if [[ -n "${autocycler_compress_threads}" ]]; then
    compress_threads="${autocycler_compress_threads}"
else
    compress_threads="${threads}"
    if (( compress_threads > 100 )); then
        compress_threads=100
        log_info "Autocycler compress supports at most 100 threads; using 100 of requested ${threads}."
    fi
fi

# --- Utility: portable float comparison (replaces bc -l dependency) ----------
_float_gt() {
    awk -v v1="$1" -v v2="$2" 'BEGIN{exit (v1+0 > v2+0) ? 0 : 1}' 2>/dev/null
}

_normalize_genome_size_bp() {
    local value="${1,,}"
    local multiplier=1
    case "${value}" in
        *k) multiplier=1000; value="${value%?}" ;;
        *m) multiplier=1000000; value="${value%?}" ;;
        *g) multiplier=1000000000; value="${value%?}" ;;
    esac
    [[ "${value}" =~ ^[0-9]+([.][0-9]+)?$ ]] || return 1
    awk -v value="${value}" -v multiplier="${multiplier}" \
        'BEGIN { printf "%.0f\n", value * multiplier }'
}

# ==============================================================================
# STEP 9 — Autocycler Assembly
# ==============================================================================

log_step "Step 9: Autocycler assembly pipeline"

# --- Resume logic: skip completed assembly sub-stages (cluster/trim/dnaapler) -
_resume_target="${PROKARYONT_RESUME:-}"
if [[ "${_resume_target}" == "cluster" || "${_resume_target}" == "trim" || "${_resume_target}" == "dnaapler" ]]; then
    log_info ">>> RESUME from '${_resume_target}' requested — skipping completed assembly sub-stages"
fi

# --- 9a. Genome size recovery (consolidated — previously duplicated between
#     the old resume-branch and the old STEP 4a; collapses into one block
#     now that filtering always already happened in a separate script) ------
if [[ -n "${genome_size_override}" ]]; then
    genome_size="${genome_size_override}"
    log_info "9a. Using user-provided genome size: ${genome_size}"
elif [[ -f "${genome_size_dir}/mean_genome_size.txt" ]]; then
    genome_size=$(cat "${genome_size_dir}/mean_genome_size.txt")
    log_info "9a. Using weighted mean genome size from QC step: ${genome_size}"
else
    log_info "9a. Estimating genome size from input reads..."
    if [[ -n "${dry_run:-}" ]]; then
        log_info "[DRY-RUN] autocycler helper genome_size --reads ${input_reads} --threads ${threads}"
        genome_size="DRY_RUN_PLACEHOLDER"
    else
        genome_size=$(autocycler helper genome_size \
            --reads "${input_reads}" --threads "${threads}")
    fi
    log_info "Autocycler genome size: ${genome_size}"
fi
if [[ -f "${genome_size_dir}/lrge_output.txt" ]]; then
    log_info "LRGE estimate: $(cat "${genome_size_dir}/lrge_output.txt")"
fi

# Pre-flight: measure available depth for the selected subsampling strategy.
_genome_size_bp=$(_normalize_genome_size_bp "${genome_size}" 2>/dev/null || echo 0)
_avail_bases=$(seqkit stats -T "${input_reads}" 2>/dev/null \
    | awk -F'\t' 'NR==2{print $5}' || echo 0)
_avail_bases="${_avail_bases:-0}"
if [[ "${subsampler}" == "autocycler" && "${_genome_size_bp}" =~ ^[1-9][0-9]*$ ]]; then
    _min_required=$(awk -v gs="${_genome_size_bp}" -v sc="${subsample_count}" \
        -v md="${min_read_depth:-25}" 'BEGIN{printf "%.0f", gs * sc * md}')
    if awk -v av="${_avail_bases}" -v req="${_min_required}" \
        'BEGIN{exit (av+0 < req+0) ? 0 : 1}' 2>/dev/null; then
        log_warn "Available bases (${_avail_bases}) may be below the minimum required for ${subsample_count} subsamples at ${min_read_depth:-25}× depth (${_min_required} estimated). Assembly may be incomplete."
    fi
fi

# --- 9b. Subsample reads ---
log_info "9b. Subsampling reads (method=${subsampler}, count=${subsample_count})..."
if [[ "${subsampler}" == "autocycler" ]]; then
    subsample_flags=(
        --reads "${input_reads}"
        --out_dir subsampled_reads
        --genome_size "${genome_size}"
        --count "${subsample_count}"
        --seed "${subsample_seed}"
    )
    [[ -n "${min_read_depth}" ]] && subsample_flags+=(--min_read_depth "${min_read_depth}")
    run_cmd autocycler subsample "${subsample_flags[@]}"
else
    if [[ ! "${_genome_size_bp}" =~ ^[1-9][0-9]*$ || ! "${_avail_bases}" =~ ^[0-9]+$ || "${_avail_bases}" -eq 0 ]]; then
        log_error "Rasusa subsampling requires a numeric genome size and readable input base count."
    fi
    _total_depth=$(awk -v bases="${_avail_bases}" -v genome="${_genome_size_bp}" \
        'BEGIN { printf "%.3f", bases / genome }')
    if [[ -n "${rasusa_coverage}" ]]; then
        _rasusa_target="${rasusa_coverage}"
    else
        _rasusa_target=$(awk -v depth="${_total_depth}" -v minimum="${min_read_depth:-25}" \
            'BEGIN {
                ratio = 4 * depth / minimum
                if (ratio <= 1) exit 1
                printf "%.1f", minimum * (log(ratio) / log(2)) / 2
            }') || log_error "Input depth (${_total_depth}x) is too low to auto-calculate a positive Rasusa subset target. Set --rasusa-coverage explicitly or use native Autocycler subsampling."
    fi

    log_info "Subsampler: Rasusa"
    log_info "Input depth: ${_total_depth}x"
    log_info "Rasusa subset target: ${_rasusa_target}x"
    log_info "Subsets: ${subsample_count}"
    log_info "Base seed: ${subsample_seed}"
    if _float_gt "${_rasusa_target}" "${_total_depth}"; then
        log_warn "Rasusa target coverage (${_rasusa_target}x) exceeds input depth (${_total_depth}x); --strict will reject the run. Lower --rasusa-coverage."
    fi
    _sampling_fraction=$(awk -v target="${_rasusa_target}" -v total="${_total_depth}" \
        'BEGIN { if (total > 0) printf "%.3f", target / total; else print 0 }')
    if _float_gt "${_sampling_fraction}" 0.5; then
        log_warn "Rasusa subsets will be highly overlapping at the requested coverage. Native 'autocycler subsample' is preferable when maximal subset independence is important."
    fi

    mkdir -p subsampled_reads
    for ((i = 1; i <= subsample_count; i++)); do
        printf -v _sample_index '%02d' "${i}"
        _sample_seed=$(( subsample_seed + i - 1 ))
        run_cmd rasusa reads \
            --coverage "${_rasusa_target}" \
            --genome-size "${genome_size}" \
            --seed "${_sample_seed}" \
            --strict \
            --output "subsampled_reads/sample_${_sample_index}.fastq" \
            "${input_reads}"
    done

    if [[ -z "${dry_run:-}" ]]; then
        _rasusa_metadata="subsampled_reads/rasusa_subsample.tsv"
        printf 'sample\tmethod\tseed\ttarget_coverage\tread_count\tbases\tobserved_coverage\tn50\n' > "${_rasusa_metadata}"
        for ((i = 1; i <= subsample_count; i++)); do
            printf -v _sample_index '%02d' "${i}"
            _sample_seed=$(( subsample_seed + i - 1 ))
            _sample_file="subsampled_reads/sample_${_sample_index}.fastq"
            read -r _sample_count _sample_bases _sample_n50 <<< "$(seqkit stats -a -T "${_sample_file}" 2>/dev/null | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="num_seqs")c=i;if($i=="sum_len")b=i;if($i=="N50")n=i}} NR==2{print $c, $b, $n}')"
            _observed_coverage=$(awk -v bases="${_sample_bases:-0}" -v genome="${_genome_size_bp}" 'BEGIN { printf "%.2f", bases / genome }')
            printf 'sample_%s\trasusa\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${_sample_index}" "${_sample_seed}" "${_rasusa_target}" \
                "${_sample_count:-0}" "${_sample_bases:-0}" "${_observed_coverage}" "${_sample_n50:-0}" \
                >> "${_rasusa_metadata}"
        done
    else
        log_info "[DRY-RUN] Would write subsampled_reads/rasusa_subsample.tsv from seqkit stats."
    fi
fi

# --- 9c. Build job lists ---
if [[ -f "${assemblies_dir}/joblog.tsv" || -f "${assemblies_dir}/joblog_canu.tsv" ]]; then
    log_info "9c. Existing joblogs found. Reusing job lists for --resume..."
else
    log_info "9c. Building assembly job lists..."
    mkdir -p "${assemblies_dir}"
    rm -f "${assemblies_dir}/jobs.txt" "${assemblies_dir}/jobs_canu.txt"

    # --min_depth_rel: accept assemblies from subsamples with >= this fraction
    # of target depth. Default 0.1 (10%). Lower values tolerate sparser
    # subsamples but may produce fragmented assemblies.
    _canu_selected=false
    for assembler in "${assemblers[@]}"; do
        if [[ "${assembler}" == "canu" ]]; then
            _canu_selected=true
            continue
        fi
        extra="${assembler_args[$assembler]:-}"
        for i in $(seq -f "%02g" 1 "${subsample_count}"); do
            cmd="autocycler helper ${assembler} --reads subsampled_reads/sample_${i}.fastq --out_prefix ${assemblies_dir}/${assembler}_${i} --threads ${threads_per_job} --genome_size ${genome_size} --read_type ${read_type} --min_depth_rel ${min_depth_rel}"
            [[ -n "${extra}" ]] && cmd+=" -- ${extra}"
            echo "${cmd}" >> "${assemblies_dir}/jobs.txt"
        done
    done

    if [[ "${_canu_selected}" == "true" ]]; then
        for i in $(seq -f "%02g" 1 "${subsample_count}"); do
            cmd="autocycler helper canu --reads subsampled_reads/sample_${i}.fastq --out_prefix ${assemblies_dir}/canu_${i} --threads ${canu_threads_per_job} --genome_size ${genome_size} --read_type ${read_type} --min_depth_rel ${min_depth_rel}"
            [[ -n "${canu_extra_args}" ]] && cmd+=" -- ${canu_extra_args}"
            echo "${cmd}" >> "${assemblies_dir}/jobs_canu.txt"
        done
    fi

    log_info "Jobs: $(wc -l < "${assemblies_dir}/jobs.txt" 2>/dev/null || echo 0) general, $(wc -l < "${assemblies_dir}/jobs_canu.txt" 2>/dev/null || echo 0) Canu"
fi # End joblog guard

# --- 9d. Run assemblies ---
log_step "9d. Running assemblies..."

if [[ -f "${assemblies_dir}/jobs_canu.txt" && -s "${assemblies_dir}/jobs_canu.txt" ]]; then
    log_info "Running Canu (${canu_parallel_jobs} jobs × ${canu_threads_per_job} threads)..."
    set +e
    run_cmd parallel --nice 19 --bar --jobs "${canu_parallel_jobs}" \
        --joblog "${assemblies_dir}/joblog_canu.tsv" \
        --resume --resume-failed \
        --results "${assemblies_dir}/logs" --timeout 16h \
        < "${assemblies_dir}/jobs_canu.txt"
    parallel_exit_canu=$?
    set -e
else
    parallel_exit_canu=0
fi

if [[ ${parallel_exit_canu} -ne 0 ]]; then
    if [[ ${parallel_exit_canu} -eq 255 ]]; then
        log_error "GNU parallel failed catastrophically (exit 255). Aborting."
    fi
    actual_fails=$(awk -F'\t' 'NR>1 && $7!=0 {count++} END {print count+0}' \
        "${assemblies_dir}/joblog_canu.tsv" 2>/dev/null || echo 0)
    _sigkill_count=$(awk -F'\t' 'NR>1 && $8==9 {count++} END {print count+0}' \
        "${assemblies_dir}/joblog_canu.tsv" 2>/dev/null || echo 0)
    if [[ ${_sigkill_count} -gt 0 ]]; then
        log_warn "${_sigkill_count} Canu job(s) killed by SIGKILL (OOM likely). Consider reducing --canu-parallel-jobs."
    fi
    if [[ ${parallel_exit_canu} -ge 128 && ${actual_fails} -eq 0 ]]; then
        log_error "GNU parallel killed by signal $(( parallel_exit_canu - 128 )) (likely OOM). Aborting."
    elif [[ ${actual_fails} -gt 0 ]]; then
        log_warn "${actual_fails} Canu job(s) failed. Surviving assemblies will continue downstream."
        awk -F'\t' 'NR>1 && $7!=0 {print "  -> Failed: " $9}' "${assemblies_dir}/joblog_canu.tsv" >&2
    else
        log_error "Unknown GNU parallel failure (exit ${parallel_exit_canu}, 0 logged failures). Aborting."
    fi
fi

if [[ -f "${assemblies_dir}/jobs.txt" && -s "${assemblies_dir}/jobs.txt" ]]; then
    log_info "Running general assemblers (${parallel_jobs} jobs × ${threads_per_job} threads)..."
    set +e
    run_cmd parallel --nice 19 --bar --jobs "${parallel_jobs}" \
        --joblog "${assemblies_dir}/joblog.tsv" \
        --resume --resume-failed \
        --results "${assemblies_dir}/logs" --timeout 8h \
        < "${assemblies_dir}/jobs.txt"
    parallel_exit_general=$?
    set -e
else
    parallel_exit_general=0
fi

if [[ ${parallel_exit_general} -ne 0 ]]; then
    if [[ ${parallel_exit_general} -eq 255 ]]; then
        log_error "GNU parallel failed catastrophically (exit 255). Aborting."
    fi
    actual_fails=$(awk -F'\t' 'NR>1 && $7!=0 {count++} END {print count+0}' \
        "${assemblies_dir}/joblog.tsv" 2>/dev/null || echo 0)
    _sigkill_count=$(awk -F'\t' 'NR>1 && $8==9 {count++} END {print count+0}' \
        "${assemblies_dir}/joblog.tsv" 2>/dev/null || echo 0)
    if [[ ${_sigkill_count} -gt 0 ]]; then
        log_warn "${_sigkill_count} general assembler job(s) killed by SIGKILL (OOM likely). Consider reducing --parallel-jobs."
    fi
    if [[ ${parallel_exit_general} -ge 128 && ${actual_fails} -eq 0 ]]; then
        log_error "GNU parallel killed by signal $(( parallel_exit_general - 128 )) (likely OOM). Aborting."
    elif [[ ${actual_fails} -gt 0 ]]; then
        log_warn "${actual_fails} general assembler job(s) failed. Surviving assemblies will continue downstream."
        awk -F'\t' 'NR>1 && $7!=0 {print "  -> Failed: " $9}' "${assemblies_dir}/joblog.tsv" >&2
    else
        log_error "Unknown GNU parallel failure (exit ${parallel_exit_general}, 0 logged failures). Aborting."
    fi
fi

# --- 9e. Apply weighting ---
log_info "9e. Applying Autocycler weighting..."
shopt -s nullglob
if [[ -z "${dry_run:-}" ]]; then
    tmp_weight=$(mktemp)
    for f in "${assemblies_dir}"/plassembler*.fasta; do
        sed '/Autocycler_cluster_weight=/! s/circular=[Tt][Rr][Uu][Ee]/circular=True Autocycler_cluster_weight=3/I' "$f" > "$tmp_weight" && mv "$tmp_weight" "$f"
    done
    for f in "${assemblies_dir}"/canu*.fasta "${assemblies_dir}"/flye*.fasta "${assemblies_dir}"/hifiasm*.fasta; do
        sed '/Autocycler_consensus_weight=/! s/^>.*$/& Autocycler_consensus_weight=2/' "$f" > "$tmp_weight" && mv "$tmp_weight" "$f"
    done
    rm -f "$tmp_weight"
else
    log_info "[DRY-RUN] Would apply Autocycler weighting tags to plassembler/canu/flye/hifiasm FASTA headers."
fi
shopt -u nullglob

# --- 9f. Cleanup ---
run_cmd rm -f subsampled_reads/*.fastq

# ==============================================================================
# 9g. CURATION POINT 1 — Inspect assemblies
# ==============================================================================

failed_canu=0
[[ -f "${assemblies_dir}/joblog_canu.tsv" ]] && failed_canu=$(awk -F'\t' 'NR>1 && $7!=0' "${assemblies_dir}/joblog_canu.tsv" 2>/dev/null | wc -l || echo 0)
failed_general=0
[[ -f "${assemblies_dir}/joblog.tsv" ]] && failed_general=$(awk -F'\t' 'NR>1 && $7!=0' "${assemblies_dir}/joblog.tsv" 2>/dev/null | wc -l || echo 0)
empty_count=0
[[ -d "${assemblies_dir}" ]] && empty_count=$(find "${assemblies_dir}" -name "*.fasta" -size 0 2>/dev/null | wc -l || echo 0)

advice="DIAGNOSTICS:
  Failed jobs: Canu=${failed_canu}, General=${failed_general}
  Empty FASTA files: ${empty_count}

  Assembly sizes:"

for f in "${assemblies_dir}"/*.fasta; do
    [[ -f "$f" ]] && advice+="
    $(basename "$f"): $(grep -v '^>' "$f" | tr -d '\r\n' | wc -c) bp"
done

advice+="

ACTION: Delete empty or broken FASTA files, then continue."

manual_curation_pause "Inspect assemblies before clustering" "${advice}"
touch .success_assembly

if [[ "${_resume_target}" != "trim" && "${_resume_target}" != "dnaapler" ]]; then

# --- 9h. Compress & Cluster ---
log_step "9h. Compress & cluster"
run_cmd autocycler compress -i "${assemblies_dir}" -a "${autocycler_dir}" --threads "${compress_threads}"

_cluster_log="${autocycler_dir}/cluster_capture.log"
if [[ -z "${dry_run:-}" ]]; then
    autocycler cluster -a "${autocycler_dir}" 2>&1 | tee "${_cluster_log}"
else
    log_info "[DRY-RUN] autocycler cluster -a ${autocycler_dir}"
    mkdir -p "${autocycler_dir}"
    touch "${_cluster_log}"
fi

# --- Cluster QC summary ---
log_info "--- Cluster QC summary ---"
awk '
/^Cluster [0-9]+:/ { cid=$0; next }
/cluster distance:/  { dist=$0 }
/passed QC$/         { printf "%s  %s  → passed QC\n", cid, dist; dist="" }
/failed QC:/         { printf "%s  %s  → %s\n", cid, dist, $0; dist="" }
' "${_cluster_log}" | while IFS= read -r line; do log_info "  $line"; done

_pass_count=$(grep -c "passed QC$" "${_cluster_log}" || true)
_fail_count=$(grep -c "failed QC:" "${_cluster_log}" || true)
log_info "Clusters passed: ${_pass_count}, failed QC: ${_fail_count}"
rm -f "${_cluster_log}"

if [[ -f "${autocycler_dir}/clustering/clustering.newick" ]]; then
    log_info "Running ETE3 sanity checks on clustering tree..."
    python3 - "${autocycler_dir}/clustering/clustering.newick" << 'PYEOF' || true
import sys, re
try:
    from ete3 import Tree
    tree_file = sys.argv[1]
    t = Tree(tree_file)

    # Compute median branch length for a data-driven outlier threshold
    dists = sorted([n.dist for n in t.traverse() if not n.is_root() and n.dist > 0])
    if dists:
        median_dist = dists[len(dists) // 2]
        outlier_threshold = median_dist * 5  # 5× median as outlier cutoff
    else:
        outlier_threshold = float('inf')

    for node in t.traverse("postorder"):
        if not node.is_leaf() and node.dist > outlier_threshold:
            print(f"WARNING: High branch length ({node.dist:.4f}) exceeds "
                  f"5× median ({median_dist:.4f}) — potential outlier cluster.",
                  file=sys.stderr)
        if not node.is_leaf() and len(node.get_leaves()) > 3:
            # Extract assembler name: everything before the last _NN suffix
            assemblers = [re.sub(r'_\d+.*$', '', l.name.strip("'\"")) for l in node.get_leaves()]
            if len(set(assemblers)) == 1:
                print(f"WARNING: Clade with {len(assemblers)} leaves all from "
                      f"assembler '{assemblers[0]}' — possible miscluster.",
                      file=sys.stderr)
except ImportError:
    print("ETE3 not installed — skipping cluster sanity check.", file=sys.stderr)
except Exception as e:
    print(f"ETE3 parsing error: {e}", file=sys.stderr)
PYEOF
fi

# ==============================================================================
# 9i. CURATION POINT 2 — Inspect clustering
# ==============================================================================

cluster_advice="CLUSTER SUMMARY:
$(column -t -s$'\t' "${autocycler_dir}/clustering/clustering.tsv" 2>/dev/null || cat "${autocycler_dir}/clustering/clustering.tsv" 2>/dev/null || echo "  (clustering.tsv not found)")

GUIDE:
  Large cluster (~genome size) = chromosome → qc_pass/
  Small clusters (1-200 kb)    = plasmids   → qc_pass/
  Tiny clusters (<500 bp)      = noise      → qc_fail/
  
  Bandage integration: If enabled, load GFAs in ${autocycler_dir}/clustering/qc_pass/*/1_untrimmed.gfa to inspect nodes.

ACTION: Move misclassified clusters between qc_pass/ and qc_fail/, then continue."

manual_curation_pause "Inspect clustering results" "${cluster_advice}"
touch .success_cluster
fi # End resume skip for cluster stage

# --- Hybrid dotplot function ---
# Small clusters (< 2 MB GFA): Autocycler dotplot (fast, all-vs-all pairwise)
# Large clusters (>= 2 MB GFA): MUMmer nucmer + mummerplot (chromosome-scale)
generate_dotplots() {
    local cluster_dir="$1"
    # NOTE: cname is inherited from the caller's loop scope (not passed as arg)
    local gfa_untrimmed="${cluster_dir}/1_untrimmed.gfa"
    local gfa_trimmed="${cluster_dir}/2_trimmed.gfa"
    local gfa_seq_len

    if [[ -n "${dry_run:-}" ]]; then
        log_info "  [DRY-RUN] Skipping generate_dotplots for ${cluster_dir}"
        return 0
    fi

    gfa_seq_len=$(awk '/^S/{sum+=length($3)} END{print sum+0}' "${gfa_untrimmed}" 2>/dev/null || echo 0)

    # 2 Mbp boundary between plasmid-scale and chromosome-scale clusters
    if [[ ${gfa_seq_len} -lt 2000000 ]]; then
        log_info "  Plasmid-scale cluster (${gfa_seq_len} bp) — using Autocycler dotplot"
        run_cmd autocycler dotplot -i "${gfa_untrimmed}" -o "${cluster_dir}/1_untrimmed.png"
        if [[ -s "${gfa_trimmed}" ]]; then
            run_cmd autocycler dotplot -i "${gfa_trimmed}" -o "${cluster_dir}/2_trimmed.png"
        else
            log_warn "  ${cname}: 2_trimmed.gfa missing or empty — trim may have excluded all sequences. Skipping trimmed dotplot."
        fi
    elif command -v nucmer &>/dev/null && command -v mummerplot &>/dev/null; then
        log_info "  Chromosome-scale cluster (${gfa_seq_len} bp) — using MUMmer (nucmer + mummerplot)"
        for tag in 1_untrimmed 2_trimmed; do
            local gfa="${cluster_dir}/${tag}.gfa"
            if [[ ! -s "${gfa}" ]]; then
                log_warn "  ${cname}: ${tag}.gfa missing or empty — skipping MUMmer dotplot for this tag."
                continue
            fi
            local fasta="${cluster_dir}/${tag}.fasta"
            # Convert GFA to FASTA for nucmer
            run_cmd autocycler gfa2fasta -i "${gfa}" -o "${fasta}"
            # Self-vs-self alignment to reveal overlaps and structural issues
            run_cmd nucmer --maxmatch --nosimplify \
                "${fasta}" "${fasta}" \
                -p "${cluster_dir}/${tag}_selfmap"
            run_cmd mummerplot --png --large \
                -p "${cluster_dir}/${tag}" \
                "${cluster_dir}/${tag}_selfmap.delta"
            # Clean up intermediate files
            rm -f "${cluster_dir}/${tag}_selfmap.delta" \
                  "${cluster_dir}/${tag}_selfmap.cluster" \
                  "${cluster_dir}/${tag}.fasta" \
                  "${cluster_dir}/${tag}.filt" \
                  "${cluster_dir}/${tag}.gp"
        done
    else
        log_warn "  Large cluster but nucmer/mummerplot not found — skipping dotplots"
    fi
}

# PROKARYONT_RESUME=dnaapler: skip trim & resolve (9j) and assume autocycler_out/
# already contains valid 5_final.gfa files. This is intended for cases where 03
# completed through combine (9l) but 04 (dnaapler) needs a clean re-run.
# --- 9j. Trim & Resolve ---
if [[ "${_resume_target}" != "dnaapler" ]]; then
log_step "9j. Trim & resolve"
trim_flags=()
[[ -n "${trim_mad}" ]]          && trim_flags+=(--mad "${trim_mad}")
[[ -n "${trim_min_identity}" ]] && trim_flags+=(--min_identity "${trim_min_identity}")

# Initialise per-cluster metadata TSV (joined with depth data later)
_cluster_meta="${autocycler_dir}/cluster_metadata.tsv"
echo -e "cluster\tcluster_distance\tmedian_len\tmad\tallowed_range\tse_trimmed\thp_trimmed\tkept\texcluded\tanchors\tunique_bridges\tconflicting_bridges\tconsensus_unitigs\tconsensus_length" \
    > "${_cluster_meta}"

shopt -s nullglob
_clusters=("${autocycler_dir}"/clustering/qc_pass/cluster_*)
shopt -u nullglob

if [[ ${#_clusters[@]} -eq 0 && -n "${dry_run:-}" ]]; then
    _clusters=("${autocycler_dir}/clustering/qc_pass/cluster_1")
fi

for c in "${_clusters[@]}"; do
    cname=$(basename "$c")
    log_info "Processing ${cname}..."

    # --- Trim with capture ---
    _trim_log="${c}/trim_capture.log"
    if [[ -z "${dry_run:-}" ]]; then
        autocycler trim -c "$c" ${trim_flags[@]+"${trim_flags[@]}"} 2>&1 | tee "${_trim_log}"
    else
        log_info "[DRY-RUN] autocycler trim -c $c ${trim_flags[*]+${trim_flags[*]}}"
        mkdir -p "$(dirname "${_trim_log}")"
        touch "${_trim_log}"
    fi

    # Parse trim stats
    _se_count=$(awk '/^Trim start-end overlaps/,/^Trim hairpin overlaps/{if(/trimmed to/)c++} END{print c+0}' "${_trim_log}")
    _hp_count=$(awk '/^Trim hairpin overlaps/,/^Exclude outliers/{if(/trimmed to/)c++} END{print c+0}' "${_trim_log}")
    _median=$(grep "^Median sequence length:" "${_trim_log}" | awk '{print $NF}' || echo "?")
    _mad_val=$(grep "^Median absolute deviation:" "${_trim_log}" | awk '{print $NF}' || echo "?")
    _range=$(grep "^Allowed length range:" "${_trim_log}" | awk '{print $NF}' || echo "?")
    _excluded=$(grep -c ": excluded$" "${_trim_log}" || echo 0)
    _kept=$(grep -c ": kept$" "${_trim_log}" || echo 0)

    log_info "  Trim: ${_se_count} start-end, ${_hp_count} hairpin trimmed"
    log_info "  MAD filter: median=${_median}, MAD=${_mad_val}, range=${_range}"
    log_info "  Contigs: ${_kept} kept, ${_excluded} excluded"
    rm -f "${_trim_log}"

    # --- Dotplots ---
    generate_dotplots "$c"

    # --- Resolve with capture ---
    _resolve_log="${c}/resolve_capture.log"
    if [[ -z "${dry_run:-}" ]]; then
        autocycler resolve -c "$c" 2>&1 | tee "${_resolve_log}"
    else
        log_info "[DRY-RUN] autocycler resolve -c $c"
        touch "${_resolve_log}"
    fi

    _anchors=$(grep 'anchor unitigs found' "${_resolve_log}" \
        | awk '{
            for(i=1; i<=NF; i++) {
                if($i == "anchor" && i > 1) {
                    val = $(i-1);
                    gsub(/,/, "", val);
                    if(val ~ /^[0-9]+$/) { print val; exit }
                }
            }
        }' || echo "?")

    _unique=$(grep 'Unique bridges:' "${_resolve_log}" \
        | awk '{
            for(i=1; i<=NF; i++) {
                if($i == "Unique" && $(i+1) == "bridges:" && (i+1) < NF) {
                    val = $(i+2);
                    gsub(/,/, "", val);
                    if(val ~ /^[0-9]+$/) { print val; exit }
                }
            }
        }' || echo "?")

    _conflict=$(grep 'Conflicting bridges:' "${_resolve_log}" \
        | awk '{
            for(i=1; i<=NF; i++) {
                if($i == "Conflicting" && $(i+1) == "bridges:" && (i+1) < NF) {
                    val = $(i+2);
                    gsub(/,/, "", val);
                    if(val ~ /^[0-9]+$/) { print val; exit }
                }
            }
        }' || echo "?")

    _final_unitigs=$(grep 'unitig' "${_resolve_log}" | tail -1 \
        | awk '{
            for(i=1; i<=NF; i++) {
                if($i ~ /^unitigs?$/ && i > 1) {
                    val = $(i-1);
                    gsub(/,/, "", val);
                    if(val ~ /^[0-9]+$/) { print val; exit }
                }
            }
        }' || echo "?")

    _total_len=$(grep 'total length:' "${_resolve_log}" | tail -1 \
        | awk '{
            for(i=1; i<=NF; i++) {
                if($i == "length:" && i > 1 && $(i-1) == "total" && i < NF) {
                    val = $(i+1);
                    gsub(/,/, "", val);
                    if(val ~ /^[0-9]+$/) { print val; exit }
                }
            }
        }' || echo "?")

    log_info "  Resolve: ${_anchors} anchors, ${_unique} unique / ${_conflict} conflicting bridges"
    log_info "  Result: ${_final_unitigs} unitig(s), total length: ${_total_len}"
    rm -f "${_resolve_log}"

    # --- Recover cluster distance from clustering.tsv ---
    _cluster_dist=$(awk -F'\t' -v cn="${cname}" '$1==cn{print $2}' \
        "${autocycler_dir}/clustering/clustering.tsv" 2>/dev/null || echo "?")

    # --- Append to cluster metadata ---
    echo -e "${cname}\t${_cluster_dist}\t${_median}\t${_mad_val}\t${_range}\t${_se_count}\t${_hp_count}\t${_kept}\t${_excluded}\t${_anchors}\t${_unique}\t${_conflict}\t${_final_unitigs}\t${_total_len}" \
        >> "${_cluster_meta}"
done
touch .success_trim
fi # End resume skip for trim stage

# ==============================================================================
# 9k. CURATION POINT 3 — Inspect dotplots
# ==============================================================================

dotplot_advice="DOTPLOT FILES:"
for c in "${autocycler_dir}"/clustering/qc_pass/cluster_*; do
    [[ -f "$c/1_untrimmed.png" ]] && dotplot_advice+="
  $(basename "$c"): $c/1_untrimmed.png, $c/2_trimmed.png"
done
dotplot_advice+="

LOOK FOR: clean diagonal lines; overlap triangles removed after trimming.
RED FLAGS: jagged diagonal, off-diagonal blocks, triangle still present.

ACTION: Move bad clusters to qc_fail/, then continue."

manual_curation_pause "Inspect dotplots" "${dotplot_advice}"

# --- 9l. Combine ---
log_step "9l. Combining consensus assembly"

shopt -s nullglob
_gfas=("${autocycler_dir}"/clustering/qc_pass/cluster_*/5_final.gfa)
shopt -u nullglob

_combine_log="${autocycler_dir}/combine_capture.log"
if [[ ${#_gfas[@]} -eq 0 ]]; then
    if [[ -z "${dry_run:-}" ]]; then
        log_error "No 5_final.gfa files found in qc_pass/. All clusters may have failed QC."
    else
        log_warn "[DRY-RUN] No 5_final.gfa files found in qc_pass/. Using dry-run placeholders."
        _gfas=("${autocycler_dir}/clustering/qc_pass/cluster_1/5_final.gfa")
    fi
fi

combine_cmd=(
    autocycler combine
    -a "${autocycler_dir}"
    -i "${_gfas[@]}"
    --reads "${combine_reads}"
    --depth_kmer "${autocycler_combine_depth_kmer}"
    --threads "${combine_threads}"
)
if [[ -z "${dry_run:-}" ]]; then
    "${combine_cmd[@]}" 2>&1 | tee "${_combine_log}"
else
    run_cmd "${combine_cmd[@]}"
    touch "${_combine_log}"
fi

# Log combine summary: topology and resolved status
log_info "--- Combine summary ---"
{ grep -E "(circular|linear|total length|fully resolved)" "${_combine_log}" || true; } \
    | sed 's/^/  /' | while IFS= read -r line; do log_info "$line"; done
rm -f "${_combine_log}"

# Consensus file stays inside autocycler_out/ — downstream scripts reference it there
autocycler_consensus="${autocycler_dir}/consensus_assembly.fasta"

if [[ -z "${dry_run:-}" && ! -s "${autocycler_consensus}" ]]; then
    log_error "Autocycler combine did not produce a non-empty ${autocycler_consensus}."
fi

# --- 9m. Characterize every consensus contig against PLSDB ---
if [[ "${skip_plsdb_screen}" == "true" ]]; then
    log_info "9m. PLSDB screen skipped by configuration."
    if [[ -f "$(pwd)/plassembler_summary.tsv" ]]; then
        log_warn "An existing plassembler_summary.tsv was not refreshed because the PLSDB screen was skipped."
    fi
else
    log_step "9m. Screening consensus contigs against PLSDB with Plassembler assembled"
    plassembler_plsdb_dir="${autocycler_dir}/plassembler_plsdb"
    plassembler_summary_source="${plassembler_plsdb_dir}/plassembler_summary.tsv"
    plassembler_summary="$(pwd)/plassembler_summary.tsv"
    plassembler_cmd=(
        plassembler assembled
        --no_copy_numbers
        --database "${plassembler_db}"
        --input_plasmids "${autocycler_consensus}"
        --outdir "${plassembler_plsdb_dir}"
        --prefix plassembler
        --threads "${threads}"
        --force
    )

    if [[ -n "${dry_run:-}" ]]; then
        run_cmd "${plassembler_cmd[@]}"
        log_info "[DRY-RUN] Would validate one PLSDB summary row per consensus contig and publish ${plassembler_summary}."
    else
        _consensus_checksum_before=$(cksum "${autocycler_consensus}" | awk '{print $1 ":" $2}')
        "${plassembler_cmd[@]}"
        [[ -s "${plassembler_summary_source}" ]] || log_error "Plassembler completed without a non-empty ${plassembler_summary_source}."

        _consensus_contigs=$(grep -c '^>' "${autocycler_consensus}" || true)
        if [[ "${_consensus_contigs}" -lt 1 ]]; then
            log_error "The Autocycler consensus contains no FASTA records."
        fi
        _summary_rows=$(awk 'NR > 1 && $0 !~ /^[[:space:]]*$/ { rows++ } END { print rows + 0 }' "${plassembler_summary_source}")
        if [[ "${_summary_rows}" -ne "${_consensus_contigs}" ]]; then
            log_error "Incomplete PLSDB summary: ${_summary_rows} data rows for ${_consensus_contigs} consensus contigs."
        fi

        _consensus_checksum_after=$(cksum "${autocycler_consensus}" | awk '{print $1 ":" $2}')
        if [[ "${_consensus_checksum_before}" != "${_consensus_checksum_after}" ]]; then
            log_error "The authoritative Autocycler consensus changed during the PLSDB screen; refusing to publish the summary."
        fi
        cp "${plassembler_summary_source}" "${plassembler_summary}"
        log_info "PLSDB summary: ${plassembler_summary} (${_summary_rows} contigs)"
    fi
fi

# --- 9n. Metrics ---
if [[ -f "subsampled_reads/subsample.yaml" ]]; then
    run_cmd cp "subsampled_reads/subsample.yaml" .
fi
if [[ -f "subsampled_reads/rasusa_subsample.tsv" ]]; then
    run_cmd cp "subsampled_reads/rasusa_subsample.tsv" .
fi

# Extract initial read stats (needs -a for N50 column)
log_info "Collecting input read metrics via seqkit stats..."
read -r in_count in_bases in_n50 <<< $(seqkit stats -a -T "${input_reads}" 2>/dev/null \
    | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="num_seqs")c=i;if($i=="sum_len")b=i;if($i=="N50")n=i}} NR==2{print $c, $b, $n}' \
    || echo "0 0 0")
in_count=${in_count:-0}
in_bases=${in_bases:-0}
in_n50=${in_n50:-0}

run_cmd bash -c 'autocycler table > "$1"' _ metrics.tsv
run_cmd bash -c 'autocycler table -a "$1" -n "$2" >> "$3"' \
    _ "${autocycler_dir}" "${sample_name}" metrics.tsv

# Smart column insertion: avoid duplicating input_read columns
if [[ -z "${dry_run:-}" ]]; then
    _header=$(head -1 metrics.tsv)
    if echo "${_header}" | grep -q "input_read_count"; then
        # Columns exist — fill the empty values in the data row
        awk -F'\t' -v c="${in_count}" -v b="${in_bases}" -v n="${in_n50}" '
            BEGIN { OFS="\t" }
            NR==1 {
                for (i=1; i<=NF; i++) {
                    if ($i == "input_read_count") ci = i
                    if ($i == "input_read_bases") bi = i
                    if ($i == "input_read_n50")   ni = i
                }
                print; next
            }
            NR==2 {
                if (ci && $ci == "") $ci = c
                if (bi && $bi == "") $bi = b
                if (ni && $ni == "") $ni = n
                print; next
            }
            { print }
        ' metrics.tsv > metrics.tsv.tmp && mv metrics.tsv.tmp metrics.tsv
        log_info "Filled existing input_read columns in metrics.tsv"
    else
        # Columns don't exist: prepend them (older Autocycler versions)
        sed -i "1 s/^/input_read_count\tinput_read_bases\tinput_read_n50\t/" metrics.tsv
        sed -i "2 s/^/${in_count}\t${in_bases}\t${in_n50}\t/" metrics.tsv
        log_info "Prepended input_read columns to metrics.tsv"
    fi
else
    log_info "[DRY-RUN] Bypassing metrics.tsv post-processing."
fi

# ==============================================================================
# STEP 9o — Read-depth assessment
# ==============================================================================

log_step "9o. Read-depth assessment"
depth_report="$(pwd)/contig_depths.tsv"

if [[ -s "${depth_report}" ]]; then
    log_info "Found existing depth report. Skipping read-depth assessment..."
else
    # Map read type to minimap2 preset
    case "${read_type}" in
        ont_r9|ont_r10)  mm2_preset="map-ont"  ;;
        pacbio_clr)      mm2_preset="map-pb"   ;;
        pacbio_hifi)     mm2_preset="map-hifi" ;;
        *)               mm2_preset="map-ont"  ;;
    esac

    log_info "Mapping input reads to consensus (preset: ${mm2_preset})..."
    _mm2_threads=$(( threads * 3 / 4 ))
    _sort_threads=$(( threads / 4 ))
    [[ ${_mm2_threads} -lt 1 ]] && _mm2_threads=1
    [[ ${_sort_threads} -lt 1 ]] && _sort_threads=1

    run_cmd bash -c '
        set -o pipefail
        minimap2 -t "$1" -a -x "$2" "$3" "$4" \
            | samtools sort -@ "$5" -o consensus_mapped.bam
    ' _ "${_mm2_threads}" "${mm2_preset}" "${autocycler_consensus}" "${input_reads}" "${_sort_threads}"

    run_cmd samtools index consensus_mapped.bam

    # Compute per-contig depths and relative depth vs longest contig (chromosome)
    log_info "Computing per-contig depths..."
    run_cmd bash -c '
        depth_report="$1"
        echo -e "contig\tlength_bp\tavg_depth\trelative_depth" > "${depth_report}"

        # Collect per-contig average depth via samtools coverage (one line per contig)
        samtools coverage consensus_mapped.bam \
            | awk '\''{
                if (NR == 1) next          # skip header
                name  = $1
                len   = $3                 # endpos (= contig length for full coverage)
                depth = $7                 # meandepth
                d[name]   = depth
                l[name]   = len
            }
            END {
                max_len = 0; chr_name = ""; chr_depth = 1
                for (c in l) {
                    if (l[c] > max_len) {
                        max_len  = l[c]
                        chr_name = c
                        chr_depth = d[c]
                    }
                }
                for (c in l) {
                    printf "%s\t%d\t%.1f\t%.2f\n", c, l[c], d[c], d[c] / chr_depth
                }
            }'\'' \
            | sort -t$'\''\t'\'' -k2 -rn >> "${depth_report}"
    ' _ "${depth_report}"

    log_info "Depth report: ${depth_report}"

    run_cmd rm -f consensus_mapped.bam consensus_mapped.bam.bai
fi

# ==============================================================================
# STEP 9p — Contig depth summary
# ==============================================================================

log_step "9p. Building contig depth summary"
_cluster_meta="${autocycler_dir}/cluster_metadata.tsv"

if [[ -f "${depth_report}" ]]; then
    log_info "--- Contig depth summary ---"
    while IFS=$'\t' read -r contig length avg_depth rel_depth; do
        [[ "${contig}" == "contig" ]] && continue
        flag=""
        if _float_gt "${rel_depth}" 1.5; then
            flag=" ← elevated copy number"
        fi
        log_info "  ${contig} (${length} bp): ${avg_depth}× avg, ${rel_depth}× relative${flag}"
    done < "${depth_report}"
fi

log_step "Assembly complete."
log_info "Consensus: ${autocycler_consensus}"
if [[ "${skip_plsdb_screen}" != "true" ]]; then
    log_info "PLSDB:    $(pwd)/plassembler_summary.tsv"
fi
log_info "Metrics:   metrics.tsv"
log_info "Depths:    ${depth_report}"
log_info "Cluster metadata: ${_cluster_meta}"

# RagTag scaffolding has been removed from this script and deferred to a
# dedicated post-assembly scaffolding script (not yet implemented). Avoid a
# silent regression for anyone already using --ragtag:
if [[ -n "${PROKARYONT_RAGTAG:-}" ]]; then
    log_warn "RagTag scaffolding has moved out of 03_autocycler_assemble.sh pending a dedicated post-assembly script (not yet implemented). --ragtag is currently a no-op here."
fi
