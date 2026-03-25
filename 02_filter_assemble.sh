#!/usr/bin/env bash
# ==============================================================================
# 02_filter_assemble.sh — Read filtering & Autocycler assembly pipeline
# ==============================================================================
# Usage:
#   bash 02_filter_assemble.sh --input-fastq reads.fq.gz --read-type ont_r10
#
# Optional flags:
#   --config FILE           Path to pipeline.conf (values override defaults)
#   --threads N             Number of threads (default: 128)
#   --sample-name NAME      Sample name for metrics (default: MyBacteria)
#   --min-length N          Filtlong minimum read length (default: 200)
#   --min-qscore N          Seqkit minimum average read quality (default: 7)
#   --keep-percent N        Filtlong keep percent (default: 90)
#   --genome-size SIZE      Override genome size (skip re-estimation)
#   --subsample-count N     Number of read subsamples (default: 4)
#   --min-read-depth N      Min subset read depth for subsampling (default: 25)
#   --parallel-jobs N       Concurrent assembler jobs (default: 4)
#   --canu-parallel-jobs N  Concurrent Canu jobs (default: 2)
#   --skip-curation         Skip manual curation pauses
#   --dry-run               Print commands without executing
#   --help                  Show this help
# ==============================================================================

source "$(dirname "$0")/00_setup.sh"

# --- Defaults ----------------------------------------------------------------
threads="${threads:-128}"
input_fastq="${input_fastq:-}"
read_type="${read_type:-ont_r10}"
sample_name="${sample_name:-MyBacteria}"
filtlong_min_length="${filtlong_min_length:-200}"
min_qscore="${min_qscore:-7}"
filtlong_keep_percent="${filtlong_keep_percent:-90}"
genome_size_override="${genome_size_override:-}"
subsample_count="${subsample_count:-4}"
min_read_depth="${min_read_depth:-}"
parallel_jobs="${parallel_jobs:-4}"
canu_parallel_jobs="${canu_parallel_jobs:-2}"
skip_curation="${skip_curation:-}"
canu_extra_args="${canu_extra_args:-}"
config_file="${config_file:-}"

# Assembler extra args (associative array — set via config file or before sourcing)
declare -A assembler_args 2>/dev/null || true
: "${assembler_args[flye]:=}"
: "${assembler_args[metamdbg]:=}"
: "${assembler_args[miniasm]:=}"
: "${assembler_args[raven]:=}"
: "${assembler_args[plassembler]:=}"
: "${assembler_args[hifiasm]:=}"
: "${assembler_args[wtdbg2]:=}"
: "${assembler_args[myloasm]:=}"

# --- Usage -------------------------------------------------------------------
usage() {
    echo "Usage: $(basename "$0") --input-fastq FILE [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --input-fastq FILE          Raw FASTQ file (gzipped ok)"
    echo ""
    echo "Optional:"
    echo "  --config FILE               Path to pipeline.conf"
    echo "  --read-type TYPE            Read type: ont_r9|ont_r10|pacbio_clr|pacbio_hifi (default: ont_r10)"
    echo "  --threads N                 Number of threads (default: 128)"
    echo "  --sample-name NAME          Sample name for metrics (default: MyBacteria)"
    echo "  --min-length N              Filtlong min read length in bp (default: 200)"
    echo "  --min-qscore N              Seqkit min average read quality (default: 7)"
    echo "  --keep-percent N            Filtlong keep percent (default: 90)"
    echo "  --genome-size SIZE          Override genome size (skip re-estimation)"
    echo "  --subsample-count N         Number of read subsamples (default: 4)"
    echo "  --min-read-depth N          Min subset read depth for subsampling (default: 25)"
    echo "  --parallel-jobs N           Concurrent assembler jobs (default: 4)"
    echo "  --canu-parallel-jobs N      Concurrent Canu jobs (default: 2)"
    echo "  --skip-curation             Skip manual curation pauses"
    echo "  --dry-run                   Print commands without executing"
    echo "  --help                      Show this help"
    exit 0
}

# --- Parse arguments ---------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-fastq)        input_fastq="$2"; shift 2 ;;
        --read-type)          read_type="$2"; shift 2 ;;
        --config)             config_file="$2"; shift 2 ;;
        --threads)            threads="$2"; shift 2 ;;
        --sample-name)        sample_name="$2"; shift 2 ;;
        --min-length)         filtlong_min_length="$2"; shift 2 ;;
        --min-qscore)         min_qscore="$2"; shift 2 ;;
        --keep-percent)       filtlong_keep_percent="$2"; shift 2 ;;
        --genome-size)        genome_size_override="$2"; shift 2 ;;
        --subsample-count)    subsample_count="$2"; shift 2 ;;
        --min-read-depth)     min_read_depth="$2"; shift 2 ;;
        --parallel-jobs)      parallel_jobs="$2"; shift 2 ;;
        --canu-parallel-jobs) canu_parallel_jobs="$2"; shift 2 ;;
        --skip-curation)      skip_curation=true; shift ;;
        --dry-run)            dry_run=true; shift ;;
        --help|-h)            usage ;;
        *) log_error "Unknown flag: $1. Use --help for usage." ;;
    esac
done

# --- Load config (CLI flags parsed above take priority) ----------------------
[[ -n "${config_file}" ]] && load_config "${config_file}"

# --- Validate ----------------------------------------------------------------
require_arg "--input-fastq" "${input_fastq}"
require_file "${input_fastq}" ""

require_tool filtlong
require_tool seqkit
require_tool NanoPlot
require_tool autocycler
require_tool parallel

for tool in canu flye metaMDBG miniasm minipolish minimap2 raven \
            plassembler hifiasm wtdbg2 wtpoa-cns myloasm \
            nucmer mummerplot; do
    command -v "${tool}" &>/dev/null || log_warn "'${tool}' not found — its jobs will fail."
done

# --- Derived values ----------------------------------------------------------
threads_per_job=$(( threads / parallel_jobs ))
canu_threads_per_job=$(( threads / canu_parallel_jobs ))
qc_dir="$(pwd)/01_qc"
genome_size_dir="$(pwd)/02_genome_size"
assemblies_dir="$(pwd)/assemblies"
autocycler_dir="$(pwd)/autocycler_out"
filtered_reads="$(pwd)/filtered_input.fastq.gz"
autocycler_consensus="$(pwd)/autocycler_consensus.fasta"

# ==============================================================================
# STEP 2 — Filtering with Filtlong
# ==============================================================================

log_step "Step 2: Filtering reads (min_qscore=${min_qscore}, min_length=${filtlong_min_length}, keep=${filtlong_keep_percent}%)"

run_cmd bash -c 'seqkit seq -Q "$1" "$2" 2>/dev/null | filtlong --min_length "$3" --keep_percent "$4" - | ${GZIP_BIN} > "$5"' \
    _ "${min_qscore}" "${input_fastq}" "${filtlong_min_length}" "${filtlong_keep_percent}" "${filtered_reads}"

log_info "Running NanoPlot on filtered reads..."
mkdir -p "${qc_dir}"
run_cmd NanoPlot --threads "${threads}" -c green \
    --fastq "${filtered_reads}" \
    -o "${qc_dir}/NanoPlot_FiltPol" \
    --loglength --plots dot --N50

log_info ">>> CHECK: Compare ${qc_dir}/NanoPlot_FiltPol/ to ${qc_dir}/NanoPlot_sample/"

# ==============================================================================
# STEP 4 — Autocycler Assembly
# ==============================================================================

log_step "Step 4: Autocycler assembly pipeline"

# --- 4a. Genome size estimation ---
if [[ -n "${genome_size_override}" ]]; then
    genome_size="${genome_size_override}"
    log_info "4a. Using user-provided genome size: ${genome_size}"
else
    log_info "4a. Estimating genome size..."
    if [[ -n "${dry_run:-}" ]]; then
        log_info "[DRY-RUN] autocycler helper genome_size --reads ${filtered_reads} --threads ${threads}"
        genome_size="DRY_RUN_PLACEHOLDER"
    else
        genome_size=$(autocycler helper genome_size \
            --reads "${filtered_reads}" --threads "${threads}")
    fi
    log_info "Autocycler genome size: ${genome_size}"
fi

if [[ -f "${genome_size_dir}/lrge_output.txt" ]]; then
    log_info "LRGE estimate: $(cat "${genome_size_dir}/lrge_output.txt")"
fi

# --- 4b. Subsample reads ---
log_info "4b. Subsampling reads (count=${subsample_count})..."
subsample_flags=(
    --reads "${filtered_reads}"
    --out_dir subsampled_reads
    --genome_size "${genome_size}"
    --count "${subsample_count}"
)
[[ -n "${min_read_depth}" ]] && subsample_flags+=(--min_read_depth "${min_read_depth}")
run_cmd autocycler subsample "${subsample_flags[@]}"

# --- 4c. Build job lists ---
log_info "4c. Building assembly job lists..."
mkdir -p "${assemblies_dir}"
rm -f "${assemblies_dir}/jobs.txt" "${assemblies_dir}/jobs_canu.txt"

for assembler in $(echo "${!assembler_args[@]}" | tr ' ' '\n' | sort); do
    extra="${assembler_args[$assembler]}"
    for i in $(seq -f "%02g" 1 "${subsample_count}"); do
        cmd="autocycler helper ${assembler} --reads subsampled_reads/sample_${i}.fastq --out_prefix ${assemblies_dir}/${assembler}_${i} --threads ${threads_per_job} --genome_size ${genome_size} --read_type ${read_type} --min_depth_rel 0.1"
        [[ -n "${extra}" ]] && cmd+=" -- ${extra}"
        echo "${cmd}" >> "${assemblies_dir}/jobs.txt"
    done
done

for i in $(seq -f "%02g" 1 "${subsample_count}"); do
    cmd="autocycler helper canu --reads subsampled_reads/sample_${i}.fastq --out_prefix ${assemblies_dir}/canu_${i} --threads ${canu_threads_per_job} --genome_size ${genome_size} --read_type ${read_type} --min_depth_rel 0.1"
    [[ -n "${canu_extra_args}" ]] && cmd+=" -- ${canu_extra_args}"
    echo "${cmd}" >> "${assemblies_dir}/jobs_canu.txt"
done

log_info "Jobs: $(wc -l < "${assemblies_dir}/jobs.txt") general, $(wc -l < "${assemblies_dir}/jobs_canu.txt") Canu"

# --- 4d. Run assemblies ---
log_step "4d. Running assemblies..."

log_info "Running Canu (${canu_parallel_jobs} jobs × ${canu_threads_per_job} threads)..."
set +e
run_cmd nice -n 19 parallel --jobs "${canu_parallel_jobs}" \
    --joblog "${assemblies_dir}/joblog_canu.tsv" \
    --results "${assemblies_dir}/logs" --timeout 16h \
    < "${assemblies_dir}/jobs_canu.txt"
parallel_exit_canu=$?
set -e

if [[ ${parallel_exit_canu} -ne 0 ]]; then
    if [[ ${parallel_exit_canu} -ge 1 && ${parallel_exit_canu} -le 254 ]]; then
        log_warn "${parallel_exit_canu} Canu job(s) failed. Surviving assemblies will continue downstream."
        awk -F'\t' 'NR>1 && $7!=0 {print "  -> Failed: " $9}' "${assemblies_dir}/joblog_canu.tsv" >&2
    else
        log_error "GNU parallel failed catastrophically (exit code ${parallel_exit_canu}). Aborting."
    fi
fi

log_info "Running general assemblers (${parallel_jobs} jobs × ${threads_per_job} threads)..."
set +e
run_cmd nice -n 19 parallel --jobs "${parallel_jobs}" \
    --joblog "${assemblies_dir}/joblog.tsv" \
    --results "${assemblies_dir}/logs" --timeout 8h \
    < "${assemblies_dir}/jobs.txt"
parallel_exit_general=$?
set -e

if [[ ${parallel_exit_general} -ne 0 ]]; then
    if [[ ${parallel_exit_general} -ge 1 && ${parallel_exit_general} -le 254 ]]; then
        log_warn "${parallel_exit_general} general assembler job(s) failed. Surviving assemblies will continue downstream."
        awk -F'\t' 'NR>1 && $7!=0 {print "  -> Failed: " $9}' "${assemblies_dir}/joblog.tsv" >&2
    else
        log_error "GNU parallel failed catastrophically (exit code ${parallel_exit_general}). Aborting."
    fi
fi

# --- 4e. Apply weighting ---
log_info "4e. Applying Autocycler weighting..."
shopt -s nullglob
tmp_weight=$(mktemp)
for f in "${assemblies_dir}"/plassembler*.fasta; do
    sed '/Autocycler_cluster_weight=/! s/circular=[Tt][Rr][Uu][Ee]/circular=True Autocycler_cluster_weight=3/I' "$f" > "$tmp_weight" && mv "$tmp_weight" "$f"
done
for f in "${assemblies_dir}"/canu*.fasta "${assemblies_dir}"/flye*.fasta "${assemblies_dir}"/hifiasm*.fasta; do
    sed '/Autocycler_consensus_weight=/! s/^>.*$/& Autocycler_consensus_weight=2/' "$f" > "$tmp_weight" && mv "$tmp_weight" "$f"
done
run_cmd rm -f "$tmp_weight"
shopt -u nullglob

# --- 4f. Cleanup ---
run_cmd rm -f subsampled_reads/*.fastq

# ==============================================================================
# CURATION POINT 1 — Inspect assemblies
# ==============================================================================

failed_canu=$(awk -F'\t' 'NR>1 && $7!=0' "${assemblies_dir}/joblog_canu.tsv" 2>/dev/null | wc -l)
failed_general=$(awk -F'\t' 'NR>1 && $7!=0' "${assemblies_dir}/joblog.tsv" 2>/dev/null | wc -l)
empty_count=$(find "${assemblies_dir}" -name "*.fasta" -size 0 2>/dev/null | wc -l)

advice="DIAGNOSTICS:
  Failed jobs: Canu=${failed_canu}, General=${failed_general}
  Empty FASTA files: ${empty_count}

  Assembly sizes:"

for f in "${assemblies_dir}"/*.fasta; do
    [[ -f "$f" ]] && advice+="
    $(basename "$f"): $(grep -v '^>' "$f" | tr -d '\n' | wc -c) bp"
done

advice+="

ACTION: Delete empty or broken FASTA files, then continue."

manual_curation_pause "Inspect assemblies before clustering" "${advice}"

# --- 4h. Compress & Cluster ---
log_step "4h. Compress & cluster"
run_cmd autocycler compress -i "${assemblies_dir}" -a "${autocycler_dir}"
run_cmd autocycler cluster -a "${autocycler_dir}"

# ==============================================================================
# CURATION POINT 2 — Inspect clustering
# ==============================================================================

cluster_advice="CLUSTER SUMMARY:
$(column -t -s, "${autocycler_dir}/clustering/summary.csv" 2>/dev/null || cat "${autocycler_dir}/clustering/summary.csv" 2>/dev/null || echo "  (summary.csv not found)")

GUIDE:
  Large cluster (~genome size) = chromosome → qc_pass/
  Small clusters (1-200 kb)    = plasmids   → qc_pass/
  Tiny clusters (<500 bp)      = noise      → qc_fail/

ACTION: Move misclassified clusters between qc_pass/ and qc_fail/, then continue."

manual_curation_pause "Inspect clustering results" "${cluster_advice}"

# --- Hybrid dotplot function ---
# Small clusters (< 2 MB GFA): Autocycler dotplot (fast, all-vs-all pairwise)
# Large clusters (>= 2 MB GFA): MUMmer nucmer + mummerplot (chromosome-scale)
generate_dotplots() {
    local cluster_dir="$1"
    local gfa_untrimmed="${cluster_dir}/1_untrimmed.gfa"
    local gfa_trimmed="${cluster_dir}/2_trimmed.gfa"
    local gfa_size
    gfa_size=$(wc -c < "${gfa_untrimmed}")

    if [[ ${gfa_size} -lt 2000000 ]]; then
        # Plasmid-scale: Autocycler dotplot (purpose-built for trim verification)
        log_info "  Small cluster — using Autocycler dotplot"
        run_cmd autocycler dotplot -i "${gfa_untrimmed}" -o "${cluster_dir}/1_untrimmed.png"
        run_cmd autocycler dotplot -i "${gfa_trimmed}"   -o "${cluster_dir}/2_trimmed.png"
    elif command -v nucmer &>/dev/null && command -v mummerplot &>/dev/null; then
        # Chromosome-scale: MUMmer (efficient on multi-Mbp sequences)
        log_info "  Large cluster — using MUMmer (nucmer + mummerplot)"
        for tag in 1_untrimmed 2_trimmed; do
            local gfa="${cluster_dir}/${tag}.gfa"
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

# --- 4j. Trim & Resolve ---
log_step "4j. Trim & resolve"
for c in "${autocycler_dir}"/clustering/qc_pass/cluster_*; do
    log_info "Processing $(basename "$c")..."
    run_cmd autocycler trim -c "$c"
    generate_dotplots "$c"
    run_cmd autocycler resolve -c "$c"
done

# ==============================================================================
# CURATION POINT 3 — Inspect dotplots
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

# --- 4l. Combine ---
log_step "4l. Combining consensus assembly"
run_cmd autocycler combine -a "${autocycler_dir}" \
    -i "${autocycler_dir}"/clustering/qc_pass/cluster_*/5_final.gfa

run_cmd cp "${autocycler_dir}/consensus_assembly.fasta" "${autocycler_consensus}"

# --- 4o. Metrics ---
run_cmd bash -c 'autocycler table > "$1"' _ metrics.tsv
run_cmd bash -c 'autocycler table -a "$1" -n "$2" >> "$3"' \
    _ "${autocycler_dir}" "${sample_name}" metrics.tsv

log_step "Filtering & assembly complete."
log_info "Consensus: ${autocycler_consensus}"
log_info "Metrics:   metrics.tsv"
