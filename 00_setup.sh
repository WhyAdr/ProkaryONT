#!/usr/bin/env bash
# Shared logging, validation, configuration, dry-run, and stage-contract helpers.

set -euo pipefail

_prokaryont_setup_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
stage_contract_py="${stage_contract_py:-${_prokaryont_setup_dir}/utils/stage_contract.py}"
_stage_python_bin="${PROKARYONT_PYTHON:-}"

# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------

_log_ts() { date "+%Y-%m-%d %H:%M:%S"; }

# Older stages expect pipeline.log by default.  It is initialized lazily so
# merely sourcing this file, and especially a strict dry-run, never writes.
log_file="${log_file:-$(pwd)/pipeline.log}"
_log_ready=""

_append_log() {
    local msg="$1"
    [[ -z "${log_file:-}" || -n "${dry_run:-}" ]] && return 0
    if [[ -z "${_log_ready}" ]]; then
        if ! touch "${log_file}" &>/dev/null; then
            log_file="/dev/null"
        fi
        _log_ready=true
    fi
    printf '%s\n' "${msg}" >> "${log_file}"
}

init_stage_logging() {
    log_file="$1"
    _log_ready=""
    [[ -n "${dry_run:-}" ]] && return 0
    mkdir -p "$(dirname "${log_file}")"
    if ! touch "${log_file}" &>/dev/null; then
        log_file="/dev/null"
    fi
    _log_ready=true
}

log_step() {
    local msg="[$(_log_ts)] [STEP] $*"
    printf '\033[1;32m%s\033[0m\n' "${msg}"
    _append_log "${msg}"
}

log_info() {
    local msg="[$(_log_ts)] [INFO] $*"
    printf '\033[0;36m%s\033[0m\n' "${msg}"
    _append_log "${msg}"
}

log_warn() {
    local msg="[$(_log_ts)] [WARN] $*"
    printf '\033[1;33m%s\033[0m\n' "${msg}" >&2
    _append_log "${msg}"
}

log_error() {
    local msg="[$(_log_ts)] [ERROR] $*"
    printf '\033[1;31m%s\033[0m\n' "${msg}" >&2
    _append_log "${msg}"
    exit 1
}

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------

require_tool() {
    local tool="$1"
    command -v "${tool}" &>/dev/null || log_error "'${tool}' not found in \$PATH. Please install it."
}

require_file() {
    local filepath="$1"
    local hint="${2:-}"
    if [[ ! -f "${filepath}" ]]; then
        local msg="Required file not found: ${filepath}"
        [[ -n "${hint}" ]] && msg+="\n       Hint: run '${hint}' first."
        log_error "${msg}"
    fi
}

require_dir() {
    local dirpath="$1"
    local hint="${2:-}"
    if [[ ! -d "${dirpath}" ]]; then
        local msg="Required directory not found: ${dirpath}"
        [[ -n "${hint}" ]] && msg+="\n       Hint: run '${hint}' first."
        log_error "${msg}"
    fi
}

require_arg() {
    local flag="$1"
    local value="${2:-}"
    [[ -n "${value}" ]] || log_error "Missing required argument: ${flag}. Use --help for usage."
}

require_option_value() {
    local flag="$1"
    local value="${2:-}"
    [[ -n "${value}" && "${value}" != --* ]] || log_error "${flag} requires a value. Use --help for usage."
}

require_positive_int() {
    local label="$1"
    local value="$2"
    [[ "${value}" =~ ^[1-9][0-9]*$ ]] || log_error "${label} must be a positive integer."
}

require_nonnegative_int() {
    local label="$1"
    local value="$2"
    [[ "${value}" =~ ^[0-9]+$ ]] || log_error "${label} must be a non-negative integer."
}

require_number_range() {
    local label="$1"
    local value="$2"
    local lower="$3"
    local upper="$4"
    local lower_inclusive="${5:-true}"
    local upper_inclusive="${6:-true}"
    [[ "${value}" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] || log_error "${label} must be numeric."
    awk -v value="${value}" -v lower="${lower}" -v upper="${upper}" \
        -v li="${lower_inclusive}" -v ui="${upper_inclusive}" \
        'BEGIN {
            lower_ok = (li == "true") ? value >= lower : value > lower
            upper_ok = (ui == "true") ? value <= upper : value < upper
            exit (lower_ok && upper_ok) ? 0 : 1
        }' || log_error "${label} must be in the range ${lower_inclusive} ${lower} to ${upper} ${upper_inclusive}."
}

require_boolean() {
    local label="$1"
    local value="$2"
    [[ "${value}" == "true" || "${value}" == "false" ]] || log_error "${label} must be true or false."
}

require_safe_sample_name() {
    local value="$1"
    [[ "${value}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || \
        log_error "--sample-name must contain only letters, digits, '.', '_', or '-', and start with a letter or digit."
}

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

load_config() {
    local config_file="$1"
    [[ -f "${config_file}" ]] || log_error "Config file not found: ${config_file}"
    log_info "Loading config from: ${config_file}"

    local line key value
    while IFS= read -r line || [[ -n "${line}" ]]; do
        [[ "${line}" =~ ^[[:space:]]*# ]] && continue
        [[ -z "${line//[[:space:]]/}" ]] && continue
        [[ "${line}" == *"="* ]] || {
            log_warn "Ignoring malformed configuration line without '=': ${line}"
            continue
        }
        key="${line%%=*}"
        value="${line#*=}"
        key="$(printf '%s' "${key}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
        value="$(printf '%s' "${value}" | sed 's/[[:space:]]#.*$//;s/^[[:space:]]*//;s/[[:space:]]*$//')"
        if [[ ( "${value}" == \"*\" && "${value}" == *\" ) || ( "${value}" == \'*\' && "${value}" == *\' ) ]]; then
            value="${value:1:${#value}-2}"
        fi
        if [[ ! "${key}" =~ ^[a-z][a-z0-9_]*$ ]]; then
            log_warn "Ignoring invalid or restricted configuration key: ${key}"
            continue
        fi
        if [[ -z "${!key:-}" ]]; then
            declare -g "${key}=${value}"
        fi
    done < "${config_file}"
}

# -----------------------------------------------------------------------------
# Manual curation and dry-run
# -----------------------------------------------------------------------------

manual_curation_pause() {
    local title="$1"
    local advice="$2"
    local skip="${skip_curation:-}"
    printf '\n========================================================================\n'
    log_warn "MANUAL CURATION POINT: ${title}"
    printf '========================================================================\n\n%s\n\n' "${advice}"
    if [[ "${skip}" == "true" ]]; then
        log_info "Curation skipped (--skip-curation is set)."
    else
        printf '%s\n' '------------------------------------------------------------------------'
        read -rp "Press Enter to continue (or Ctrl+C to abort)... "
        printf '\n'
    fi
}

dry_run="${dry_run:-}"

command_string() {
    local cmd_str
    cmd_str=$(printf "%q " "$@")
    printf '%s' "${cmd_str% }"
}

run_cmd() {
    local cmd_str
    cmd_str="$(command_string "$@")"
    if [[ -n "${dry_run:-}" ]]; then
        log_info "[DRY-RUN] ${cmd_str}"
    else
        log_info "[CMD] ${cmd_str}"
        "$@"
    fi
}

# -----------------------------------------------------------------------------
# Portable paths, identities, signatures, and completion manifests
# -----------------------------------------------------------------------------

stage_python() {
    if [[ -n "${_stage_python_bin}" ]]; then
        printf '%s' "${_stage_python_bin}"
        return
    fi
    if [[ -n "${MSYSTEM:-}" ]] && command -v python &>/dev/null && python -c 'import sys' &>/dev/null; then
        _stage_python_bin="python"
    elif command -v python3 &>/dev/null && python3 -c 'import sys' &>/dev/null; then
        _stage_python_bin="python3"
    elif command -v python &>/dev/null && python -c 'import sys' &>/dev/null; then
        _stage_python_bin="python"
    else
        log_error "A working Python 3 interpreter is required."
    fi
    printf '%s' "${_stage_python_bin}"
}

canonical_path() {
    local resolved
    resolved="$("$(stage_python)" "${stage_contract_py}" canonical "$1")"
    if command -v cygpath &>/dev/null && [[ "${resolved}" =~ ^[A-Za-z]:[\\/] ]]; then
        cygpath -u "${resolved}"
    else
        printf '%s' "${resolved}"
    fi
}

input_identity() {
    local path="$1"
    local hash_mode="${2:-false}"
    if [[ "${hash_mode}" == "true" ]]; then
        "$(stage_python)" "${stage_contract_py}" identity --sha256 "${path}"
    else
        "$(stage_python)" "${stage_contract_py}" identity "${path}"
    fi
}

stable_signature() {
    "$(stage_python)" "${stage_contract_py}" signature
}

validate_fastq_nonempty() {
    "$(stage_python)" "${stage_contract_py}" validate-fastq "$1" >/dev/null
}

tool_version() {
    local tool="$1"
    local output
    output="$("${tool}" --version 2>&1 || true)"
    [[ -n "${output}" ]] || output="$("${tool}" -V 2>&1 || true)"
    output="${output%%$'\n'*}"
    printf '%s' "${output:-version-unavailable}"
}

tool_path() {
    command -v "$1" 2>/dev/null || printf 'not-found'
}

manifest_signature() {
    local manifest="$1"
    awk -F '\t' '$1 == "signature" { print $2; exit }' "${manifest}"
}

validate_completion_manifest() {
    local manifest="$1"
    [[ -f "${manifest}" ]] || return 1
    grep -q $'^expected\t' "${manifest}" || return 1
    local kind path
    while IFS=$'\t' read -r kind path; do
        case "${kind}" in
            file) [[ -f "${path}" ]] || return 1 ;;
            nonempty) [[ -s "${path}" ]] || return 1 ;;
            dir) [[ -d "${path}" ]] || return 1 ;;
            fastq) validate_fastq_nonempty "${path}" || return 1 ;;
            *) return 1 ;;
        esac
    done < <(awk -F '\t' '$1 == "expected" { print $2 "\t" $3 }' "${manifest}")
}

write_completion_manifest() {
    local manifest="$1"
    local stage="$2"
    local signature="$3"
    local run_id="$4"
    shift 4
    local tmp_manifest="${manifest}.tmp.$$"
    {
        printf 'schema\tprokaryont.stage-completion.v1\n'
        printf 'stage\t%s\n' "${stage}"
        printf 'signature\t%s\n' "${signature}"
        printf 'run_id\t%s\n' "${run_id}"
        printf 'completed_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        local spec kind path key value
        for spec in "$@"; do
            if [[ "${spec}" == field:* ]]; then
                key="${spec#field:}"
                value="${key#*=}"
                key="${key%%=*}"
                printf 'field\t%s\t%s\n' "${key}" "${value}"
                continue
            fi
            kind="${spec%%:*}"
            path="${spec#*:}"
            printf 'expected\t%s\t%s\n' "${kind}" "${path}"
        done
    } > "${tmp_manifest}"
    mv -f -- "${tmp_manifest}" "${manifest}"
}

infer_sample_name() {
    local name
    name="$(basename "$1")"
    name="${name%.gz}"
    name="${name%.fastq}"
    name="${name%.fq}"
    name="${name//[^A-Za-z0-9._-]/_}"
    printf '%s' "${name:-sample}"
}

detect_threads() {
    local candidate="${SLURM_CPUS_PER_TASK:-}"
    if [[ "${candidate}" =~ ^[1-9][0-9]*$ ]]; then
        printf '%s' "${candidate}"
        return
    fi
    candidate="$(nproc 2>/dev/null || true)"
    if [[ "${candidate}" =~ ^[1-9][0-9]*$ ]]; then
        local cgroup_threads=""
        if [[ -r /sys/fs/cgroup/cpu.max ]]; then
            read -r quota period < /sys/fs/cgroup/cpu.max || true
            if [[ "${quota:-}" =~ ^[1-9][0-9]*$ && "${period:-}" =~ ^[1-9][0-9]*$ ]]; then
                cgroup_threads=$((quota / period))
                ((cgroup_threads >= 1)) || cgroup_threads=1
            fi
        elif [[ -r /sys/fs/cgroup/cpu/cpu.cfs_quota_us && -r /sys/fs/cgroup/cpu/cpu.cfs_period_us ]]; then
            quota="$(cat /sys/fs/cgroup/cpu/cpu.cfs_quota_us)"
            period="$(cat /sys/fs/cgroup/cpu/cpu.cfs_period_us)"
            if [[ "${quota}" =~ ^[1-9][0-9]*$ && "${period}" =~ ^[1-9][0-9]*$ ]]; then
                cgroup_threads=$((quota / period))
                ((cgroup_threads >= 1)) || cgroup_threads=1
            fi
        fi
        if [[ "${cgroup_threads}" =~ ^[1-9][0-9]*$ && "${cgroup_threads}" -lt "${candidate}" ]]; then
            candidate="${cgroup_threads}"
        fi
        printf '%s' "${candidate}"
    else
        printf '4'
    fi
}

detect_memory_gb() {
    local memory_mb=""
    if [[ "${SLURM_MEM_PER_NODE:-}" =~ ^[1-9][0-9]*$ ]]; then
        memory_mb="${SLURM_MEM_PER_NODE}"
    elif [[ "${SLURM_MEM_PER_CPU:-}" =~ ^[1-9][0-9]*$ && "${SLURM_CPUS_PER_TASK:-}" =~ ^[1-9][0-9]*$ ]]; then
        memory_mb=$((SLURM_MEM_PER_CPU * SLURM_CPUS_PER_TASK))
    elif [[ -r /proc/meminfo ]]; then
        memory_mb="$(awk '/^MemAvailable:/ { printf "%d", $2 / 1024; exit }' /proc/meminfo)"
    fi
    local cgroup_bytes=""
    if [[ -r /sys/fs/cgroup/memory.max ]]; then
        cgroup_bytes="$(cat /sys/fs/cgroup/memory.max)"
    elif [[ -r /sys/fs/cgroup/memory/memory.limit_in_bytes ]]; then
        cgroup_bytes="$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes)"
    fi
    if [[ "${cgroup_bytes}" =~ ^[1-9][0-9]*$ ]]; then
        local cgroup_mb=$((cgroup_bytes / 1024 / 1024))
        if ((cgroup_mb >= 1)) && { [[ ! "${memory_mb}" =~ ^[1-9][0-9]*$ ]] || ((cgroup_mb < memory_mb)); }; then
            memory_mb="${cgroup_mb}"
        fi
    fi
    [[ "${memory_mb}" =~ ^[1-9][0-9]*$ ]] || return 1
    # Bound Meryl to 70% of the visible/allotted memory and retain at least 1 GB.
    local memory_gb=$((memory_mb * 70 / 100 / 1024))
    ((memory_gb >= 1)) || memory_gb=1
    printf '%s' "${memory_gb}"
}

get_gzip_cmd() {
    if command -v pigz &>/dev/null; then
        printf 'pigz -p %s' "${threads:-4}"
    else
        printf 'gzip'
    fi
}
