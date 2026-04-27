#!/usr/bin/env bash
# ==============================================================================
# 00_setup.sh — Shared function library for the pipeline
# ==============================================================================
# Source this file from every other script:
#   source "$(dirname "$0")/00_setup.sh"
#
# Provides: logging, validation, config loading, manual curation, dry-run.
# Does NOT set any user-configurable variables — those come from CLI flags
# or the config file.
# ==============================================================================

set -euo pipefail

# ==============================================================================
# LOGGING
# ==============================================================================

_log_ts() { date "+%Y-%m-%d %H:%M:%S"; }

log_step() {
    local msg="[$(_log_ts)] [STEP] $*"
    echo -e "\033[1;32m${msg}\033[0m"
    [[ -n "${log_file:-}" ]] && echo "${msg}" >> "${log_file}"
}

log_info() {
    local msg="[$(_log_ts)] [INFO] $*"
    echo -e "\033[0;36m${msg}\033[0m"
    [[ -n "${log_file:-}" ]] && echo "${msg}" >> "${log_file}"
}

log_warn() {
    local msg="[$(_log_ts)] [WARN] $*"
    echo -e "\033[1;33m${msg}\033[0m" >&2
    [[ -n "${log_file:-}" ]] && echo "${msg}" >> "${log_file}"
}

log_error() {
    local msg="[$(_log_ts)] [ERROR] $*"
    echo -e "\033[1;31m${msg}\033[0m" >&2
    [[ -n "${log_file:-}" ]] && echo "${msg}" >> "${log_file}"
    return 1 2>/dev/null || exit 1
}

# ==============================================================================
# VALIDATION
# ==============================================================================

require_tool() {
    local tool="$1"
    if ! command -v "${tool}" &>/dev/null; then
        log_error "'${tool}' not found in \$PATH. Please install it."
    fi
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
    # Usage: require_arg "--flag-name" "$value"
    local flag="$1"
    local value="$2"
    if [[ -z "${value}" ]]; then
        log_error "Missing required argument: ${flag}. Use --help for usage."
    fi
}

# ==============================================================================
# CONFIG FILE LOADING
# ==============================================================================

load_config() {
    # Usage: load_config "/path/to/pipeline.conf"
    # Reads key=value pairs, ignoring comments (#) and blank lines.
    # Only sets variables that are currently unset or empty, so CLI flags
    # (parsed before this call) take priority.
    local config_file="$1"
    if [[ ! -f "${config_file}" ]]; then
        log_error "Config file not found: ${config_file}"
    fi
    log_info "Loading config from: ${config_file}"

    while IFS= read -r line || [[ -n "$line" ]]; do
        # Skip comments and blank lines
        [[ "$line" =~ ^[[:space:]]*# ]] && continue
        [[ -z "${line//[[:space:]]/}" ]] && continue

        # Extract key=value
        local key="${line%%=*}"
        local value="${line#*=}"

        # Remove leading/trailing whitespace, inline comments, and quotes
        key="$(echo "${key}" | xargs)"
        value="$(echo "${value}" | sed 's/[[:space:]]#.*$//')"
        value="$(echo "${value}" | sed 's/^["'\'']*//;s/["'\'']*$//' | xargs)"

        # Only set if not already set (CLI flags take priority)
        if [[ -z "${!key:-}" ]]; then
            declare -g "${key}=${value}"
        fi
    done < "${config_file}"
}

# ==============================================================================
# MANUAL CURATION
# ==============================================================================

manual_curation_pause() {
    local title="$1"
    local advice="$2"
    local skip="${skip_curation:-}"

    echo ""
    echo "========================================================================"
    log_warn "MANUAL CURATION POINT: ${title}"
    echo "========================================================================"
    echo ""
    echo "${advice}"
    echo ""

    if [[ -n "${skip}" ]]; then
        log_info "Curation skipped (--skip-curation is set)."
    else
        echo "------------------------------------------------------------------------"
        read -rp "Press Enter to continue (or Ctrl+C to abort)... "
        echo ""
    fi
}

# ==============================================================================
# DRY-RUN
# ==============================================================================

dry_run="${dry_run:-}"

run_cmd() {
    if [[ -n "${dry_run}" ]]; then
        log_info "[DRY-RUN] $*"
    else
        "$@"
    fi
}

# ==============================================================================
# INIT
# ==============================================================================

log_file="${log_file:-$(pwd)/pipeline.log}"
touch "${log_file}"

# ==============================================================================
# COMPRESSION UTILITIES
# ==============================================================================

if command -v pigz &>/dev/null; then
    export GZIP_BIN="pigz -p ${threads:-4}"
else
    export GZIP_BIN="gzip"
fi
