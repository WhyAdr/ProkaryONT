#!/usr/bin/env bash
# ==============================================================================
# prokaryont — Single Orchestrator Entrypoint for the ProkaryONT Pipeline
# ==============================================================================
#
# This wrapper provides a subcommand interface over the pipeline scripts.
# It does NOT redefine logging — 00_setup.sh (sourced by each script) owns that.
#
# Communication with pipeline scripts uses exported PROKARYONT_* env vars.
# Child scripts inherit these automatically via bash's export mechanism:
#   PROKARYONT_RESUME     — resume target stage (cluster, trim, dnaapler)
#   PROKARYONT_BANDAGE    — whether to hint at Bandage inspection
#   PROKARYONT_RAGTAG     — reference FASTA for optional RagTag scaffolding
#   PROKARYONT_ASSEMBLERS — comma-separated assembler override list
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

# ==============================================================================
# USAGE
# ==============================================================================
usage() {
    cat << EOF
Usage: prokaryont <subcommand> [options]

Subcommands:
  assemble      Run the full assembly pipeline (Stages 1-3)
  taxonomy      Run taxonomy assignment (04_taxonomy.sh)
  annotate      Run Bakta annotation and QA (05_annotate_assess.sh)

Global Options:
  --config FILE         Path to pipeline.conf (required for most runs)
  --resume-from STAGE   Resume pipeline from a specific stage (cluster, trim, dnaapler)
                        Existing outputs at that stage will be overwritten.
  --ragtag REF_FASTA    Run RagTag scaffolding post-assembly
  --bandage             Hint at Bandage GFA inspection during curation points
  --assemblers LIST     Comma-separated list of assemblers to use
                        (default: flye,canu,hifiasm,raven,miniasm,metamdbg,
                         myloasm,plassembler,nextdenovo,wtdbg2)
  --help                Show this help message
EOF
    exit 0
}

if [[ $# -eq 0 ]]; then
    usage
fi

SUBCOMMAND=$1
shift

# --- Parse Global Options ---
CONFIG_FILE=""
RESUME_FROM=""
RAGTAG_REF=""
RUN_BANDAGE=false
ASSEMBLERS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)       CONFIG_FILE="$2"; shift 2 ;;
        --resume-from)  RESUME_FROM="$2"; shift 2 ;;
        --ragtag)       RAGTAG_REF="$2"; shift 2 ;;
        --bandage)      RUN_BANDAGE=true; shift ;;
        --assemblers)   ASSEMBLERS="$2"; shift 2 ;;
        --help|-h)      usage ;;
        *)              echo "ERROR: Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -n "${CONFIG_FILE}" && ! -f "${CONFIG_FILE}" ]]; then
    echo "ERROR: Config file not found: ${CONFIG_FILE}" >&2
    exit 1
fi

# Validate resume target
if [[ -n "${RESUME_FROM}" ]]; then
    case "${RESUME_FROM}" in
        cluster|trim|dnaapler) ;;
        *) echo "ERROR: Invalid --resume-from target '${RESUME_FROM}'. Valid: cluster, trim, dnaapler" >&2; exit 1 ;;
    esac
fi

# Export env vars — child scripts (sourcing 00_setup.sh) read these directly
export PROKARYONT_RESUME="${RESUME_FROM}"
export PROKARYONT_BANDAGE="${RUN_BANDAGE}"
export PROKARYONT_RAGTAG="${RAGTAG_REF}"
export PROKARYONT_ASSEMBLERS="${ASSEMBLERS}"

# ==============================================================================
# EXECUTION
# ==============================================================================

case "${SUBCOMMAND}" in
    assemble)
        [[ -z "${CONFIG_FILE}" ]] && { echo "ERROR: --config is required for 'assemble'" >&2; exit 1; }
        bash "${SCRIPT_DIR}/run_all.sh" --config "${CONFIG_FILE}"
        ;;
    taxonomy)
        [[ -z "${CONFIG_FILE}" ]] && { echo "ERROR: --config is required for 'taxonomy'" >&2; exit 1; }
        bash "${SCRIPT_DIR}/04_taxonomy.sh" --config "${CONFIG_FILE}"
        ;;
    annotate)
        [[ -z "${CONFIG_FILE}" ]] && { echo "ERROR: --config is required for 'annotate'" >&2; exit 1; }
        bash "${SCRIPT_DIR}/05_annotate_assess.sh" --config "${CONFIG_FILE}"
        ;;
    *)
        echo "ERROR: Unknown subcommand '${SUBCOMMAND}'. Use 'prokaryont --help'" >&2
        exit 1
        ;;
esac
