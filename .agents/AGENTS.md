# ProkaryONT Project Rules

Behavioral rules for editing, testing, and staging files in the ProkaryONT
bacterial assembly pipeline repository.

## 1. Pipeline Architecture

The pipeline is split into four core stages plus downstream analysis scripts:

| Stage | Script | Role |
|-------|--------|------|
| 1 | `01_qc_estimate.sh` | QC reports (NanoPlot, Fastcat, SNIKT) & genome size estimation |
| 2 | `02_preprocess_filter.sh` | Adapter trimming, end-cropping, low-complexity filtering, length/quality filtering |
| 3 | `03_autocycler_assemble.sh` | Multi-assembler subsampling, assembly, consensus resolution |
| 4 | `04_polish_orient.sh` | Dorado polishing & Dnaapler reorientation |

Downstream scripts are numbered `05_` through `07e_`. Keep numbering consistent
when adding or renaming scripts, and update cross-references in comments,
`pipeline.conf`, and `run_all.sh`.

`03_autocycler_assemble.sh` supports standalone execution via `--reads` or
`--input-fastq` to bypass Stage 2 with an arbitrary FASTQ file.

## 2. Chopper Trimming Rules

- **Use `--cutoff`** (the official Chopper flag) for all quality-aware trimming
  modes. Do NOT use `--quality`.
- **Supported `--trim-approach` values**: `fixed-crop`, `trim-by-quality`,
  `best-read-segment`, `split-by-low-quality`.
- When modifying the preprocessing stage, ensure `--chopper-trim-cutoff` and
  `--split-window` are forwarded through `run_all.sh` and documented in
  `pipeline.conf`.

## 3. Fastcat Lint Rules

- The subcommand is `fastcat lint` (also aliased as `fastcat fastlint`).
- Three parameters: `--threshold` (DUST T, default 20), `--window` (W, default
  64), `--max-proportion` (max masked fraction, default 0.95).
- Pipeline CLI flags are `--lint-threshold`, `--lint-window`,
  `--lint-max-proportion`.
- Rejected read headers are logged to `01_qc/02_lint_rejected.log`.

## 4. Dry-Run Verification

- **Setup**: `python utils/generate_mock_environment.py` creates mock binaries,
  dummy FASTQ files, and required directories.
- **Run**: `export PATH="$(pwd)/mock_bin:$PATH"` then execute scripts with
  `--dry-run`.
- **Cleanup**: `python utils/generate_mock_environment.py --clean` after every
  test session.
- When adding new tool dependencies to scripts, also add them to the `TOOLS`
  list in `utils/generate_mock_environment.py`.

## 5. Git Hygiene

- `mock_bin/` and `pod5_dir/` are not in `.gitignore` — never stage them.
- Do not stage draft design documents (e.g., `*_implementation_plan*.md`,
  `*_overview.md`).
- Prepend new entries to `changelogs.txt`. Only edit existing entries to correct
  factual errors.
- Always run `git status` before committing to verify only intended files are
  staged.
