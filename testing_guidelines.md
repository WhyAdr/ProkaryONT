# Pipeline Testing & Verification Guidelines

This document outlines how to safely test and dry-run the ProkaryONT assembly pipeline. Since the pipeline depends on a heavy stack of external command-line tools (e.g., `nanoplot`, `fastcat`, `seqkit`, `autocycler`), a mock testing utility is provided to simulate runs without installing these dependencies.

---

## 🛠️ The Testing Utility: `generate_mock_environment.py`

Located at [generate_mock_environment.py](file:///d:/W/ProkaryONT/utils/generate_mock_environment.py), this Python script automates the creation and cleanup of test fixtures, including fake executable binaries (which return exit code `0`) and valid gzipped sequence structures.

### 1. Setup the Mock Environment
Before testing, run the python script from the root of the repository to generate all necessary directories, mock files, and executables:
```bash
python utils/generate_mock_environment.py
```
This command:
*   Creates a `mock_bin/` folder populated with fake scripts for all pipeline dependencies.
*   Creates directories like `pod5_dir/` and `autocycler_out/`.
*   Generates a valid dummy FASTQ structure in `input.fastq.gz` and `filtered_input.fastq.gz`.
*   Touches intermediate files (such as `consensus_assembly.fasta` and `polished_assembly.fasta`) required by the sequential script validation checks (`require_file`).

### 2. Run Pipeline Dry-Run tests
To run the dry-run verification, add the `mock_bin/` folder to your path and execute the orchestration scripts in dry-run mode:

#### On Linux / macOS / Git Bash:
```bash
# Add mock bin folder to Path
export PATH="$(pwd)/mock_bin:$PATH"

# Run the complete orchestrator
bash run_all.sh --config pipeline.conf --dry-run --skip-curation
```

#### Running Standalone Steps:
You can also run specific scripts independently:
```bash
# Verify Autocycler assembly alone on a custom read path
bash 03_autocycler_assemble.sh --reads input.fastq.gz --dry-run --skip-curation
```

### 3. Clean up the Environment
Once testing is finished, remove all mock directories and files by running the clean flag:
```bash
python utils/generate_mock_environment.py --clean
```
This guarantees that no test artifacts (FASTQ, fasta, logs, or metrics tables) are left behind in the workspace.

### 4. Verify the SDUST assessment artifacts

Stage 1 profiles SDUST low-complexity burden without modifying the input FASTQ.
After a non-dry Stage 1 run, inspect `01_qc/fastcat_sdust/` and confirm that
the partition invariant holds:

```text
total_reads = zero_sdust_reads + sdust_positive_reads
```

The machine-readable outputs are `sdust_summary.tsv`,
`sdust_fraction_hist.tsv`, `sdust_positive_reads.tsv`, and
`sdust_high_burden_reads.tsv`. The reported masked fractions originate from
Fastcat stderr and are rounded by Fastcat; they must not be converted into
claimed exact masked-base counts.

---

## 🚫 Git Policies for Test Artifacts

To prevent local test junk from polluting the remote repository:
1.  **Never stage mock binaries**: The `mock_bin/` directory is for local dry-run testing only.
2.  **Ignored extensions**: The [.gitignore](file:///d:/W/ProkaryONT/.gitignore) configuration is set up to automatically block intermediate file types (like `*.log`, `*.bam`, `*.fastq.gz`, `*.fasta`, and directories like `autocycler_out/`, `assemblies/`, and `01_qc/`).
3.  **Untracked clean files**: Always verify `git status` before committing to confirm that only core script updates are staged.
