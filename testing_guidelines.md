# Pipeline Testing & Verification Guidelines

This document outlines how to safely test and dry-run the ProkaryONT assembly pipeline. Since the pipeline depends on a heavy stack of external command-line tools (e.g., `nanoplot`, `fastcat`, `seqkit`, `autocycler`), a mock testing utility is provided to simulate runs without installing these dependencies.

---

## 🛠️ The Testing Utility: `generate_mock_environment.py`

Located at [generate_mock_environment.py](utils/generate_mock_environment.py), this Python script automates the creation and cleanup of test fixtures. The Stage 1/2 mocks validate supported subcommands and options and produce representative reports rather than returning success for arbitrary syntax.

### 1. Setup the Mock Environment
Before testing, run the python script from the root of the repository to generate all necessary directories, mock files, and executables:
```bash
python utils/generate_mock_environment.py
```
This command:
*   Creates a `mock_bin/` folder populated with fake scripts for all pipeline dependencies.
*   Creates directories like `pod5_dir/`, `autocycler_out/`, and `plassembler_db/`.
*   Generates a valid dummy FASTQ structure in `input.fastq.gz` and `filtered_input.fastq.gz`.
*   Touches intermediate files (such as `consensus_assembly.fasta` and `polished_assembly.fasta`) required by the sequential script validation checks (`require_file`).

### 2. Run modular dry-run tests
Add the `mock_bin/` folder and mock Plassembler database to your environment,
then invoke each active script separately. There is no active master
orchestrator.

#### On Linux / macOS / Git Bash:
```bash
# Add mock bin folder to Path
export PATH="$(pwd)/mock_bin:$PATH"
export PLASSEMBLER_DB="$(pwd)/plassembler_db"

# Examples: validate the changed active-stage command contracts
rm -f filtered_input.fastq.gz
bash 01_qc_estimate.sh \
  --input-fastq input.fastq.gz \
  --memory 1 \
  --dry-run
bash 02_preprocess_filter.sh --input-fastq input.fastq.gz --dry-run > no_keep.log 2>&1
! grep -q -- '--keep_percent' no_keep.log
bash 02_preprocess_filter.sh --input-fastq input.fastq.gz --keep-percent 90 --dry-run > keep90.log 2>&1
grep -q -- '--keep_percent 90' keep90.log
bash 03_autocycler_assemble.sh \
  --reads input.fastq.gz \
  --read-type ont_r10 \
  --genome-size 5000000 \
  --assemblers flye,miniasm,myloasm \
  --subsample-count 2 \
  --compress-threads 100 \
  --combine-threads 100 \
  --skip-curation \
  --dry-run
bash 04_polish_orient.sh \
  --assembly autocycler_out/consensus_assembly.fasta \
  --pod5-dir pod5_dir/ \
  --device cuda:1 \
  --read-group sample_rg \
  --double-polish \
  --skip-curation \
  --dry-run
```

The Stage 3 dry-run must show both `autocycler compress` and `autocycler
cluster` with the same `--max_contigs` placeholder. Because dry-run mode does
not count real assembler FASTAs, it deliberately labels the effective value
`DRY_RUN_UNAVAILABLE`.

The Stage 4 dry-run must show two `dorado polish` commands containing
`--RG sample_rg` when `--double-polish` and `--read-group sample_rg` are used.
Header membership is validated only in a real run after each BAM is created;
dry-run checks validate the generated command contract.

### 3. Run deterministic contig-policy tests

The standard-library test suites exercise strict FASTA/PAF parsing,
single-target containment, exact equal-length duplicate handling, conditional
minimap2 invocation, prepared-copy preservation, guard decisions, successful
PAF deletion, and fingerprint reuse:

```bash
python3 -m unittest discover -s tests -v
python3 -m py_compile dedup_contained.py tests/test_dedup_contained.py \
  tests/test_stage3_contig_policy.py tests/test_documentation_cli_contracts.py \
  tests/test_stage1_stage2_hardening.py utils/generate_mock_environment.py \
  utils/stage_contract.py
bash -n 01_qc_estimate.sh 02_preprocess_filter.sh \
  03_autocycler_assemble.sh 04_polish_orient.sh 05_taxonomy.sh
```

`utils/generate_mock_environment.py` also creates
`mock_contig_policy/{normal,fragmented}.fasta` and a representative
`contained.paf` for manual inspection. These fixtures and the mock minimap2 are
structural tests only; they are not biological validation.

For a successful fragmented-assembly fixture, confirm that
`assembly_assessment/dedup/*.events.tsv`, `*.minimap2.log`, and `*.dedup.log`
remain while the corresponding `*.paf` has been deleted. On failure, inspect
the reported `assembly_assessment.tmp.<pid>/dedup/` directory before cleanup.

### 4. Clean up the Environment
Once testing is finished, remove all mock directories and files by running the clean flag:
```bash
python utils/generate_mock_environment.py --clean
```
This guarantees that no test artifacts (FASTQ, fasta, logs, or metrics tables) are left behind in the workspace.

### 5. Verify the SDUST assessment artifacts

Stage 1 profiles SDUST low-complexity burden without modifying the input FASTQ.
After a non-dry Stage 1 run, inspect `01_qc/fastcat_sdust/` and confirm that
the partition invariant holds:

```text
total_reads = zero_sdust_reads + sdust_positive_reads
```

The machine-readable outputs are `sdust_summary.tsv`,
`sdust_fraction_hist.tsv`, `sdust_positive_reads.tsv`, and
`sdust_high_burden_reads.tsv`. The histogram distinguishes `exact_zero` from
`positive_reported` values that Fastcat rounded to `0.00`. Reported masked
fractions originate from Fastcat stderr and must not be converted into claimed
exact masked-base counts. `sdust_summary.tsv` separately records the rejection
count from a count-only pass at the exact configured Stage 2 threshold.

### 6. Verify strict dry-run and completion contracts

For Stages 1 and 2, snapshot the target directory before and after `--dry-run`;
the file tree must be identical. A real run writes its completion manifest
only after every expected output validates. Repeating an identical run should
resume immediately, while changing an input metadata field, effective option,
or tool version must fail closed until `--force` is supplied. Run
`tests/test_stage1_stage2_hardening.py` for deterministic coverage of these
contracts, input/output alias rejection, metadata-versus-SHA identity, atomic
FASTQ validation, and every optional Stage 2 dry-run branch.

---

## 🚫 Git Policies for Test Artifacts

To prevent local test junk from polluting the remote repository:
1.  **Never stage mock binaries**: The `mock_bin/` directory is for local dry-run testing only.
2.  **Ignored extensions**: The [.gitignore](.gitignore) configuration is set up to automatically block intermediate file types (like `*.log`, `*.bam`, `*.fastq.gz`, `*.fasta`, and directories like `autocycler_out/`, `assemblies/`, and `01_qc/`).
3.  **Untracked clean files**: Always verify `git status` before committing to confirm that only core script updates are staged.
