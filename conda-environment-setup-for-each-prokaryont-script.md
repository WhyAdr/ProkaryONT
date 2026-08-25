# Conda Environment Setup for Each ProkaryONT Stage

Tested, modular Conda/Mamba environment specifications for every active
ProkaryONT stage, audited against the `require_tool` and `command -v` checks
in the current `main` branch scripts.

---

## 1. Environment Strategy

### Why per-stage environments?

ProkaryONT's five stages span C/C++/Rust binaries, R graphing libraries,
Java (Canu), Perl (NECAT, barrnap), Python 3.10+, and CUDA GPU runtimes
(Dorado).  A single monolithic environment will hit unresolvable solver
conflicts on most HPC systems.  Each environment below is independently
solvable and tested.

### Channel configuration (run once)

```bash
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels nanoporetech   # required for fastcat
conda config --set channel_priority strict
```

> [!IMPORTANT]
> **`fastcat` is distributed on the `nanoporetech` channel**, not on
> `bioconda`.  Every environment that includes `fastcat` must list
> `nanoporetech` as a channel or the solver will fail silently.

### Solver recommendation

Use [`mamba`](https://mamba.readthedocs.io/) (or `micromamba`) instead of
`conda` for all environment creation commands below.

---

## 2. Ground-truth dependency map

The table below is derived mechanically from `require_tool` calls (hard
failures) and `command -v` checks (soft warnings) in each script.

| Tool | Binary name checked | Stage 1 | Stage 2 | Stage 3 | Stage 4 | Stage 5 | Notes |
|---|---|:---:|:---:|:---:|:---:|:---:|---|
| NanoPlot | `NanoPlot` | ● | ● | | | | Hard requirement in both |
| fastcat | `fastcat` | ● | ◐ | | | | Hard in S1; conditional on `--enable-fastcat-lint` in S2 |
| Porechop_ABI | `porechop_abi` | ● | ◐ | | | | Hard in S1; conditional on `--enable-porechopabi-trim` in S2 |
| SNIKT | `snikt.R` | ● | ● | | | | Hard requirement (Bioconda may install as `snikt`; see symlink note in §3) |
| LRGE | `lrge` | ● | | | | | |
| Raven | `raven` | ● | | | | | |
| Meryl | `meryl` | ● | | | | | |
| Filtlong | `filtlong` | | ● | | | | |
| SeqKit | `seqkit` | | ● | ● | | | Hard in S2 and S3 |
| Fastplong | `fastplong` | | ● | | | | |
| Chopper | `chopper` | | ◐ | | | | Conditional on `--enable-chopper-trim` |
| Autocycler | `autocycler` | | | ● | | | |
| GNU Parallel | `parallel` | | | ● | | | |
| minimap2 | `minimap2` | | | ● | | | Also needed for iLesta & Miniasm helpers |
| samtools | `samtools` | | | ● | ● | | |
| Rasusa | `rasusa` | | | ◐ | | | Conditional on `--subsampler rasusa` |
| Plassembler | `plassembler` | | | ◐ | | | Conditional unless `--skip-plsdb-screen` |
| nucmer | `nucmer` | | | ◦ | | | Soft; needed for dotplots |
| mummerplot | `mummerplot` | | | ◦ | | | Soft; needs gnuplot at runtime |
| Dorado | `dorado` | | | | ● | | Not conda; see §6 |
| Dnaapler | `dnaapler` | | | | ● | | |
| GTDB-Tk | `gtdbtk` | | | | | ● | |
| MLST | `mlst` | | | | | ◦ | Soft; skipped if missing |
| Barrnap | `barrnap` | | | | | ◦ | Soft; skipped if missing |
| pigz | `pigz` | ○ | ○ | ○ | ○ | | Falls back to gzip if missing |

**Legend:** ● hard requirement  ◐ conditional hard  ◦ soft/optional  ○ implicit via `get_gzip_cmd`

### Per-assembler tool dependencies (Stage 3 only)

These are checked with `command -v` warnings, not hard failures.
Install only the ones matching your `--assemblers` list.

| Assembler | Binary(ies) checked | Bioconda package | Notes |
|---|---|---|---|
| flye | `flye` | `flye` | |
| canu | `canu` | `canu` | Java (`openjdk>=11`) required at runtime |
| hifiasm | `hifiasm` | `hifiasm` | |
| raven | `raven` | `raven-assembler` | |
| miniasm | `miniasm`, `minipolish`, `minimap2`, `racon` | `miniasm`, `minipolish`, `racon` | |
| metamdbg | `metaMDBG` | `metamdbg` | Note: lowercase conda package, uppercase binary |
| myloasm | `myloasm` | `myloasm` | |
| lja | `lja` | `lja` | |
| redbean (wtdbg2) | `wtdbg2` | `wtdbg` | `wtdbg2` is the binary; `wtdbg` is the conda package |
| necat | `necat` | `necat` | Perl-based; version `0.0.1_update20200803` on bioconda |
| nextdenovo | `nextdenovo` | `nextdenovo` | Bioconda installs `nextDenovo` (capital D); see symlink note in §5 |
| plassembler | `plassembler` | `plassembler` | Also an assembler helper, separate from PLSDB screen |
| ilesta | `Ilesta`, `minipolish`, `minimap2`, `racon` | **not on bioconda** | See iLesta note below |

> [!WARNING]
> **iLesta is not available as a bioconda package.** Install it manually from
> source (GitHub) and ensure the `Ilesta` binary is on `$PATH`.  Its
> additional dependencies (`minipolish`, `minimap2`, `racon`) are available
> on bioconda and should be installed in the Stage 3 environment.

---

## 3. Stage 1 — QC & Genome Size Estimation

**Script:** [`01_qc_estimate.sh`](file:///d:/W/ProkaryONT/01_qc_estimate.sh)

### `env_stage1_qc.yml`

```yaml
name: prokaryont-stage1-qc
channels:
  - conda-forge
  - bioconda
  - nanoporetech          # required for fastcat
  - defaults
dependencies:
  - python>=3.10,<3.12
  - pigz
  # --- QC & diagnostics ---
  - nanoplot>=1.42.0
  - fastcat>=0.18.0       # from nanoporetech channel
  - porechop_abi>=0.5.0
  - snikt                 # pulls in R, seqtk, and R libraries automatically
  # --- Genome size estimation ---
  - lrge>=0.1.0
  - raven-assembler>=1.8.3
  - meryl>=1.4.1
```

> [!IMPORTANT]
> **`snikt.R` Binary Name Compatibility**: `01_qc_estimate.sh` checks for
> `snikt.R` via `require_tool snikt.R`.  Bioconda packages sometimes install the
> binary entrypoint as `snikt` instead of `snikt.R`.  If `command -v snikt.R`
> fails after activating your environment, create a symlink in your conda `bin/`:
> ```bash
> [[ ! -f "${CONDA_PREFIX}/bin/snikt.R" && -f "${CONDA_PREFIX}/bin/snikt" ]] && \
>   ln -s "${CONDA_PREFIX}/bin/snikt" "${CONDA_PREFIX}/bin/snikt.R"
> ```

> [!NOTE]
> The bioconda `snikt` package automatically resolves its R dependencies
> (`r-tidyverse`, `r-gridextra`, `r-docopt`, `r-lubridate`) and the
> system tool `seqtk`.  Do not manually list individual R packages unless
> you are installing SNIKT from source.

> [!NOTE]
> Stage 1 does **not** use `seqkit`, `filtlong`, or `chopper`.  Stage 1 only
> requires `NanoPlot`, `fastcat`, `porechop_abi`, `snikt.R`, `lrge`, `raven`, and `meryl`.

### Verification

```bash
mamba env create -f env_stage1_qc.yml
conda activate prokaryont-stage1-qc

# Ensure snikt.R is available as an executable name
[[ ! -f "${CONDA_PREFIX}/bin/snikt.R" && -f "${CONDA_PREFIX}/bin/snikt" ]] && \
  ln -s "${CONDA_PREFIX}/bin/snikt" "${CONDA_PREFIX}/bin/snikt.R"

NanoPlot --version
fastcat --version
porechop_abi --version
snikt.R --help 2>&1 | head -1
lrge --help 2>&1 | head -1
raven --version
meryl --version

bash 01_qc_estimate.sh --help
```

---

## 4. Stage 2 — Preprocessing & Filtering

**Script:** [`02_preprocess_filter.sh`](file:///d:/W/ProkaryONT/02_preprocess_filter.sh)

Stages 1 and 2 share significant overlap (`NanoPlot`, `porechop_abi`,
`snikt.R`, `fastcat`) and can be merged into one environment if preferred.

### `env_stage2_preproc.yml`

```yaml
name: prokaryont-stage2-preproc
channels:
  - conda-forge
  - bioconda
  - nanoporetech          # required for fastcat
  - defaults
dependencies:
  - python>=3.10,<3.12
  - pigz
  # --- Always required ---
  - filtlong>=0.2.1
  - seqkit>=2.8.0
  - nanoplot>=1.42.0
  - snikt                 # unconditional require_tool in script
  - fastplong>=0.2.2
  # --- Conditional (install all; cheaper than debugging) ---
  - porechop_abi>=0.5.0   # required when --enable-porechopabi-trim
  - chopper>=0.7.0        # required when --enable-chopper-trim
  - fastcat>=0.18.0       # required when --enable-fastcat-lint
```

> [!TIP]
> Even though `porechop_abi`, `chopper`, and `fastcat` are conditional on
> opt-in CLI flags, install all three.  The `require_tool` check runs at
> script startup *after* flag parsing, so missing any enabled tool causes
> an immediate abort with no partial output.

### Verification

```bash
mamba env create -f env_stage2_preproc.yml
conda activate prokaryont-stage2-preproc

# Ensure snikt.R is available as an executable name
[[ ! -f "${CONDA_PREFIX}/bin/snikt.R" && -f "${CONDA_PREFIX}/bin/snikt" ]] && \
  ln -s "${CONDA_PREFIX}/bin/snikt" "${CONDA_PREFIX}/bin/snikt.R"

filtlong --version
seqkit version
chopper --version
fastplong --version
NanoPlot --version
snikt.R --help 2>&1 | head -1

bash 02_preprocess_filter.sh --help
```

---

## 5. Stage 3 — Multi-Assembler & Autocycler

**Script:** [`03_autocycler_assemble.sh`](file:///d:/W/ProkaryONT/03_autocycler_assemble.sh)

This is the largest environment.  The assembler list is user-selectable
via `--assemblers`, so only install what you intend to run.  The
**core workflow utilities** are always required.

### `env_stage3_assemble.yml`

```yaml
name: prokaryont-stage3-assemble
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python>=3.10,<3.12
  - pigz
  # --- Core workflow (always required) ---
  - autocycler>=0.5.0
  - parallel>=20230122
  - seqkit>=2.8.0
  - minimap2>=2.28
  - samtools>=1.19
  # --- Subsampling (install both; rasusa is small) ---
  - rasusa>=0.8.0
  # --- PLSDB characterization ---
  - plassembler>=1.6.0
  # --- Dotplot utilities (optional, soft-checked) ---
  - mummer4>=4.0.0rc1
  - gnuplot                # runtime dep of mummerplot
  # --- Runtime for Java-based tools (Canu) ---
  - openjdk>=11            # prevents JVM crashes during Canu overlap/gatekeeper jobs
  # --- Helper assemblers (install only what you use) ---
  - flye>=2.9.3
  - canu>=2.2
  - hifiasm>=0.19.8
  - raven-assembler>=1.8.3
  - miniasm>=0.3_r179
  - minipolish>=0.1.3
  - racon>=1.5.0
  - metamdbg>=1.1
  - myloasm>=0.1.0
  - lja>=0.2
  - wtdbg>=2.5
  - necat>=0.0.1
  - nextdenovo>=2.5.2
  # NOTE: iLesta (Ilesta binary) is not on bioconda.
  # Install from source and ensure Ilesta is on $PATH.
  # Its deps (minipolish, minimap2, racon) are listed above.
```

> [!IMPORTANT]
> **`nextdenovo` Case Sensitivity on Linux**: Bioconda packages install the
> binary as `nextDenovo` (with a capital `D`).  Stage 3's pre-flight loop checks
> `command -v nextdenovo`.  On case-sensitive Linux filesystems, create a symlink:
> ```bash
> [[ ! -f "${CONDA_PREFIX}/bin/nextdenovo" && -f "${CONDA_PREFIX}/bin/nextDenovo" ]] && \
>   ln -s "${CONDA_PREFIX}/bin/nextDenovo" "${CONDA_PREFIX}/bin/nextdenovo"
> ```

> [!IMPORTANT]
> **Plassembler PLSDB database** must be downloaded separately.  The
> script's database search order is:
> 1. `--plassembler-db` CLI argument
> 2. `plassembler_db` in `pipeline.conf`
> 3. `$PLASSEMBLER_DB` environment variable
> 4. `$CONDA_PREFIX/plassembler_db` directory
>
> ```bash
> conda activate prokaryont-stage3-assemble
> plassembler download -d /databases/plassembler_db
> export PLASSEMBLER_DB=/databases/plassembler_db
> ```

> [!NOTE]
> `edlib` is statically linked into the Autocycler binary and does not
> need to be installed separately via conda.

### Solver tips for large environments

If the solver cannot resolve all assemblers simultaneously:

```bash
# Option A: Split assemblers into a separate overlay environment
mamba create -n prokaryont-s3-core \
  -c conda-forge -c bioconda \
  autocycler parallel seqkit minimap2 samtools rasusa plassembler mummer4 gnuplot openjdk=11

mamba create -n prokaryont-s3-assemblers \
  -c conda-forge -c bioconda \
  flye canu hifiasm raven-assembler miniasm minipolish racon \
  metamdbg myloasm lja wtdbg necat nextdenovo

# Option B: Use --assemblers to limit the set and install only what you need
```

### Verification

```bash
mamba env create -f env_stage3_assemble.yml
conda activate prokaryont-stage3-assemble

# Create nextdenovo lowercase symlink if needed
[[ ! -f "${CONDA_PREFIX}/bin/nextdenovo" && -f "${CONDA_PREFIX}/bin/nextDenovo" ]] && \
  ln -s "${CONDA_PREFIX}/bin/nextDenovo" "${CONDA_PREFIX}/bin/nextdenovo"

autocycler --version
parallel --version | head -1
seqkit version
minimap2 --version
samtools --version | head -1
rasusa --version
plassembler --version
java -version 2>&1 | head -1

# Verify selected assemblers
for asm in flye canu hifiasm raven miniasm wtdbg2 necat nextdenovo; do
  command -v "$asm" && echo "  $asm: OK" || echo "  $asm: MISSING"
done

bash 03_autocycler_assemble.sh --help
```

---

## 6. Stage 4 — Dorado Polishing & Dnaapler

**Script:** [`04_polish_orient.sh`](file:///d:/W/ProkaryONT/04_polish_orient.sh)

### Dorado: Install from ONT releases, not conda

Dorado bundles its own CUDA runtime and cuDNN libraries.  Installing it
from conda will almost certainly break GPU compatibility with your HPC
driver stack.  Always use the official standalone binary:

```bash
# 1. Download (check https://github.com/nanoporetech/dorado/releases for latest)
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-<VERSION>-linux-x64.tar.gz
tar -xzf dorado-<VERSION>-linux-x64.tar.gz
export PATH="$(pwd)/dorado-<VERSION>-linux-x64/bin:$PATH"

# 2. Download basecalling models
dorado download --model sup
# Or a specific chemistry model:
# dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0

# 3. Verify GPU access
nvidia-smi
dorado basecaller --help
dorado polish --help
```

### `env_stage4_polish.yml`

```yaml
name: prokaryont-stage4-polish
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python>=3.10,<3.12
  - samtools>=1.19
  - pigz
  - dnaapler>=0.7.0
  # Dnaapler runtime dependencies:
  - blast>=2.14.0         # needed by dnaapler for start-gene identification
  - pyrodigal>=3.0        # needed by dnaapler for gene prediction
```

> [!WARNING]
> Dorado is **not** listed in the YAML.  It must be available on `$PATH`
> from the standalone install above.  The script's `require_tool dorado`
> check will fail immediately if it is missing.

### Verification

```bash
mamba env create -f env_stage4_polish.yml
conda activate prokaryont-stage4-polish

dorado --version          # from standalone install
samtools --version | head -1
dnaapler --version

bash 04_polish_orient.sh --help
```

---

## 7. Stage 5 — Taxonomy Classification

**Script:** [`05_taxonomy.sh`](file:///d:/W/ProkaryONT/05_taxonomy.sh)

### `env_stage5_taxonomy.yml`

```yaml
name: prokaryont-stage5-taxonomy
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  # GTDB-Tk pins its own Python version; let the solver decide
  - gtdbtk>=2.3.2
  # Optional but recommended (script skips gracefully if missing)
  - mlst>=2.23.0
  - barrnap>=0.9
  # GTDB-Tk transitive dependencies (usually pulled automatically):
  - hmmer>=3.3.2
  - prodigal>=2.6.3
  - fastani>=1.33
  - mash>=2.3
  - pplacer
```

> [!IMPORTANT]
> **GTDB-Tk reference database** (~85 GB decompressed) must be downloaded
> separately:
> ```bash
> conda activate prokaryont-stage5-taxonomy
>
> # Use the bundled download helper:
> download-db.sh /databases/gtdbtk_r220
>
> # Point the pipeline to it (set in pipeline.conf or shell profile):
> export GTDBTK_DATA_PATH=/databases/gtdbtk_r220
> ```

> [!TIP]
> GTDB-Tk v2.3+ requires Python ≤3.10 due to `pplacer` constraints.
> Do not force `python>=3.11` in this environment.

### Verification

```bash
mamba env create -f env_stage5_taxonomy.yml
conda activate prokaryont-stage5-taxonomy

gtdbtk --version
mlst --version
barrnap --version

bash 05_taxonomy.sh --help
```

---

## 8. Combined Quick-Setup Script

Save as `setup_all_envs.sh` and run on your HPC login node:

```bash
#!/usr/bin/env bash
set -euo pipefail

command -v mamba &>/dev/null && PM="mamba" || PM="conda"
echo "Using: ${PM}"

# --- Stage 1: QC & Genome Size ---
${PM} create -y -n prokaryont-stage1-qc \
  -c conda-forge -c bioconda -c nanoporetech \
  python=3.10 pigz nanoplot fastcat porechop_abi snikt lrge raven-assembler meryl

# Post-install symlink for snikt.R in Stage 1
_s1_prefix="$(conda info --base)/envs/prokaryont-stage1-qc"
[[ -f "${_s1_prefix}/bin/snikt" && ! -f "${_s1_prefix}/bin/snikt.R" ]] && \
  ln -s "${_s1_prefix}/bin/snikt" "${_s1_prefix}/bin/snikt.R"

# --- Stage 2: Preprocess & Filter ---
${PM} create -y -n prokaryont-stage2-preproc \
  -c conda-forge -c bioconda -c nanoporetech \
  python=3.10 pigz filtlong seqkit nanoplot snikt fastplong \
  porechop_abi chopper fastcat

# Post-install symlink for snikt.R in Stage 2
_s2_prefix="$(conda info --base)/envs/prokaryont-stage2-preproc"
[[ -f "${_s2_prefix}/bin/snikt" && ! -f "${_s2_prefix}/bin/snikt.R" ]] && \
  ln -s "${_s2_prefix}/bin/snikt" "${_s2_prefix}/bin/snikt.R"

# --- Stage 3: Autocycler Assembly ---
# Core utilities + Java + Assemblers:
${PM} create -y -n prokaryont-stage3-assemble \
  -c conda-forge -c bioconda \
  python=3.10 pigz autocycler parallel seqkit minimap2 samtools \
  rasusa plassembler mummer4 gnuplot openjdk=11 \
  flye canu hifiasm raven-assembler miniasm minipolish racon \
  metamdbg myloasm lja wtdbg necat nextdenovo

# Post-install symlink for nextdenovo lowercase alias
_s3_prefix="$(conda info --base)/envs/prokaryont-stage3-assemble"
[[ -f "${_s3_prefix}/bin/nextDenovo" && ! -f "${_s3_prefix}/bin/nextdenovo" ]] && \
  ln -s "${_s3_prefix}/bin/nextDenovo" "${_s3_prefix}/bin/nextdenovo"

# --- Stage 4: Polish & Orient ---
# NOTE: Dorado must be installed separately from ONT releases.
${PM} create -y -n prokaryont-stage4-polish \
  -c conda-forge -c bioconda \
  python=3.10 samtools pigz dnaapler blast pyrodigal

# --- Stage 5: Taxonomy ---
${PM} create -y -n prokaryont-stage5-taxonomy \
  -c conda-forge -c bioconda \
  gtdbtk mlst barrnap hmmer prodigal fastani mash pplacer

echo "=== All environments created ==="
echo ""
echo "Remaining manual steps:"
echo "  1. Install Dorado from https://github.com/nanoporetech/dorado/releases"
echo "  2. Install iLesta from source if using --assemblers ilesta"
echo "  3. Download Plassembler PLSDB:  plassembler download -d /path/to/db"
echo "  4. Download GTDB-Tk database:   download-db.sh /path/to/gtdbtk_db"
```

---

## 9. HPC / Slurm Integration

### Stage 3 example (compute node)

```bash
#!/bin/bash
#SBATCH --job-name=prokaryont_s3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --partition=compute

source ~/.bashrc
conda activate prokaryont-stage3-assemble

export PLASSEMBLER_DB=/databases/plassembler_db

bash 03_autocycler_assemble.sh \
  --config pipeline.conf \
  --reads filtered_input.fastq.gz \
  --read-type ont_r10 \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --assemblers flye,canu,hifiasm,raven,miniasm,lja,redbean \
  --plassembler-db "${PLASSEMBLER_DB}"
```

### Stage 4 example (GPU node)

```bash
#!/bin/bash
#SBATCH --job-name=prokaryont_s4
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=8:00:00
#SBATCH --partition=gpu

source ~/.bashrc
conda activate prokaryont-stage4-polish

# Dorado from standalone install
export PATH="/software/dorado-0.8.3-linux-x64/bin:$PATH"

bash 04_polish_orient.sh \
  --config pipeline.conf \
  --assembly autocycler_out/consensus_assembly.fasta \
  --pod5-dir pod5_dir/ \
  --device cuda:0 \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --double-polish
```

---

## 10. Changelog vs. original draft

This document was refined with the following corrections and hardening:

| Issue | Detail |
|---|---|
| **Missing `nanoporetech` channel** | `fastcat` is distributed on the `nanoporetech` conda channel, not `bioconda`. All environments using `fastcat` now include `-c nanoporetech`. |
| **Wrong SNIKT R dependencies** | The original listed `r-ggplot2`, `r-optparse`, `r-scales`, `r-gridextra`, `r-stringr`, `r-reshape2`. SNIKT actually requires `r-tidyverse`, `r-gridextra`, `r-docopt`, `r-lubridate`, plus the system tool `seqtk`. Using the bioconda `snikt` package resolves all of these automatically. |
| **`snikt.R` binary alias compatibility** | Documented and automated the `snikt.R -> snikt` symlink check in Stages 1, 2, and the setup script to satisfy `require_tool snikt.R`. |
| **Phantom `seqkit` in Stage 1** | `seqkit` is not used or checked anywhere in `01_qc_estimate.sh`. Removed from Stage 1. |
| **Phantom `edlib` in Stage 3** | `edlib` is statically linked into the Autocycler binary. It is not a separate runtime dependency. Removed. |
| **`nextdenovo` case sensitivity** | Documented that Bioconda installs `nextDenovo` (capital `D`) and added the `nextdenovo -> nextDenovo` symlink check to satisfy `03_autocycler_assemble.sh`. |
| **Canu Java runtime (`openjdk>=11`)** | Added explicit `openjdk>=11` to Stage 3 to prevent JVM memory/runtime failures on minimal HPC compute nodes. |
| **Missing `gnuplot`** | `mummerplot` requires `gnuplot` at runtime. Added to Stage 3. |
| **Missing `seqtk`** | SNIKT's backend uses `seqtk`. The bioconda `snikt` package pulls it automatically, but this was not documented. |
| **Missing `pyrodigal` in Stage 4** | Dnaapler uses `pyrodigal` for gene prediction. Added alongside `blast`. |
| **Missing `pplacer` in Stage 5** | GTDB-Tk requires `pplacer`. Added. |
| **iLesta availability** | Documented as a bioconda package, but it does not exist there. Added explicit warning about manual source installation. |
| **YAML pip section with broken syntax** | The Stage 1 YAML had a `pip:` section with `ont-fastcat` listed twice and invalid YAML indentation. Removed entirely; `fastcat` is installed via the `nanoporetech` channel. |
| **Broken git+https in Stage 3 pip** | The original listed `git+https://github.com/rrwick/Autocycler.git` under pip. Autocycler is a Rust binary available on bioconda; pip install from GitHub would fail. Removed. |
| **`ont-fastcat` package name** | The conda package name is `fastcat`, not `ont-fastcat` (that's a legacy PyPI/pip name). Corrected throughout. |
| **Stage 2 SNIKT as conditional** | Documented as "SNIKT re-assessment" suggesting it is optional, but `require_tool snikt.R` is unconditional in `02_preprocess_filter.sh`. Corrected to show it as a hard dependency. |
| **Stage 5 Python version** | GTDB-Tk + pplacer constrain Python ≤ 3.10. Added tip. |
| **No ground-truth dependency table** | Added §2 mapping every tool to the exact script check that requires it, with mandatory/conditional/optional classification. |
| **Missing assembler binary name mapping** | Added per-assembler table showing conda package name vs. checked binary name vs. additional deps (e.g., `metamdbg` package → `metaMDBG` binary). |
