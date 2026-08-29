# Conda Environment Setup for Each ProkaryONT Stage

These modular Conda/Mamba specifications are **audited against current script
checks** in the active ProkaryONT stages. They have not been claimed as solved,
installed, or exercised on representative data. Perform and record that
acceptance work on the target HPC system as described in
[HPC acceptance record](#hpc-acceptance-record).

## Environment strategy

Separate stage environments isolate the pipeline's C/C++, Rust, R, Java,
Python, and GPU-facing toolchains. Create the committed specifications one at a
time so that each solve and verification record remains attributable to one
stage.

### Channel configuration

The YAML files commit their channel order. If configuring Conda globally as
well, use strict channel priority:

~~~bash
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels nanoporetech
conda config --set channel_priority strict
~~~

<code>fastcat</code> is requested from the <code>nanoporetech</code> channel
in the Stage 1 and Stage 2 YAMLs. Use <code>mamba</code> or
<code>micromamba</code> when available, but record the exact solver and version
used during HPC acceptance.

## Dependency contracts

The following categories are deliberately separate. A package appearing in a
YAML does not by itself make it a direct script-level contract.

### Direct pipeline preflight checks

| Stage | Unconditional hard checks | Conditional hard checks | Soft or optional checks |
|---|---|---|---|
| 1 | <code>NanoPlot</code>, <code>lrge</code>, <code>raven</code>, <code>meryl</code>, <code>fastcat</code>, <code>porechop_abi</code>, <code>snikt.R</code> | - | <code>pigz</code> is used when present; otherwise <code>gzip</code> |
| 2 | <code>filtlong</code>, <code>seqkit</code>, <code>NanoPlot</code>, <code>snikt.R</code>, <code>fastplong</code> | <code>porechop_abi</code>, <code>chopper</code>, and <code>fastcat</code> when their modes are enabled | <code>pigz</code> is used when present; otherwise <code>gzip</code> |
| 3 | <code>autocycler</code>, <code>parallel</code>, <code>seqkit</code>, <code>minimap2</code>, <code>samtools</code>, <code>python3</code> | <code>rasusa</code> for Rasusa subsampling; <code>plassembler</code> unless the PLSDB screen is skipped | <code>nucmer</code> and <code>mummerplot</code> enable optional dotplots; <code>pigz</code> is optional |
| 4 | <code>dorado</code>, <code>samtools</code>, <code>dnaapler</code> | - | <code>pigz</code> is optional |
| 5 | <code>gtdbtk</code> | - | <code>mlst</code> and <code>barrnap</code> are skipped when absent |

### Stage 3 selected-backend checks

Stage 3 applies three distinct layers for every requested assembler:

1. <code>autocycler helper ASSEMBLER --help</code> is a hard compatibility
   gate. A selected helper task missing from the installed Autocycler release
   aborts the stage before assembly.
2. Backend executables are checked with <code>command -v</code>. Missing
   executables produce preflight warnings, not immediate hard failures.
3. A missing or broken backend then causes its GNU Parallel assembler job to
   fail when that job is executed.

The script default is the installable five-assembler ensemble
<code>flye,canu,hifiasm,miniasm,myloasm</code>. The complete canonical helper
catalog remains selectable explicitly:

| Helper task | Backend executable checks | Package/install note |
|---|---|---|
| <code>flye</code> | <code>flye</code> | <code>flye</code> |
| <code>canu</code> | <code>canu</code> | <code>canu</code>; Java runtime required |
| <code>hifiasm</code> | <code>hifiasm</code> | <code>hifiasm</code> |
| <code>miniasm</code> | <code>miniasm</code>, <code>minipolish</code>, <code>minimap2</code>, <code>racon</code> | corresponding Bioconda packages |
| <code>myloasm</code> | <code>myloasm</code> | <code>myloasm</code> |
| <code>ilesta</code> | <code>Ilesta</code>, <code>minipolish</code>, <code>minimap2</code>, <code>racon</code> | install iLesta manually |
| <code>lja</code> | <code>lja</code> | <code>lja</code> |
| <code>raven</code> | <code>raven</code> | <code>raven-assembler</code> |
| <code>metamdbg</code> | <code>metaMDBG</code> | <code>metamdbg</code> |
| <code>plassembler</code> | <code>plassembler</code> | also used by the PLSDB screen |
| <code>nextdenovo</code> | <code>nextdenovo</code> | package may expose <code>nextDenovo</code>; see Stage 3 note |
| <code>redbean</code> | <code>wtdbg2</code> | <code>wtdbg</code>; <code>wtdbg2</code> is accepted as an input alias and normalized to <code>redbean</code> |
| <code>necat</code> | <code>necat</code> | <code>necat</code> |

### Supporting package and runtime dependencies

The YAMLs also declare selected supporting dependencies that are not checked
directly by the stage scripts:

- Python runtimes and <code>pigz</code> for Stages 1-4.
- <code>openjdk</code> for Canu.
- <code>gnuplot</code> for <code>mummerplot</code>.
- <code>blast</code> and <code>pyrodigal</code> for Dnaapler.
- <code>hmmer</code>, <code>prodigal</code>, <code>fastani</code>,
  <code>mash</code>, and <code>pplacer</code> for the GTDB-Tk runtime.

### External databases and assets

Conda environments do not provide the raw reads, POD5 input, Dorado models,
Plassembler database, or GTDB-Tk reference data. Record those assets and their
versions separately from the environment solve.

## Stage 1 - QC and genome-size estimation

**Script:** [<code>01_qc_estimate.sh</code>](01_qc_estimate.sh)

**Environment:** [<code>envs/env_stage1_qc.yml</code>](envs/env_stage1_qc.yml)

~~~bash
mamba env create -f envs/env_stage1_qc.yml
conda activate prokaryont-stage1-qc
~~~

Stage 1 checks for <code>snikt.R</code>. If the installed recipe exposes only
<code>snikt</code>, create the compatibility symlink after inspecting both
paths:

~~~bash
[[ ! -f "${CONDA_PREFIX}/bin/snikt.R" && -f "${CONDA_PREFIX}/bin/snikt" ]] && \
  ln -s "${CONDA_PREFIX}/bin/snikt" "${CONDA_PREFIX}/bin/snikt.R"
~~~

Static command check:

~~~bash
bash 01_qc_estimate.sh -h
~~~

## Stage 2 - preprocessing and filtering

**Script:** [<code>02_preprocess_filter.sh</code>](02_preprocess_filter.sh)

**Environment:** [<code>envs/env_stage2_preproc.yml</code>](envs/env_stage2_preproc.yml)

~~~bash
mamba env create -f envs/env_stage2_preproc.yml
conda activate prokaryont-stage2-preproc
~~~

Apply the same inspected <code>snikt.R</code> compatibility symlink if needed.
The YAML includes all three conditional tools so the corresponding opt-in modes
can be enabled without changing environments.

Static command check:

~~~bash
bash 02_preprocess_filter.sh -h
~~~

## Stage 3 - multi-assembler Autocycler workflow

**Script:** [<code>03_autocycler_assemble.sh</code>](03_autocycler_assemble.sh)

**Environment:** [<code>envs/env_stage3_assemble.yml</code>](envs/env_stage3_assemble.yml)

~~~bash
mamba env create -f envs/env_stage3_assemble.yml
conda activate prokaryont-stage3-assemble
~~~

The committed YAML installs the core workflow, conditional Rasusa and
Plassembler modes, optional MUMmer dotplots, and the default assembler set:

~~~text
flye,canu,hifiasm,miniasm,myloasm
~~~

Install only the extra backends needed for a non-default
<code>--assemblers</code> selection. For example:

~~~bash
mamba install -n prokaryont-stage3-assemble \
  -c conda-forge -c bioconda \
  raven-assembler metamdbg lja wtdbg necat nextdenovo
~~~

iLesta is not included in the YAML; install it separately and confirm that
<code>Ilesta</code>, <code>minipolish</code>, <code>minimap2</code>, and
<code>racon</code> are all on <code>PATH</code>.

If a NextDenovo installation exposes only <code>nextDenovo</code>, create a
lowercase compatibility symlink after inspecting both paths:

~~~bash
[[ ! -f "${CONDA_PREFIX}/bin/nextdenovo" && -f "${CONDA_PREFIX}/bin/nextDenovo" ]] && \
  ln -s "${CONDA_PREFIX}/bin/nextDenovo" "${CONDA_PREFIX}/bin/nextdenovo"
~~~

Download the Plassembler database separately:

~~~bash
plassembler download -d /databases/plassembler_db
export PLASSEMBLER_DB=/databases/plassembler_db
~~~

Stage 3 resolves the database in this order: <code>--plassembler-db</code>,
<code>plassembler_db</code> in <code>pipeline.conf</code>,
<code>PLASSEMBLER_DB</code>, then
<code>$CONDA_PREFIX/plassembler_db</code>.

Static command check:

~~~bash
bash 03_autocycler_assemble.sh -h
~~~

## Stage 4 - Dorado polishing and Dnaapler orientation

**Script:** [<code>04_polish_orient.sh</code>](04_polish_orient.sh)

**Environment:** [<code>envs/env_stage4_polish.yml</code>](envs/env_stage4_polish.yml)

~~~bash
mamba env create -f envs/env_stage4_polish.yml
conda activate prokaryont-stage4-polish
~~~

Dorado is a direct hard requirement but is intentionally absent from the YAML.
Install an official standalone Dorado build that matches the target HPC
architecture and GPU driver, add it to <code>PATH</code>, and install the
required basecalling model.

For a BAM with multiple read groups, use the non-interactive
<code>--read-group ID</code> option. The script validates the exact ID against
the aligned BAM header and passes <code>--RG ID</code> to both polish passes:

~~~bash
bash 04_polish_orient.sh \
  --assembly autocycler_out/consensus_assembly.fasta \
  --pod5-dir pod5_dir/ \
  --read-group sample_rg \
  --double-polish
~~~

Static command check:

~~~bash
bash 04_polish_orient.sh -h
~~~

## Stage 5 - taxonomy classification

**Script:** [<code>05_taxonomy.sh</code>](05_taxonomy.sh)

**Environment:** [<code>envs/env_stage5_taxonomy.yml</code>](envs/env_stage5_taxonomy.yml)

~~~bash
mamba env create -f envs/env_stage5_taxonomy.yml
conda activate prokaryont-stage5-taxonomy
~~~

The YAML lets GTDB-Tk constrain its compatible Python runtime. Download the
GTDB-Tk reference database separately and record its release:

~~~bash
download-db.sh /databases/gtdbtk
export GTDBTK_DATA_PATH=/databases/gtdbtk
~~~

Static command check:

~~~bash
bash 05_taxonomy.sh -h
~~~

## Sequential environment creation

Run these commands one by one and stop to review each solve:

~~~bash
mamba env create -f envs/env_stage1_qc.yml
mamba env create -f envs/env_stage2_preproc.yml
mamba env create -f envs/env_stage3_assemble.yml
mamba env create -f envs/env_stage4_polish.yml
mamba env create -f envs/env_stage5_taxonomy.yml
~~~

Do not describe the environments as tested merely because the YAML parses or
solves. A representative-data run is a separate validation level.

## HPC acceptance record

For every completed solve/install validation, record:

1. OS, architecture, and exact Conda/Mamba version.
2. Solve date, strict-priority setting, and effective channel order.
3. An exact environment export or lockfile committed or archived with the run.
4. Tool versions actually verified.
5. Whether each stage was only checked with <code>-h/--help</code> or exercised
   on representative data.

A compact record can use this structure:

~~~text
Stage:
OS/architecture:
Conda/Mamba version:
Solve date:
Channel order and priority:
Export or lockfile:
Verified tool versions:
Validation level: --help only | representative data
Dataset/run reference:
~~~

## Slurm invocation examples

Stage 3 with the committed default ensemble:

~~~bash
#!/bin/bash
#SBATCH --job-name=prokaryont_s3
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=24:00:00

source ~/.bashrc
conda activate prokaryont-stage3-assemble
export PLASSEMBLER_DB=/databases/plassembler_db

bash 03_autocycler_assemble.sh \
  --config pipeline.conf \
  --reads filtered_input.fastq.gz \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --skip-curation \
  --plassembler-db "${PLASSEMBLER_DB}"
~~~

Stage 4 with a reproducible read-group selection:

~~~bash
#!/bin/bash
#SBATCH --job-name=prokaryont_s4
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --gres=gpu:1
#SBATCH --time=8:00:00

source ~/.bashrc
conda activate prokaryont-stage4-polish
export PATH="/software/dorado/bin:$PATH"

bash 04_polish_orient.sh \
  --config pipeline.conf \
  --assembly autocycler_out/consensus_assembly.fasta \
  --pod5-dir pod5_dir/ \
  --device cuda:0 \
  --read-group sample_rg \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --double-polish
~~~
