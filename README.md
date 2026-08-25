# ProkaryONT

ProkaryONT is a step-by-step Oxford Nanopore bacterial genome assembly
workflow. Each active stage is invoked separately so its outputs can be
reviewed before the next stage starts. There is no active master orchestrator.

## Active workflow

| Stage | Script | Purpose |
|---|---|---|
| 1 | `01_qc_estimate.sh` | Read QC and genome-size estimation |
| 2 | `02_preprocess_filter.sh` | Optional trimming plus quality/length filtering |
| 3 | `03_autocycler_assemble.sh` | Subsampling, multi-assembler Autocycler resolution, combine depth annotation, and PLSDB characterization |
| 4 | `04_polish_orient.sh` | Dorado polishing and Dnaapler reorientation |
| 5 | `05_taxonomy.sh` | GTDB-Tk classification, MLST, and 16S extraction |

The former `prokaryont.sh`, `run_all.sh`, Stage 6 assessment, and Stage 7
annotation stack are retained under `archived/` for historical reference. They
are not active entrypoints.

## Core requirements

- QC/filtering: `NanoPlot`, `fastcat`, `porechop_abi`, `snikt.R`, `seqkit`,
  `filtlong`, `fastplong`, and optional `chopper`
- Assembly: current `autocycler`, GNU `parallel`, `samtools`, `minimap2`, and
  the selected assembler helpers
- Optional Rasusa mode: a Rasusa release exposing `rasusa reads`
- PLSDB screen: `plassembler` plus its database
- Polishing/orientation: `dorado`, `samtools`, and `dnaapler`
- Taxonomy: `gtdbtk`; `mlst` and `barrnap` are optional

The canonical Stage 3 helper set is:

```text
flye,canu,hifiasm,ilesta,lja,raven,miniasm,metamdbg,myloasm,plassembler,nextdenovo,redbean,necat
```

`wtdbg2` remains accepted as a compatibility alias and is normalized to the
current Autocycler helper task `redbean`. iLesta additionally requires
`Ilesta`, `minipolish`, `minimap2`, and `racon`; LJA requires `lja`.

## Environment setup

For comprehensive per-stage Conda/Mamba YAML environment files, Slurm batch
templates, and package verification steps, see:

👉 [`conda-environment-setup-for-each-prokaryont-script.md`](conda-environment-setup-for-each-prokaryont-script.md)

### Essential Conda channel configuration

```bash
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels nanoporetech   # required for fastcat
conda config --set channel_priority strict
```

> [!TIP]
> Use modular per-stage environments (`prokaryont-stage1-qc` through
> `prokaryont-stage5-taxonomy`) to prevent solver conflicts across C++, R,
> Java, and CUDA runtimes. Dorado should be installed via its official
> standalone release rather than Conda to ensure host GPU driver compatibility.

## Step-by-step execution

Copy and edit `pipeline.conf`, then run each stage explicitly. CLI values take
precedence over config values.

### 1. QC and genome-size estimation

```bash
bash 01_qc_estimate.sh \
  --config pipeline.conf \
  --input-fastq input.fastq.gz \
  --sequencing-summary sequencing_summary.txt
```

Review `01_qc/` and `02_genome_size/` before preprocessing.

### 2. Preprocess and filter

```bash
bash 02_preprocess_filter.sh \
  --config pipeline.conf \
  --input-fastq input.fastq.gz
```

Filtlong `--keep_percent` is opt-in. Set `filtlong_keep_percent=N` or pass
`--keep-percent N`; otherwise only the configured minimum-length filter is
applied by Filtlong.

### 3. Assemble with Autocycler

Select assemblers directly on the Stage 3 CLI:

```bash
bash 03_autocycler_assemble.sh \
  --config pipeline.conf \
  --reads filtered_input.fastq.gz \
  --read-type ont_r10 \
  --assemblers flye,canu,hifiasm,ilesta,lja,raven,miniasm,metamdbg,myloasm,plassembler,nextdenovo,redbean,necat \
  --plassembler-db /path/to/plassembler_db
```

The default subsampler remains native Autocycler. Rasusa can generate the same
`subsampled_reads/sample_NN.fastq` contract with deterministic distinct seeds:

```bash
bash 03_autocycler_assemble.sh \
  --reads filtered_input.fastq.gz \
  --read-type ont_r10 \
  --genome-size 5000000 \
  --assemblers flye,lja,ilesta \
  --subsampler rasusa \
  --subsample-count 4 \
  --subsample-seed 13 \
  --rasusa-coverage 60 \
  --plassembler-db /path/to/plassembler_db
```

Without `--rasusa-coverage`, Stage 3 derives an ensemble-oriented target from
input depth and Autocycler's minimum-depth formula. Independently sampled
Rasusa subsets can overlap substantially at high requested coverage; native
Autocycler subsampling remains preferable when maximal subset independence is
important.

`autocycler compress` and `autocycler combine` are each capped at 100 threads.
Override their independent caps with `--compress-threads N` and
`--combine-threads N` (`1..100`). Combine receives the full Stage 3 input reads
by default for final depth annotation; use `--combine-reads PATH` only for an
intentional alternative. The depth k-mer defaults to 19 and can be changed with
`--combine-depth-kmer` to an odd value from 11 through 31.

Immediately after combine, Stage 3 runs a non-mutating PLSDB similarity screen
against every Autocycler consensus contig and publishes
`plassembler_summary.tsv`. This is characterization evidence, not a plasmid
classification or replacement assembly. `autocycler_out/consensus_assembly.fasta`
remains authoritative. Use `--skip-plsdb-screen` only when the screen is
intentionally omitted.

### 4. Polish and orient

```bash
bash 04_polish_orient.sh \
  --config pipeline.conf \
  --assembly autocycler_out/consensus_assembly.fasta \
  --pod5-dir pod5_dir/ \
  --device cuda:0
```

`--device` is passed only to `dorado polish`. It is reused for the optional
second polish pass and does not constrain Dorado basecalling or alignment.

### 5. Taxonomy

```bash
bash 05_taxonomy.sh \
  --config pipeline.conf \
  --assembly dnaapler_reoriented.fasta \
  --gtdbtk-db /path/to/gtdbtk_db
```

## Important outputs

```text
01_qc/                                      Stage 1 and post-filter QC
02_genome_size/                             Genome-size estimates
assemblies/                                 Standardized helper assemblies and job logs
autocycler_out/consensus_assembly.fasta     Authoritative Stage 3 consensus
autocycler_out/consensus_assembly.gfa       Combined graph with depth metadata
autocycler_out/consensus_assembly.yaml      Combined metadata
autocycler_out/plassembler_plsdb/           PLSDB screen provenance
plassembler_summary.tsv                     One PLSDB result row per consensus contig
rasusa_subsample.tsv                        Rasusa subset provenance, when selected
metrics.tsv                                 Autocycler metrics
contig_depths.tsv                           Independent mapping-based depth report
polished_assembly.fasta                     Dorado-polished assembly
dnaapler_reoriented.fasta                   Reoriented final assembly
07_taxonomy/                                Stage 5 outputs
```

## Laptop-safe mock checks

The mock environment supports syntax and generated-command checks without
claiming real bioinformatics validation:

```bash
python utils/generate_mock_environment.py
export PATH="$(pwd)/mock_bin:$PATH"
export PLASSEMBLER_DB="$(pwd)/plassembler_db"
# run stage commands with --dry-run
python utils/generate_mock_environment.py --clean
```

Real assembler, database, and GPU validation should be performed on the target
HPC environment with representative reads and the production databases.

## License

Apache License 2.0; see `LICENSE`.
