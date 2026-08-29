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
- Assembly: current `autocycler`, GNU `parallel`, `seqkit`, `samtools`,
  `minimap2`, and the selected assembler helpers
- Optional Rasusa mode: a Rasusa release exposing `rasusa reads`
- PLSDB screen: `plassembler` plus its database
- Polishing/orientation: `dorado`, `samtools`, and `dnaapler`
- Taxonomy: `gtdbtk`; `mlst` and `barrnap` are optional

Stage 3 defaults to this practical five-assembler set:

```text
flye,canu,hifiasm,miniasm,myloasm
```

The complete canonical helper catalog is:

```text
flye,canu,hifiasm,ilesta,lja,raven,miniasm,metamdbg,myloasm,plassembler,nextdenovo,redbean,necat
```

`wtdbg2` remains accepted as a compatibility alias and is normalized to the
current Autocycler helper task `redbean`. iLesta additionally requires
`Ilesta`, `minipolish`, `minimap2`, and `racon`; LJA requires `lja`.

For each selected assembler, Stage 3 first applies a hard
`autocycler helper ASSEMBLER --help` compatibility gate. Backend executable
discovery is then a soft warning; a missing backend will still cause its GNU
Parallel job to fail when executed.

## Environment setup

The committed per-stage Conda/Mamba files are:

- [`env_stage1_qc.yml`](envs/env_stage1_qc.yml)
- [`env_stage2_preproc.yml`](envs/env_stage2_preproc.yml)
- [`env_stage3_assemble.yml`](envs/env_stage3_assemble.yml)
- [`env_stage4_polish.yml`](envs/env_stage4_polish.yml)
- [`env_stage5_taxonomy.yml`](envs/env_stage5_taxonomy.yml)

For creation commands, dependency categories, Slurm examples, and the HPC
validation record, see:

👉 [`conda-environment-setup-for-each-prokaryont-script.md`](conda-environment-setup-for-each-prokaryont-script.md)

The specifications are audited against current script checks; solve/install
and representative-data validation remain HPC acceptance tasks.

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

Common CLI/config mappings:

| CLI option | `pipeline.conf` key |
|---|---|
| `--output-dir` | `output_dir` |
| `--output-fastq` | `output_fastq` |
| `--sample-name` | `sample_name` |
| `--tmp-dir` | `tmp_dir` |
| `--min-length` | `filtlong_min_length` |
| `--keep-percent` | `filtlong_keep_percent` |
| `--genome-size` | `genome_size_override` |
| `--assemblers` | `assemblers` |
| `--combine-reads` | `autocycler_combine_reads` |
| `--combine-depth-kmer` | `autocycler_combine_depth_kmer` |
| `--combine-threads` | `autocycler_combine_threads` |
| `--compress-threads` | `autocycler_compress_threads` |
| `--device` | `dorado_polish_device` |
| `--read-group` | `dorado_read_group` |
| `--gtdbtk-db` | `gtdbtk_data_path` |

`PROKARYONT_ASSEMBLERS` is an environment override; `assemblers=` is the
lowercase config key. Pass assembler-specific options through scalar
`*_extra_args=` config values. Config values are trimmed strings rather than
shell syntax, so complex nested quoting is not preserved.

### 1. QC and genome-size estimation

```bash
bash 01_qc_estimate.sh \
  --config pipeline.conf \
  --input-fastq input.fastq.gz
```

`--sequencing-summary` is optional; when supplied, it adds the summary-based
NanoPlot report. All other Stage 1 diagnostics remain default-on and expose
individual `--skip-*` flags. Stage 1 records LRGE and Raven independently in
`02_genome_size/genome_size_estimates.tsv`; it does not create a weighted
authoritative estimate. Review `01_qc/` and `02_genome_size/` before
preprocessing.

### 2. Preprocess and filter

```bash
bash 02_preprocess_filter.sh \
  --config pipeline.conf \
  --input-fastq input.fastq.gz
```

Filtlong `--keep_percent` is opt-in. Set `filtlong_keep_percent=N` or pass
`--keep-percent N`; otherwise SeqKit applies the probability-space mean-Q and
minimum-length filters in one pass and Filtlong is not invoked. Stage 2 writes
`02_filtering_metrics.tsv` with read/base changes, retention, length, N50,
quality, and coverage (when `--genome-size` is explicitly supplied). Zero-read
output is a hard failure; no universal 65% retention threshold is imposed. Use
`--skip-intermediate-metrics` when scratch/runtime matters more than per-step
attribution; input and final metrics are still retained.

Stages 1 and 2 preserve the historical current-directory layout by default.
Use `--output-dir`, `--output-fastq`, `--sample-name`, and `--tmp-dir` for an
explicit namespace. Each stage writes a completion manifest last and resumes
only when its input metadata, effective parameters, tool versions, and expected
outputs match. A mismatch fails closed unless `--force` is supplied. Metadata
identity is the default; `--input-sha256` opts into a full input hash.

`--dry-run` is mutation-free: it validates the CLI, inputs, enabled tool
contracts, and prints the commands, but creates no directories, temporary
files, logs, outputs, or manifests.

### 3. Assemble with Autocycler

This minimal invocation uses the default five-assembler ensemble, avoids an
interactive pause, and explicitly omits the optional database screen. It is a
first-run example, not a universally optimal biological ensemble:

```bash
bash 03_autocycler_assemble.sh \
  --config pipeline.conf \
  --reads filtered_input.fastq.gz \
  --read-type ont_r10 \
  --threads 48 \
  --parallel-jobs 4 \
  --canu-parallel-jobs 2 \
  --skip-curation \
  --skip-plsdb-screen
```

Omitting `--assemblers` selects
`flye,canu,hifiasm,miniasm,myloasm`. Pass `--assemblers` explicitly to use a
different subset of the canonical catalog, after installing both the matching
Autocycler helper tasks and backend executables.

The default subsampler remains native Autocycler. Rasusa can generate the same
`subsampled_reads/sample_NN.fastq` contract with deterministic distinct seeds:

```bash
bash 03_autocycler_assemble.sh \
  --reads filtered_input.fastq.gz \
  --read-type ont_r10 \
  --genome-size 5000000 \
  --assemblers flye,miniasm,myloasm \
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

After assembler jobs finish, Stage 3 preserves `assemblies/` and constructs the
exact Autocycler inputs under `assemblies_prepared/`. Every non-empty assembly
is validated and counted. Only an individual assembly with more than
`dedup_trigger_contigs` records (default 25) is self-aligned and assessed for
contained contigs. Automatic containment requires at least 99% single-target
coverage and 99% aggregate identity; equal-length records are collapsed only
when their complete sequences are identical, including reverse-complement
identity. No contig is removed solely for being short by default.

The pre/post counts, decisions, per-contig event TSVs, tool logs, versions,
input-policy fingerprint, and prepared-output SHA-256 manifest are written to
`assembly_assessment/`. Successful PAF
files are deleted after the event reports reconcile; failure PAFs remain in the
temporary failure directory for diagnosis. Use `--skip-contained-dedup` for an
assessment-only run.

Autocycler's `--max_contigs` is a cohort-wide mean guard, not a per-file
maximum. Stage 3 passes an explicit value to both `compress` and `cluster`. It
uses 25 when the prepared mean is at most 25 and otherwise stops before
compression with the minimum required override. An intentional override can be
provided with `--max-contigs N`; it is never raised automatically.

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
For a non-interactive read-group choice, pass `--read-group ID` or set
`dorado_read_group=ID`. The script rejects control or non-ASCII characters,
verifies the exact ID in each aligned BAM header, and passes `--RG ID` to
both Dorado polish passes. With multiple groups and no explicit ID,
`--skip-curation` retains the existing `--ignore-read-groups` behavior.

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
02_genome_size/genome_size_estimates.tsv   Independent LRGE/Raven diagnostics
02_filtering_metrics.tsv                    Stage-by-stage read/base/Q/length metrics
filtered_input.fastq.gz                     Authoritative Stage 2 to Stage 3 handoff
.prokaryont_stage1.complete.tsv             Stage 1 identity/completion contract
.prokaryont_stage2.complete.tsv             Stage 2 identity/completion contract
provenance/                                  Run-specific commands, options, tools, and versions
subsample.yaml                              Native Autocycler subsampling provenance, when present
subsampled_reads/                           Subset metadata directory; sample FASTQs are removed after assembly
assemblies/                                 Standardized helper assemblies and job logs
assemblies_prepared/                        Validated, weighted inputs supplied to Autocycler
assembly_assessment/                        Pre/post counts, event TSVs, logs, and fingerprints
autocycler_out/consensus_assembly.fasta     Authoritative Stage 3 consensus
autocycler_out/consensus_assembly.gfa       Combined graph with depth metadata
autocycler_out/consensus_assembly.yaml      Combined metadata
autocycler_out/plassembler_plsdb/           PLSDB screen provenance
plassembler_summary.tsv                     One PLSDB result row per consensus contig
rasusa_subsample.tsv                        Rasusa subset provenance, when selected
metrics.tsv                                 Autocycler metrics
contig_depths.tsv                           Independent mapping-based depth report
all_reads_w_moves.bam                       Dorado move-table basecalls; conditional cleanup
aligned.sorted.bam                          First-pass alignment; conditional cleanup
aligned_reoriented.sorted.bam               Optional second-pass alignment; conditional cleanup
polished_assembly.fasta                     Dorado-polished assembly
dnaapler_reoriented.fasta                   Reoriented final assembly
07_taxonomy/gtdbtk.bac120.summary.tsv        Principal bacterial GTDB-Tk result, when produced
07_taxonomy/gtdbtk.ar53.summary.tsv          Principal archaeal GTDB-Tk result, when produced
07_taxonomy/mlst_result.tsv                  Optional MLST result
07_taxonomy/rrna_predictions.gff3            Optional Barrnap rRNA calls
07_taxonomy/16S_sequences.fasta              Optional Barrnap sequence export
07_taxonomy/16S_only.fasta                   Filtered 16S-only sequence export
```

Without `--cleanup-bam`, the Stage 4 BAMs are retained. With cleanup enabled,
`aligned.sorted.bam` is removed after the first pass; `all_reads_w_moves.bam`
is retained only until a requested second pass completes, and the second-pass
alignment is then removed.

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
