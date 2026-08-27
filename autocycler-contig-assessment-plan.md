# Autocycler contig assessment, conditional containment cleanup, and `--max_contigs` plan

Status: implemented in the working tree; local validation passed, HPC acceptance pending
Repository baseline: `main` at `e674d07`
Real-data validation target: the user's HPC environment, not the laptop

## 1. Goal

Add a deterministic preparation step between parallel assembly and Autocycler
compression that:

1. counts FASTA records in every surviving assembly;
2. identifies only highly fragmented individual assemblies;
3. runs `dedup_contained.py` only on those flagged assemblies;
4. preserves the original assembler outputs;
5. recounts the prepared assembly cohort;
6. passes one explicitly recorded `--max_contigs` value to both
   `autocycler compress` and `autocycler cluster`; and
7. stops with actionable diagnostics instead of silently weakening Autocycler's
   fragmentation guard.

This is an input-preparation and safety feature. It must not become a general
assembly-cleaning step or silently remove plausible replicons.

## 2. Repository and upstream findings

### Baseline Stage 3 flow at `e674d07`

Before this plan was implemented, `03_autocycler_assemble.sh`:

- finishes GNU Parallel assembly jobs in step 9d;
- immediately edits FASTA headers in place to add Autocycler weights in 9e;
- reports empty files and total assembly size at the first curation pause;
- runs `autocycler compress` and `autocycler cluster` without
  `--max_contigs`, so both use their upstream default; and
- supplies `assemblies/` directly to `autocycler compress`.

The new step belongs after 9d and before header weighting. Weighting should be
applied to prepared copies, not to pristine assembler FASTAs.

### What `--max_contigs` actually means

Current Autocycler exposes `--max_contigs` on both `compress` and `cluster`,
with a default of 25 and the description "refuse to run if mean contigs per
assembly exceeds this value." The cluster implementation calculates:

```text
mean_contigs = total input sequences / number of input assemblies
```

and refuses only when `mean_contigs > max_contigs`.

Consequences:

- it is a cohort-wide mean, not the maximum contig count of any one FASTA;
- one fragmented assembly does not necessarily trip the guard;
- the value does not tune clustering or improve fragmented input; it is a
  circuit breaker intended to prevent obviously bad input from making
  clustering unreasonably expensive; and
- automatically setting it to `ceil(observed_mean)` would always let the input
  through and would largely defeat the guard.

Therefore ProkaryONT should compute and expose the evidence needed for a
decision, but should not silently auto-raise the value above 25.

### Pre-implementation `dedup_contained.py` audit

The original workspace script was a useful starting point. It read a FASTA and a
self-alignment PAF, removes short records, detects query coverage by qualifying
alignments to equal-or-longer targets, and writes retained records to stdout.

Before automatic pipeline use, the following behaviors need correction or
explicit policy:

1. **Cross-target pooling:** query intervals are currently pooled across every
   target. A query covered 50% by target A and 50% by target B can therefore be
   deleted even though it is not contained in either target. Coverage must be
   grouped by `(query, target)` and at least one single target must satisfy the
   containment rule.
2. **Aggressive defaults:** 90% coverage, 90% identity, and unconditional
   removal below 500 bp are too permissive for an automatic biological-data
   mutation. Automatic mode should default to near-exact containment
   (`min_cov=0.99`, `min_identity=0.99`) and no length-only deletion
   (`min_len=0`). The existing 0.90/0.90/500 behavior can remain available as
   an explicit aggressive configuration.
3. **Small-replicon safety:** a sequence must not be removed solely because it
   is short unless the user explicitly configures a positive minimum length.
4. **Input validation:** duplicate FASTA IDs, empty sequences, malformed FASTA,
   malformed PAF rows, unknown PAF IDs, inconsistent PAF lengths, and invalid
   thresholds must fail clearly.
5. **Auditability:** the script currently reports only totals. It must emit a
   per-record TSV containing the action, reason, matched target, coverage,
   identity, and length.
6. **Atomic use:** stdout must be captured to a temporary file, validated, and
   renamed only after success. A failed run must never leave a partial prepared
   FASTA.

## 3. Proposed policy

### 3.1 Always assess; deduplicate conditionally

Count records in all non-empty `assemblies/*.fasta` files. Empty or missing
outputs are not part of the mean and must be reported separately. A non-empty
malformed FASTA is treated as corruption and fails the preparation stage.

An individual assembly is a deduplication candidate only when:

```text
raw_contig_count > dedup_trigger_contigs
```

Use 25 as the initial default for `dedup_trigger_contigs`, because it aligns
with Autocycler's documented default guard while applying it per assembly as a
pipeline heuristic. Make the trigger configurable and label it clearly: it is
not an upstream per-assembly rule and does not prove that the assembly is bad.

Assemblies at or below the trigger must not invoke minimap2 or
`dedup_contained.py`; they are copied byte-for-byte into the prepared input
directory. This satisfies the requirement that containment cleanup run only
for fragmented assemblies.

### 3.2 Preserve source assemblies

Use separate runtime locations:

```text
assemblies/                    pristine assembler outputs and job logs
assemblies_prepared/           exact inputs supplied to Autocycler
assembly_assessment/
  assemblies.tsv              pre/post counts and decisions
  max_contigs.tsv              cohort calculation and effective value
  dedup/<assembly>.paf         temporary self-alignment evidence; retained only on failure
  dedup/<assembly>.events.tsv  per-contig keep/drop evidence
  input_manifest.sha256        source and policy fingerprint
  prepared_manifest.sha256     exact prepared-output fingerprint
```

Apply Plassembler cluster weights and Canu/Flye/Hifiasm consensus weights only
to `assemblies_prepared/*.fasta`. Do not edit `assemblies/*.fasta` in place.

### 3.3 Self-alignment command

For each candidate assembly, use minimap2's assembly preset and self-homology
behavior, with the pipeline's allocated deduplication thread count:

```bash
minimap2 -x asm5 -DP -t "$dedup_threads" assembly.fasta assembly.fasta
```

`-D` suppresses trivial perfect diagonal matches and `-P` retains all chains,
which avoids relying on the normal default limit of five secondary mappings.
PAF columns 10 and 11 supply matching and alignment-block lengths, so CIGAR
generation is not required for the current containment calculation.

Pin the exact command and minimap2 version in the assessment report. If HPC
validation shows that `-P` creates excessive PAFs for highly repetitive input,
change the mapping policy only with a regression fixture demonstrating that
contained duplicates are still found.

### 3.4 Containment decision

For each query:

1. ignore the diagonal self-hit;
2. consider a strictly longer target; collapse equal-length records only when
   their complete sequences are identical (including reverse-complement
   identity), using a stable lexical tie-break so exactly one duplicate
   survives;
3. group qualifying alignment intervals by `(query, target, strand)`;
4. merge query intervals only within that group;
5. calculate union query coverage and aggregate identity for the group; and
6. drop the query only if one group meets both thresholds.

The event report must name the winning target and all decision metrics. A query
whose coverage is distributed across multiple targets remains retained.

### 3.5 `--max_contigs` decision

After preparation, calculate:

```text
post_mean = post_total_contigs / prepared_assembly_count
required_integer = ceil(post_mean)
```

Recommended default behavior:

- If `post_mean <= 25`, set `effective_max_contigs=25` and pass it explicitly
  to both commands.
- If `post_mean > 25`, stop before compression, retain all reports, and print
  the minimum explicit override (`required_integer`). Do not silently raise the
  guard.
- If the user supplies `--max-contigs N`, require `N >= required_integer` and
  pass the same `N` to both commands. Reject zero, negative, non-integer, or
  insufficient values.

Commands become:

```bash
autocycler compress \
  -i "$prepared_assemblies_dir" \
  -a "$autocycler_dir" \
  --threads "$compress_threads" \
  --max_contigs "$effective_max_contigs"

autocycler cluster \
  -a "$autocycler_dir" \
  --max_contigs "$effective_max_contigs"
```

This lets the script decide and record whether the default is valid, while an
intentional guard increase remains a visible user decision. A truly adaptive
raise can be reconsidered only after HPC measurements establish a defensible
resource ceiling; contig count alone cannot predict safe clustering RAM/time.

## 4. CLI and configuration contract

Add Stage 3 CLI options:

```text
--max-contigs N              Explicit Autocycler compress/cluster guard.
                             Default: guarded automatic assessment using 25.
--dedup-trigger-contigs N    Run containment cleanup only above this per-file
                             record count. Default: 25.
--dedup-min-cov FLOAT        Single-target query coverage. Default: 0.99.
--dedup-min-identity FLOAT   Alignment identity. Default: 0.99.
--dedup-min-len N            Explicit length-only filter; 0 disables it.
                             Default: 0.
--dedup-threads N            Threads for one minimap2 self-alignment.
                             Default: min(global threads, 16).
--skip-contained-dedup       Count and decide max_contigs without mutating any
                             assembly; useful for diagnosis and comparison.
```

Configuration keys:

```bash
# autocycler_max_contigs=       # unset means guarded automatic assessment
dedup_trigger_contigs=25
dedup_min_cov=0.99
dedup_min_identity=0.99
dedup_min_len=0
# dedup_threads=16
skip_contained_dedup=false
```

CLI must override config. Validate numeric values before subsampling or
assembly begins. `min_cov` and `min_identity` must be in `(0, 1]`; all count and
thread values must be positive integers except `dedup_min_len`, where zero is
valid.

Do not name the option `--auto-max-contigs`: automatic mode does not raise the
guard; it assesses whether the upstream default can be used.

## 5. Stage 3 integration

### Step ordering

Refactor the current post-assembly sequence to:

1. **9e - Assess and prepare assemblies**
   - enumerate only regular, non-empty `.fasta` files;
   - count records and bases;
   - conditionally self-align/deduplicate candidates;
   - copy non-candidates unchanged;
   - validate every prepared FASTA;
   - calculate raw and prepared cohort means;
   - choose or reject `effective_max_contigs`;
   - write reports and manifest.
2. **9f - Apply Autocycler weights** to prepared copies.
3. **9g - Cleanup subsampled reads.**
4. **9h - Curation point 1.** Show per-assembly raw/post counts, IDs removed,
   total bases removed, raw/post cohort means, and the effective max-contigs
   decision. Keep the existing failed-job and empty-output summary.
5. **9i - Compress and cluster** using `assemblies_prepared/` and the same
   explicit guard on both commands.
6. Renumber the remaining internal 9x labels and comments consistently.

The curation message must distinguish three states:

- not fragmented, cleanup not run;
- fragmented, cleanup ran and changed/no-change result; and
- fragmented, cleanup skipped by user.

### Transaction and failure rules

- Build `assemblies_prepared.tmp.<pid>/`, validate it, then rename it into
  place.
- Never overwrite or delete source FASTAs in `assemblies/`.
- Fail if no valid assembly FASTAs remain.
- Fail if deduplication removes every record from a candidate assembly.
- Fail if post-dedup record count or total bases unexpectedly increases.
- Fail if the event report and pre/post record counts disagree.
- A candidate with no contained sequences is valid and is copied through with
  `changed=false`.
- Keep PAF and event evidence on failure. After a successful preparation and
  report reconciliation, delete PAFs while retaining event TSVs and tool logs.
- Dry-run prints which files would be assessed, which gate would be applied,
  and the final command shape, but must clearly label the max-contigs value as
  unavailable when no real FASTAs were counted.

### Resume behavior

Fingerprint source FASTA SHA-256 values plus all preparation settings and the
`dedup_contained.py` SHA-256. Reuse `assemblies_prepared/` only when the
fingerprint matches exactly and every prepared FASTA matches its recorded
output SHA-256; otherwise rebuild it.

- Resume from `cluster`: assessment/preparation may be reused or rebuilt, then
  compress/cluster reruns with its recorded guard.
- Resume from `trim` or later: do not change prepared inputs underneath the
  existing Autocycler graph. Verify and report the saved manifest/decision;
  treat missing legacy metadata as a warning rather than retroactively
  deduplicating the run.
- A changed max-contigs value invalidates compress/cluster reuse but does not
  require rerunning assemblers.

## 6. `dedup_contained.py` implementation work

Keep the script standalone and dependency-free (Python standard library only),
but make these changes before pipeline integration:

1. add typed argument validators and actionable exit messages;
2. parse FASTA strictly and reject duplicate first-token IDs;
3. validate all PAF mandatory fields, IDs, coordinates, and FASTA lengths;
4. group intervals by one target and strand;
5. use a deterministic equal-length duplicate tie-break;
6. separate `contained` and `short` reasons, with the short rule disabled at
   zero;
7. add `--report-tsv PATH` and include one row per input contig;
8. include selected target, covered bp, coverage, identity, and keep/drop
   reason;
9. retain stdout FASTA behavior for Unix composability; and
10. return non-zero on any validation or write error.

Suggested event TSV fields:

```text
query_id query_length action reason target_id target_length strand
covered_bp coverage identity min_cov min_identity min_len
```

## 7. Exact files in implementation scope

- `dedup_contained.py` - harden containment logic and reporting.
- `03_autocycler_assemble.sh` - assessment/preparation gate, CLI/config
  plumbing, reports, max-contigs decision, prepared input directory, resume
  behavior, and curation summary.
- `pipeline.conf` - document defaults and explicit override.
- `README.md` - explain the new post-assembly preparation step and examples.
- `testing_guidelines.md` - laptop mock checks and HPC validation boundary.
- `utils/generate_mock_environment.py` - realistic FASTA/PAF fixtures and mock
  minimap2 behavior.
- `.gitignore` - ignore `assemblies_prepared/` and `assembly_assessment/`
  explicitly (including failure PAF locations).
- `tests/test_dedup_contained.py` - standard-library unit tests.
- `tests/test_stage3_contig_policy.sh` (or equivalent Python subprocess test) -
  Stage 3 policy and generated-command tests.
- `changelogs.txt` - prepend the implemented change after validation.

No archived orchestrator or annotation scripts are reintroduced or modified.

## 8. Test matrix

### Utility tests

1. contained shorter query is dropped;
2. non-contained sequence is retained;
3. reverse-complement containment is detected;
4. exact equal-length duplicates retain the lexically chosen record
   deterministically, while merely near-identical equal-length records remain;
5. split alignments to one target can satisfy coverage (circular-origin case);
6. partial alignments spread across two targets do not satisfy containment;
7. overlapping query intervals are not double-counted;
8. identity and coverage boundary values behave exactly at the threshold;
9. `min_len=0` keeps short records; positive `min_len` reports a short drop;
10. duplicate IDs, malformed FASTA, empty sequences, malformed PAF, stale PAF
    lengths, unknown IDs, and invalid thresholds fail;
11. output order and retained headers are unchanged; and
12. report counts reconcile with output FASTA counts.

### Stage policy tests

1. all assemblies at or below 25: minimap2/dedup are never invoked;
2. one assembly above 25: only that file is processed;
3. `--skip-contained-dedup`: all inputs are copied and assessment still runs;
4. raw mean above 25, post mean at/below 25: explicit value 25 is passed to
   both compress and cluster;
5. post mean above 25 without override: stop before compress and print the
   required integer;
6. sufficient explicit override: identical value reaches both commands;
7. insufficient or invalid override: early failure;
8. failed/empty assembler outputs are excluded and reported;
9. originals remain byte-identical; only prepared headers receive weights;
10. stale manifest causes rebuild; matching manifest permits reuse;
11. resume from trim/later does not retroactively alter inputs; and
12. dry-run remains runnable with the generated mock environment.

### Laptop validation

The laptop gate is limited to cheap deterministic checks:

```bash
python3 -m unittest discover -s tests -v
python3 -m py_compile dedup_contained.py tests/test_dedup_contained.py \
  tests/test_stage3_contig_policy.py utils/generate_mock_environment.py
bash -n 03_autocycler_assemble.sh
python3 utils/generate_mock_environment.py
# run focused Stage 3 mock/policy checks
python3 utils/generate_mock_environment.py --clean
git diff --check
```

### HPC validation

Use the real HPC for the acceptance run with installed Autocycler and minimap2:

1. a normal isolate cohort where every assembly is below the trigger;
2. a cohort with one intentionally fragmented/duplicated assembly;
3. a cohort whose post-cleanup mean remains above 25;
4. an explicit override run for case 3, with wall time and peak RAM captured;
5. review every dropped contig against its named containing target, especially
   circular sequences and small plasmid candidates;
6. compare Autocycler cluster membership and consensus output with/without
   cleanup; and
7. confirm that increasing the guard on the tested hardware has acceptable
   compress/cluster time and memory before considering any future adaptive
   auto-raise mode.

Do not claim biological or performance validation from laptop mocks.

## 9. Implementation phases and gates

### Phase 1 - Harden the utility

Implement strict parsing, single-target grouping, safe defaults, and the event
report. Gate on the full utility test matrix and Python compilation.

### Phase 2 - Add preparation and policy calculation

Integrate post-assembly counting, conditional invocation, prepared copies,
manifesting, and max-contigs calculation into Stage 3. Gate on shell syntax and
focused synthetic Stage 3 tests.

### Phase 3 - CLI, config, resume, and curation

Add validated options/config, explicit override behavior, resume invalidation,
and human-readable diagnostics. Gate on generated command assertions proving
compress and cluster receive the same value.

### Phase 4 - Mock environment and documentation

Add realistic mock outputs, update user documentation and ignored runtime
paths, prepend the changelog, and run all laptop-safe checks.

### Phase 5 - HPC acceptance

Run the real-data matrix, inspect all removed sequences, record resource usage,
and tune only the per-file trigger or containment thresholds if the evidence
supports it. Do not change the default Autocycler guard based only on a run that
finishes successfully.

## 10. Acceptance criteria

Implementation is complete only when:

- every valid assembly has recorded pre/post counts;
- `dedup_contained.py` runs only for assemblies above the configured per-file
  trigger;
- containment requires one qualifying target;
- automatic mode never drops a contig solely for being short;
- original assembler FASTAs remain unchanged;
- every removal is attributable in a TSV report;
- both Autocycler commands receive the same explicit `--max_contigs` value;
- the default guard is not silently raised;
- post-mean failure provides the exact minimum override and stops before
  compression;
- resume cannot unknowingly mix stale prepared assemblies with an existing
  Autocycler graph;
- laptop unit/syntax/mock gates pass; and
- the real HPC acceptance matrix passes before the feature is described as
  validated on production data.

## 11. Primary upstream references

- Autocycler CLI source (`compress` and `cluster` `--max_contigs` definitions):
  https://github.com/rrwick/Autocycler/blob/main/src/main.rs
- Autocycler cluster source (mean-contigs check):
  https://github.com/rrwick/Autocycler/blob/main/src/cluster.rs
- Autocycler release history (the guard was introduced to catch obviously bad
  input that could make clustering hang):
  https://github.com/rrwick/Autocycler/releases
- Current Ryan Wick automated pipeline, for upstream step ordering:
  https://github.com/rrwick/Autocycler/blob/main/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick/autocycler_full.sh
- Minimap2 assembly alignment, secondary-mapping, and PAF documentation:
  https://github.com/lh3/minimap2
- Minimap2 self-homology command guidance (`-DP`):
  https://github.com/lh3/minimap2/blob/master/cookbook.md
