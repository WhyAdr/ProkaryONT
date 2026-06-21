# Filtering Approach Overview — Future Development Reference

> **Status: DEFERRED**
> This document captures the design decisions and architectural plan for
> enhancing the ProkaryONT read-filtering pipeline with Chopper, Fastplong,
> and Fastcat. Implementation is deferred until the internal filtering logic
> of each tool has been thoroughly reviewed and benchmarked against
> prokaryotic long-read datasets.

---

## 1. Current Architecture

The filtering step in
[02_filter_assemble.sh](02_filter_assemble.sh) (lines 173–193) uses a
two-stage serial approach:

```text
Raw FASTQ ──▶ Seqkit (true Phred Q-score filter) ──▶ tmpfile ──▶ Filtlong (length + priority subset)
```

### Why Seqkit before Filtlong (not redundant)

- **Filtlong's `--min_mean_q`** does _not_ perform standard average
  Phred-score filtering. It calculates a custom metric based on the Z-score
  distribution of per-base qualities, meaning a few high-quality bases can
  mask an otherwise poor read. Seqkit's `--min-qual` enforces a true
  arithmetic mean Phred threshold.
- **Filtlong cannot read from stdin.** Depending on compile-time choices it
  will silently fail or crash. The intermediate temp file is intentional and
  has been validated through trial and error.
- **Storage/memory overhead is acceptable** — the pipeline runs on HPC with
  ample scratch space.

### Post-filter QC

After filtering, NanoPlot runs on the filtered reads to produce a visual QC
report (`NanoPlot_FiltPol/`) that operators compare against the raw-read
NanoPlot report.

---

## 2. Proposed Enhancement (Future)

### Target architecture

```text
Raw FASTQ ──▶ Chopper (end-trim) ──▶ Fastplong (Q-filter + low-complexity) ──▶ Filtlong (length + priority subset)
```

With Fastcat added to Stage 1 QC for instant per-read metrics.

### Resolved Design Decisions

| Question | Decision |
| :--- | :--- |
| **Tool enforcement** | Soft dependencies — graceful fallback to existing Seqkit + Filtlong chain if Chopper/Fastplong/Fastcat are not installed |
| **Post-filter NanoPlot** | Keep running NanoPlot on filtered reads even when Fastplong produces its own QC report — run both |
| **Fastplong low-complexity** | Expose both `-y` (enable/disable) and `-Y` (threshold value) via CLI flags and config keys so the user can fine-tune |
| **Fastplong polyX trimming** | Expose `--trim_poly_x` as an opt-in CLI flag and config key |
| **Chopper crop defaults** | `--headcrop 50 --tailcrop 30` (configurable via `pipeline.conf`) |

---

## 3. Tool-by-Tool Assessment

### 3.1 Fastcat (Stage 1 QC addition)

**Purpose**: Instant per-read QC metrics before the slow NanoPlot runs.

**Risk level**: Very low — read-only, does not modify data.

**Value**: Produces structured `per_read_stats.tsv.gz` and a compact summary.
Gives the operator immediate numeric feedback (read count, total bases, N50,
mean Q) in the pipeline log. Can also be used for automated go/no-go checks
downstream.

**Complementarity with NanoPlot**: Not redundant. NanoPlot produces
publication-quality plots; Fastcat produces machine-readable metrics and runs
orders of magnitude faster. Both should coexist.

**Status**: Ready to implement. Low risk. Could be added independently.

---

### 3.2 Chopper (Stage 2, Step 2a — end-trimming)

**Purpose**: Crop low-quality bases from read ends caused by adapter
translocation and pore entry/exit artifacts.

**Key flags under review**:
- `--headcrop N` / `--tailcrop N` — Fixed-length cropping from 5′/3′ ends
- `--quality Q` — Sliding-window quality trimming (aggressive — needs
  benchmarking)
- `--minlength N` — Minimum read length after trimming
- `--maxlength N` — Maximum read length filter
- `--contam FILE` — Contamination screening (e.g., lambda spike-in)

**Concerns requiring investigation**:
- **Sliding-window quality trimming** (`--quality`): The window size and
  step-down behavior need characterization. Aggressive trimming can shorten
  reads substantially, reducing the overlap information that OLC assemblers
  rely on. For the initial implementation, we plan to use only fixed
  `--headcrop`/`--tailcrop` and _not_ enable quality-based window trimming.
- **Interaction with Filtlong**: If Chopper shortens reads below Filtlong's
  `--min_length`, those reads are silently discarded. The interaction between
  crop values and minimum length thresholds needs testing.
- **Read identity preservation**: Chopper does not modify read IDs. Verified.

**Status**: Deferred. Needs benchmarking of fixed-crop vs. window-trim
on representative prokaryotic ONT datasets to quantify read loss.

---

### 3.3 Fastplong (Stage 2, Step 2b — quality + low-complexity filtering)

**Purpose**: Replace Seqkit for Q-score filtering while adding entropy-based
low-complexity read removal and generating a QC report as a side effect.

**Key flags under review**:
- `-q N` — Minimum average quality score (true Phred mean — verify this
  matches Seqkit's behavior)
- `-y` — Enable low-complexity filter
- `-Y N` — Low-complexity threshold (default: 30 = 30% of bases are
  low-complexity → discard). The entropy calculation method needs review.
- `--trim_poly_x` — Trim poly-A/T/G/C tails
- `--length_required N` — Minimum read length
- `--html FILE` / `--json FILE` — QC report outputs

**Concerns requiring investigation**:
- **Q-score calculation method**: Does Fastplong's `-q` flag compute a true
  arithmetic mean of Phred scores (like Seqkit), or does it use a different
  formula? If it differs, we cannot simply swap it in as a Seqkit
  replacement. This is the single most critical question.
- **Low-complexity definition**: What entropy metric does Fastplong use? Is
  it Shannon entropy on dinucleotides, trinucleotides, or raw base
  composition? Prokaryotic genomes with extreme GC bias (e.g., _Streptomyces_
  at 70% GC) or repeat-rich regions (IS elements, rRNA operons) could
  trigger false positives.
- **PolyX trimming aggressiveness**: How many consecutive identical bases
  trigger trimming? ONT homopolymer errors are common — we need to ensure
  the trimmer doesn't clip legitimate biological homopolymer runs.
- **Interaction with downstream Filtlong**: If Fastplong trims reads shorter,
  the same crop-vs-min-length interaction concern as Chopper applies.

**Status**: Deferred. The Q-score calculation and low-complexity entropy
method must be verified against source code before adoption.

---

## 4. Proposed Config Additions (for future implementation)

```bash
# --- Chopper end-trimming (FUTURE — not yet active) --------------------------
# Bases to crop from the start/end of each read
chopper_headcrop=50
chopper_tailcrop=30

# --- Fastplong preprocessing (FUTURE — not yet active) ----------------------
# Enable low-complexity read filtering (entropy-based)
fastplong_low_complexity=true
# Threshold for low-complexity filter (default: 30)
# fastplong_low_complexity_threshold=30
# Enable poly-A/T/G/C tail trimming
# fastplong_trim_polyx=true
```

---

## 5. Proposed CLI Additions (for future implementation)

### 02_filter_assemble.sh

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--headcrop N` | int | 50 | Chopper: bases to crop from read start |
| `--tailcrop N` | int | 30 | Chopper: bases to crop from read end |
| `--no-low-complexity` | bool | off | Disable Fastplong low-complexity filter |
| `--low-complexity-threshold N` | int | 30 | Fastplong `-Y` threshold value |
| `--trim-polyx` | bool | off | Enable Fastplong poly-X tail trimming |

---

## 6. Backward Compatibility

All new tools are **soft dependencies**. If not installed:

| Missing Tool | Fallback Behavior |
| :--- | :--- |
| `chopper` | End-trimming step skipped; reads pass through untrimmed |
| `fastplong` | Falls back to `seqkit seq --min-qual` for Q-score filtering |
| `fastcat` | Per-read metric summary skipped; NanoPlot still runs |

The existing Seqkit + Filtlong chain remains fully functional and is the
default on environments without the new tools.

---

## 7. Next Steps (when resuming this work)

1. **Audit Fastplong source code** — Verify Q-score calculation in
   `src/qualityfilter.cpp` or equivalent. Compare output against Seqkit on
   the same dataset.
2. **Audit Chopper source code** — Characterize window-trim behavior and
   confirm fixed-crop is independent of quality scores.
3. **Benchmark on representative datasets** — Run Chopper + Fastplong +
   Filtlong vs. current Seqkit + Filtlong on 2–3 prokaryotic ONT datasets
   (varying GC content, genome size, coverage depth). Compare:
   - Read count and total bases retained at each stage
   - Final assembly contiguity (N50, number of contigs)
   - Assembly accuracy (QV from Merqury, BUSCO completeness)
4. **Implement Fastcat first** — Lowest risk, highest immediate value. Can
   be added to `01_qc_estimate.sh` independently of the filtering changes.
