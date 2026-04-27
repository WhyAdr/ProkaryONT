#!/usr/bin/env python3
"""
merge_gene_predictions.py — Merge gene predictions from multiple callers
                            into a consensus gene set.

Usage (called by 06a_predict_genes.sh):
    python3 merge_gene_predictions.py \
        --assembly genome.fasta \
        --pyrodigal-gff pyrodigal.gff \
        --glimmer-predict glimmer.predict \
        --glimmer-detail glimmer.detail \
        --genemark-gff genemark.gff \
        --merge-mode majority \
        --overlap-threshold 0.80 \
        --output-dir consensus/ \
        --log-file gene_merge.log

Outputs:
    consensus_genes.gff3      Consensus gene set in GFF3 format
    consensus_proteins.faa    Translated protein sequences
    consensus_cds.fna         Nucleotide CDS sequences
    gene_merge.log            Comprehensive decision log
"""

import argparse
import logging
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    print("ERROR: BioPython is required. Install with: pip install biopython",
          file=sys.stderr)
    sys.exit(1)


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class GenePrediction:
    """A single gene prediction from any caller."""
    seqid: str
    start: int           # 1-based, inclusive
    end: int             # 1-based, inclusive
    strand: str          # '+' or '-'
    caller: str          # 'pyrodigal', 'glimmer3', 'genemark'
    score: float = 0.0
    gene_id: str = ""
    attributes: Dict[str, str] = field(default_factory=dict)

    @property
    def length(self) -> int:
        return abs(self.end - self.start) + 1

    def overlaps(self, other: "GenePrediction", threshold: float) -> bool:
        """Check if two predictions overlap by at least threshold (reciprocal)."""
        if self.seqid != other.seqid or self.strand != other.strand:
            return False
        overlap_start = max(self.start, other.start)
        overlap_end = min(self.end, other.end)
        if overlap_start > overlap_end:
            return False
        overlap_len = overlap_end - overlap_start + 1
        recip_a = overlap_len / self.length
        recip_b = overlap_len / other.length
        return min(recip_a, recip_b) >= threshold


@dataclass
class ConsensusGene:
    """A merged consensus gene with evidence from multiple callers."""
    seqid: str
    start: int
    end: int
    strand: str
    callers: List[str]
    start_evidence: str   # 'majority', 'tiebreak', 'single'
    gene_id: str = ""
    predictions: List[GenePrediction] = field(default_factory=list)


# =============================================================================
# GFF Parsers
# =============================================================================

def parse_pyrodigal_gff(filepath: str) -> List[GenePrediction]:
    """Parse Pyrodigal GFF3 output."""
    predictions = []
    gene_idx = 0
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "CDS":
                continue
            gene_idx += 1
            attrs = _parse_gff_attributes(parts[8])
            predictions.append(GenePrediction(
                seqid=parts[0],
                start=int(parts[3]),
                end=int(parts[4]),
                strand=parts[6],
                caller="pyrodigal",
                score=float(parts[5]) if parts[5] != "." else 0.0,
                gene_id=attrs.get("ID", f"pyrodigal_{gene_idx}"),
                attributes=attrs,
            ))
    return predictions


def parse_glimmer_predict(predict_file: str,
                          detail_file: str) -> List[GenePrediction]:
    """Parse Glimmer3 .predict output.

    Glimmer3 predict format per sequence:
        >seqid
        orf00001  start  end  frame  score
    """
    predictions = []
    current_seqid = None
    gene_idx = 0

    with open(predict_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_seqid = line[1:].split()[0]
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            gene_idx += 1
            orf_id = parts[0]
            raw_start = int(parts[1])
            raw_end = int(parts[2])
            score = float(parts[4])

            # Glimmer3 uses signed coordinates: negative start = complement
            if raw_start > raw_end:
                strand = "-"
                start = raw_end
                end = raw_start
            else:
                strand = "+"
                start = raw_start
                end = raw_end

            # Clamp start to 1 (Glimmer sometimes outputs 0 or negative)
            start = max(1, start)

            predictions.append(GenePrediction(
                seqid=current_seqid,
                start=start,
                end=end,
                strand=strand,
                caller="glimmer3",
                score=score,
                gene_id=f"glimmer3_{gene_idx}",
            ))

    return predictions


def parse_genemark_gff(filepath: str) -> List[GenePrediction]:
    """Parse GeneMarkS-2 GFF output (standard GFF3 or GTF)."""
    predictions = []
    gene_idx = 0
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] not in ("CDS", "gene"):
                continue
            gene_idx += 1
            attrs = _parse_gff_attributes(parts[8])
            predictions.append(GenePrediction(
                seqid=parts[0],
                start=int(parts[3]),
                end=int(parts[4]),
                strand=parts[6],
                caller="genemark",
                score=float(parts[5]) if parts[5] != "." else 0.0,
                gene_id=attrs.get("ID",
                                  attrs.get("gene_id",
                                            f"genemark_{gene_idx}")),
                attributes=attrs,
            ))
    return predictions


def _parse_gff_attributes(attr_str: str) -> Dict[str, str]:
    """Parse GFF3 attribute column (key=value;key=value)."""
    attrs = {}
    for item in attr_str.split(";"):
        item = item.strip()
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


# =============================================================================
# Merge Logic
# =============================================================================

def group_by_locus(all_predictions: List[GenePrediction],
                   threshold: float) -> List[List[GenePrediction]]:
    """Group overlapping predictions from different callers into loci."""
    if not all_predictions:
        return []

    # Sort by seqid, then start position
    sorted_preds = sorted(all_predictions, key=lambda p: (p.seqid, p.start))

    loci: List[List[GenePrediction]] = []
    current_locus: List[GenePrediction] = [sorted_preds[0]]

    for pred in sorted_preds[1:]:
        # Check overlap with any prediction in the current locus
        merged = False
        for existing in current_locus:
            if pred.overlaps(existing, threshold):
                merged = True
                break
        if merged:
            current_locus.append(pred)
        else:
            loci.append(current_locus)
            current_locus = [pred]

    loci.append(current_locus)  # Don't forget the last group
    return loci


def merge_majority_vote(loci: List[List[GenePrediction]],
                        total_callers: int,
                        logger: logging.Logger) -> List[ConsensusGene]:
    """Majority vote: keep loci predicted by >= 2 callers (or >=1 if only 1
    caller available). Uses Pyrodigal as tiebreaker for start codons."""
    consensus = []
    min_votes = max(2, (total_callers + 1) // 2) if total_callers > 1 else 1

    for locus in loci:
        callers = list({p.caller for p in locus})
        if len(callers) < min_votes:
            for p in locus:
                logger.info(
                    "DISCARDED [majority] %s:%d-%d(%s) caller=%s "
                    "votes=%d/%d (need %d)",
                    p.seqid, p.start, p.end, p.strand,
                    p.caller, len(callers), total_callers, min_votes)
            continue

        gene = _resolve_locus(locus, callers, logger)
        consensus.append(gene)

    return consensus


def merge_union(loci: List[List[GenePrediction]],
                logger: logging.Logger) -> List[ConsensusGene]:
    """Union: keep every unique locus. Still resolve start codons by vote."""
    consensus = []
    for locus in loci:
        callers = list({p.caller for p in locus})
        gene = _resolve_locus(locus, callers, logger)
        consensus.append(gene)
    return consensus


def _resolve_locus(locus: List[GenePrediction],
                   callers: List[str],
                   logger: logging.Logger) -> ConsensusGene:
    """Resolve a single locus into one consensus gene.

    Rules:
        - Stop codon: majority vote among callers (all usually agree)
        - Start codon: majority vote, Pyrodigal tiebreaker
        - Coordinates: from the winning start to the winning stop
    """
    # All share the same seqid and strand (enforced by overlap grouping)
    seqid = locus[0].seqid
    strand = locus[0].strand

    # --- Resolve stop codon (end for +, start for -) ---
    if strand == "+":
        end_votes = defaultdict(list)
        for p in locus:
            end_votes[p.end].append(p.caller)
        winning_end = max(end_votes, key=lambda e: len(end_votes[e]))
        end = winning_end
    else:
        start_votes_stop = defaultdict(list)
        for p in locus:
            start_votes_stop[p.start].append(p.caller)
        winning_start_stop = min(start_votes_stop,
                                 key=lambda s: -len(start_votes_stop[s]))
        # For minus strand, the "stop" is the lowest coordinate

    # --- Resolve start codon ---
    if strand == "+":
        start_votes = defaultdict(list)
        for p in locus:
            start_votes[p.start].append(p.caller)
    else:
        start_votes = defaultdict(list)
        for p in locus:
            start_votes[p.end].append(p.caller)

    # Find max vote count
    max_vote = max(len(v) for v in start_votes.values())
    candidates = [pos for pos, voters in start_votes.items()
                  if len(voters) == max_vote]

    if len(candidates) == 1:
        winning_start = candidates[0]
        start_evidence = "majority"
    else:
        # Tiebreak: prefer Pyrodigal's prediction
        pyrodigal_pred = [p for p in locus if p.caller == "pyrodigal"]
        if pyrodigal_pred:
            if strand == "+":
                winning_start = pyrodigal_pred[0].start
            else:
                winning_start = pyrodigal_pred[0].end
            start_evidence = "tiebreak"
        else:
            # No Pyrodigal — pick the most upstream (smallest for +)
            winning_start = min(candidates) if strand == "+" else max(candidates)
            start_evidence = "tiebreak"

    if strand == "+":
        final_start = winning_start
        final_end = end
    else:
        final_start = winning_start_stop if 'winning_start_stop' in dir() else min(p.start for p in locus)
        final_end = winning_start

    # Build consensus gene
    gene = ConsensusGene(
        seqid=seqid,
        start=min(final_start, final_end),
        end=max(final_start, final_end),
        strand=strand,
        callers=sorted(callers),
        start_evidence=start_evidence,
        predictions=locus,
    )

    logger.info(
        "KEPT %s:%d-%d(%s) callers=%s start_evidence=%s",
        gene.seqid, gene.start, gene.end, gene.strand,
        ",".join(gene.callers), gene.start_evidence,
    )

    return gene


# =============================================================================
# Output Writers
# =============================================================================

def write_consensus_gff(genes: List[ConsensusGene],
                        outfile: str) -> None:
    """Write consensus genes as GFF3."""
    with open(outfile, "w") as f:
        f.write("##gff-version 3\n")
        for i, gene in enumerate(genes, 1):
            gene.gene_id = f"consensus_{i:05d}"
            attrs = (
                f"ID={gene.gene_id};"
                f"callers={','.join(gene.callers)};"
                f"start_evidence={gene.start_evidence};"
                f"num_callers={len(gene.callers)}"
            )
            f.write(f"{gene.seqid}\tconsensus\tCDS\t{gene.start}\t"
                    f"{gene.end}\t.\t{gene.strand}\t0\t{attrs}\n")


def extract_sequences(genes: List[ConsensusGene],
                      assembly: Dict[str, str],
                      faa_file: str,
                      fna_file: str) -> None:
    """Extract CDS nucleotide and protein sequences from assembly."""
    with open(fna_file, "w") as fna, open(faa_file, "w") as faa:
        for gene in genes:
            seq_record = assembly.get(gene.seqid)
            if seq_record is None:
                continue

            # Extract nucleotide sequence (0-based slicing)
            nuc_seq = seq_record[gene.start - 1:gene.end]

            if gene.strand == "-":
                nuc_seq = str(Seq(nuc_seq).reverse_complement())
            else:
                nuc_seq = str(nuc_seq)

            # Translate
            try:
                protein = str(Seq(nuc_seq).translate(table=11, to_stop=True))
            except Exception:
                protein = ""

            header = (f">{gene.gene_id} {gene.seqid}:{gene.start}-{gene.end}"
                      f"({gene.strand}) callers={','.join(gene.callers)}")

            fna.write(f"{header}\n{nuc_seq}\n")
            faa.write(f"{header}\n{protein}\n")


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Merge multi-caller gene predictions into consensus set")
    parser.add_argument("--assembly", required=True,
                        help="Assembly FASTA file")
    parser.add_argument("--pyrodigal-gff", required=True,
                        help="Pyrodigal GFF3 output")
    parser.add_argument("--glimmer-predict",
                        help="Glimmer3 .predict file")
    parser.add_argument("--glimmer-detail",
                        help="Glimmer3 .detail file")
    parser.add_argument("--genemark-gff",
                        help="GeneMarkS-2 GFF file")
    parser.add_argument("--merge-mode", default="majority",
                        choices=["majority", "union"],
                        help="Merge strategy (default: majority)")
    parser.add_argument("--overlap-threshold", type=float, default=0.80,
                        help="Reciprocal overlap threshold (default: 0.80)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for consensus files")
    parser.add_argument("--log-file", default="gene_merge.log",
                        help="Path for decision log")

    args = parser.parse_args()

    # --- Set up logging -------------------------------------------------------
    os.makedirs(args.output_dir, exist_ok=True)

    logger = logging.getLogger("gene_merge")
    logger.setLevel(logging.DEBUG)

    fh = logging.FileHandler(args.log_file, mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    logger.addHandler(ch)

    # --- Load predictions -----------------------------------------------------
    all_predictions: List[GenePrediction] = []
    active_callers = 0

    logger.info("Loading Pyrodigal predictions: %s", args.pyrodigal_gff)
    pyrodigal_preds = parse_pyrodigal_gff(args.pyrodigal_gff)
    all_predictions.extend(pyrodigal_preds)
    active_callers += 1
    logger.info("  Pyrodigal: %d CDS features", len(pyrodigal_preds))

    if args.glimmer_predict:
        logger.info("Loading Glimmer3 predictions: %s", args.glimmer_predict)
        glimmer_preds = parse_glimmer_predict(
            args.glimmer_predict, args.glimmer_detail or "")
        all_predictions.extend(glimmer_preds)
        active_callers += 1
        logger.info("  Glimmer3: %d CDS features", len(glimmer_preds))

    if args.genemark_gff:
        logger.info("Loading GeneMarkS-2 predictions: %s", args.genemark_gff)
        genemark_preds = parse_genemark_gff(args.genemark_gff)
        all_predictions.extend(genemark_preds)
        active_callers += 1
        logger.info("  GeneMarkS-2: %d CDS features", len(genemark_preds))

    logger.info("Total predictions: %d from %d callers",
                len(all_predictions), active_callers)

    # --- Group by locus -------------------------------------------------------
    logger.info("Grouping predictions by locus (overlap threshold: %.2f)",
                args.overlap_threshold)
    loci = group_by_locus(all_predictions, args.overlap_threshold)
    logger.info("Found %d distinct loci", len(loci))

    # --- Merge ----------------------------------------------------------------
    logger.info("Merge mode: %s", args.merge_mode)
    if args.merge_mode == "majority":
        consensus = merge_majority_vote(loci, active_callers, logger)
    else:
        consensus = merge_union(loci, logger)

    logger.info("Consensus gene set: %d genes", len(consensus))

    # --- Load assembly for sequence extraction --------------------------------
    logger.info("Loading assembly: %s", args.assembly)
    assembly_seqs = {}
    for record in SeqIO.parse(args.assembly, "fasta"):
        assembly_seqs[record.id] = str(record.seq)

    # --- Write outputs --------------------------------------------------------
    gff_out = os.path.join(args.output_dir, "consensus_genes.gff3")
    faa_out = os.path.join(args.output_dir, "consensus_proteins.faa")
    fna_out = os.path.join(args.output_dir, "consensus_cds.fna")

    write_consensus_gff(consensus, gff_out)
    extract_sequences(consensus, assembly_seqs, faa_out, fna_out)

    logger.info("Written: %s", gff_out)
    logger.info("Written: %s", faa_out)
    logger.info("Written: %s", fna_out)
    logger.info("Log:     %s", args.log_file)

    # --- Summary stats --------------------------------------------------------
    caller_counts = defaultdict(int)
    evidence_counts = defaultdict(int)
    for gene in consensus:
        for c in gene.callers:
            caller_counts[c] += 1
        evidence_counts[gene.start_evidence] += 1

    logger.info("=== MERGE SUMMARY ===")
    logger.info("Total consensus genes: %d", len(consensus))
    for caller, count in sorted(caller_counts.items()):
        logger.info("  Supported by %s: %d", caller, count)
    for evidence, count in sorted(evidence_counts.items()):
        logger.info("  Start resolved by %s: %d", evidence, count)


if __name__ == "__main__":
    main()
