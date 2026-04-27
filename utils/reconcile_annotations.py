#!/usr/bin/env python3
"""
reconcile_annotations.py — Reconcile multi-source annotations into unified output.

Usage (called by 06e_reconcile_merge.sh):
    python3 reconcile_annotations.py \
        --assembly genome.fasta \
        --consensus-gff consensus_genes.gff3 \
        --bakta-dir bakta/ \
        --emapper-dir eggnog/ \
        --interproscan-dir interproscan/ \
        --kofamscan-dir kofamscan/ \
        --diamond-dir diamond/ \
        --dram-dir dram/ \
        --rast-gto rast.gto \
        --phanotate-gff phanotate.gff3 \
        --overlap-threshold 0.80 \
        --sample-name MyBacteria \
        --output-dir final/ \
        --log-file reconciliation.log

Outputs:
    final_annotation.gff3     Unified GFF3 with all evidence in attributes
    annotation_matrix.tsv     Gene × annotation source table
    final_annotation.gbk      GenBank format
    summary.json              Per-isolate machine-readable summary
    reconciliation.log        Comprehensive decision log
"""

import argparse
import csv
import json
import logging
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
except ImportError:
    print("ERROR: BioPython is required. Install with: pip install biopython",
          file=sys.stderr)
    sys.exit(1)


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class AnnotatedGene:
    """A gene with merged annotations from all evidence sources."""
    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str

    # Gene prediction metadata
    callers: List[str] = field(default_factory=list)
    start_evidence: str = ""
    merge_source: str = ""  # 'consensus', 'bakta_only', 'rast_only', etc.
    region_type: str = ""   # '' or 'prophage'

    # Non-competing annotation layers
    product: str = "hypothetical protein"
    product_source: str = ""

    # KEGG
    kofamscan_ko: str = ""
    kofamscan_score: str = ""
    emapper_ko: str = ""

    # COG
    emapper_cog: str = ""
    bakta_cog: str = ""

    # GO
    interproscan_go: str = ""
    emapper_go: str = ""
    bakta_go: str = ""

    # Pfam
    interproscan_pfam: str = ""
    bakta_pfam: str = ""

    # EC number
    kofamscan_ec: str = ""
    interproscan_ec: str = ""
    emapper_ec: str = ""

    # Unique databases
    rast_subsystem: str = ""
    rast_role: str = ""
    dram_merops: str = ""
    dram_cazy: str = ""
    dram_vogdb: str = ""
    dram_kegg: str = ""
    dram_pfam: str = ""

    # DIAMOND/SwissProt
    swissprot_hit: str = ""
    swissprot_evalue: str = ""
    swissprot_product: str = ""

    # Bakta specifics
    bakta_product: str = ""
    bakta_dbxref: str = ""

    @property
    def length(self) -> int:
        return abs(self.end - self.start) + 1


# =============================================================================
# Parsers for each annotation source
# =============================================================================

def parse_consensus_gff(filepath: str) -> List[AnnotatedGene]:
    """Parse consensus GFF3 from merge_gene_predictions.py."""
    genes = []
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = _parse_attrs(parts[8])
            gene = AnnotatedGene(
                gene_id=attrs.get("ID", ""),
                seqid=parts[0],
                start=int(parts[3]),
                end=int(parts[4]),
                strand=parts[6],
                callers=attrs.get("callers", "").split(","),
                start_evidence=attrs.get("start_evidence", ""),
                merge_source="consensus",
            )
            genes.append(gene)
    return genes


def parse_bakta_gff(bakta_dir: str) -> List[Dict]:
    """Parse Bakta GFF3 for gene coordinates and annotations."""
    genes = []
    gff_files = list(Path(bakta_dir).glob("*.gff3"))
    if not gff_files:
        return []
    with open(gff_files[0]) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = _parse_attrs(parts[8])
            genes.append({
                "seqid": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "product": attrs.get("product", attrs.get("Name", "")),
                "dbxref": attrs.get("Dbxref", ""),
                "attrs": attrs,
            })
    return genes


def parse_emapper(emapper_dir: str) -> Dict[str, Dict]:
    """Parse eggNOG-mapper annotations file."""
    results = {}
    ann_files = list(Path(emapper_dir).glob("*.emapper.annotations"))
    if not ann_files:
        return {}
    with open(ann_files[0]) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 22:
                continue
            gene_id = parts[0]
            results[gene_id] = {
                "ko": parts[11] if parts[11] != "-" else "",
                "cog": parts[6] if parts[6] != "-" else "",
                "go": parts[9] if parts[9] != "-" else "",
                "ec": parts[10] if parts[10] != "-" else "",
                "pfam": parts[20] if len(parts) > 20 and parts[20] != "-"
                        else "",
                "product": parts[7] if parts[7] != "-" else "",
                "tax_scope": parts[16] if len(parts) > 16
                             and parts[16] != "-" else "",
            }
    return results


def parse_interproscan(ips_dir: str) -> Dict[str, Dict]:
    """Parse InterProScan TSV output."""
    results: Dict[str, Dict] = defaultdict(lambda: {
        "pfam": [], "go": [], "ec": [], "pathways": []})
    tsv_files = list(Path(ips_dir).glob("*.tsv"))
    if not tsv_files:
        return {}
    with open(tsv_files[0]) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 12:
                continue
            gene_id = row[0]
            db = row[3]        # Source database
            sig_acc = row[4]   # Signature accession
            sig_desc = row[5]  # Signature description

            if db == "Pfam":
                results[gene_id]["pfam"].append(sig_acc)

            # GO terms (column 13 if present)
            if len(row) > 13 and row[13]:
                for go_term in row[13].split("|"):
                    if go_term and go_term not in results[gene_id]["go"]:
                        results[gene_id]["go"].append(go_term)

            # Pathways (column 14 if present)
            if len(row) > 14 and row[14]:
                results[gene_id]["pathways"].append(row[14])

    # Flatten lists to strings
    flat = {}
    for gene_id, data in results.items():
        flat[gene_id] = {
            "pfam": ",".join(sorted(set(data["pfam"]))),
            "go": ",".join(sorted(set(data["go"]))),
            "pathways": ",".join(data["pathways"]),
        }
    return flat


def parse_kofamscan(kofam_dir: str) -> Dict[str, Dict]:
    """Parse KofamScan detail-tsv output."""
    results = {}
    out_files = list(Path(kofam_dir).glob("*.txt"))
    if not out_files:
        return {}
    with open(out_files[0]) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Format: significance gene_id KO threshold score E-value definition
            parts = line.split(None, 6)
            if len(parts) < 6:
                continue
            sig = parts[0]       # '*' = significant
            gene_id = parts[1]
            ko = parts[2]
            score = parts[4]
            definition = parts[6] if len(parts) > 6 else ""

            # Only keep significant hits, or best hit per gene
            if sig == "*":
                results[gene_id] = {
                    "ko": ko,
                    "score": score,
                    "definition": definition,
                }
    return results


def parse_diamond_swissprot(diamond_dir: str) -> Dict[str, Dict]:
    """Parse DIAMOND blastp TSV (outfmt 6)."""
    results = {}
    tsv_files = list(Path(diamond_dir).glob("*.tsv"))
    if not tsv_files:
        return {}
    with open(tsv_files[0]) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 7:
                continue
            gene_id = row[0]
            if gene_id not in results:  # Keep top hit only
                sseqid = row[1]
                evalue = row[4]
                # Extract product from stitle (column 7)
                stitle = row[6] if len(row) > 6 else ""
                # SwissProt stitle format: "ProtName OS=Organism ..."
                product = re.split(r"\s+OS=", stitle)[0] if stitle else ""
                results[gene_id] = {
                    "hit": sseqid,
                    "evalue": evalue,
                    "product": product,
                }
    return results


def parse_dram(dram_dir: str) -> Dict[str, Dict]:
    """Parse DRAM annotations.tsv."""
    results = {}
    ann_file = Path(dram_dir) / "annotations.tsv"
    if not ann_file.exists():
        return {}
    with open(ann_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene_id = row.get("", row.get("gene_id", ""))
            results[gene_id] = {
                "kegg": row.get("kegg_id", ""),
                "pfam": row.get("pfam_hits", ""),
                "cazy": row.get("cazy_hits", row.get("cazy_best_hit", "")),
                "merops": row.get("merops_id", row.get("merops_hits", "")),
                "vogdb": row.get("vogdb_id",
                                 row.get("vogdb_hits",
                                         row.get("vogdb_categories", ""))),
            }
    return results


def parse_phanotate_gff(filepath: str) -> List[Dict]:
    """Parse Phanotate GFF3 for prophage gene predictions."""
    genes = []
    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = _parse_attrs(parts[8])
            genes.append({
                "seqid": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "attrs": attrs,
            })
    return genes


def _parse_attrs(attr_str: str) -> Dict[str, str]:
    attrs = {}
    for item in attr_str.split(";"):
        item = item.strip()
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k] = v
    return attrs


# =============================================================================
# Coordinate Matching
# =============================================================================

def match_coordinates(consensus_genes: List[AnnotatedGene],
                      track_a_genes: List[Dict],
                      threshold: float,
                      source_name: str,
                      logger: logging.Logger) -> Dict[int, int]:
    """Match Track A genes to consensus genes by coordinate overlap.

    Returns: mapping of consensus_index → track_a_index
    """
    matches = {}
    for ci, cg in enumerate(consensus_genes):
        best_overlap = 0.0
        best_idx = -1
        for ai, ag in enumerate(track_a_genes):
            if cg.seqid != ag["seqid"] or cg.strand != ag["strand"]:
                continue
            ov_start = max(cg.start, ag["start"])
            ov_end = min(cg.end, ag["end"])
            if ov_start > ov_end:
                continue
            ov_len = ov_end - ov_start + 1
            cg_len = cg.end - cg.start + 1
            ag_len = ag["end"] - ag["start"] + 1
            recip = min(ov_len / cg_len, ov_len / ag_len)
            if recip >= threshold and recip > best_overlap:
                best_overlap = recip
                best_idx = ai

        if best_idx >= 0:
            matches[ci] = best_idx
            logger.debug(
                "MATCH %s→%s %s:%d-%d overlap=%.2f",
                consensus_genes[ci].gene_id, source_name,
                cg.seqid, cg.start, cg.end, best_overlap)

    logger.info("Matched %d/%d consensus genes to %s",
                len(matches), len(consensus_genes), source_name)
    return matches


# =============================================================================
# Product Name Resolution
# =============================================================================

PRODUCT_PRIORITY = [
    ("swissprot_product", "swissprot"),
    ("bakta_product", "bakta"),
    ("emapper_product", "eggnog"),
    ("rast_role", "rasttk"),
]


def resolve_product(gene: AnnotatedGene) -> Tuple[str, str]:
    """Resolve the best product name using priority hierarchy."""
    for attr_name, source in PRODUCT_PRIORITY:
        value = getattr(gene, attr_name, "")
        if value and value.lower() not in ("hypothetical protein", "", "-"):
            return value, source
    return "hypothetical protein", "none"


# =============================================================================
# Output Writers
# =============================================================================

def write_gff3(genes: List[AnnotatedGene], outfile: str) -> None:
    """Write final merged GFF3."""
    with open(outfile, "w") as f:
        f.write("##gff-version 3\n")
        for gene in genes:
            product, product_src = resolve_product(gene)
            attrs_parts = [
                f"ID={gene.gene_id}",
                f"product={product} [source:{product_src}]",
                f"merge_source={gene.merge_source}",
                f"callers={','.join(gene.callers)}",
                f"start_evidence={gene.start_evidence}",
            ]
            if gene.region_type:
                attrs_parts.append(f"region_type={gene.region_type}")

            # Non-competing layers
            for attr in [
                "kofamscan_ko", "kofamscan_score", "emapper_ko",
                "emapper_cog", "bakta_cog",
                "interproscan_go", "emapper_go",
                "interproscan_pfam", "bakta_pfam",
                "kofamscan_ec", "interproscan_ec", "emapper_ec",
                "rast_subsystem", "rast_role",
                "dram_merops", "dram_cazy", "dram_vogdb",
                "swissprot_hit", "swissprot_evalue",
                "bakta_dbxref",
            ]:
                value = getattr(gene, attr, "")
                if value:
                    attrs_parts.append(f"{attr}={value}")

            f.write(f"{gene.seqid}\tensemble\tCDS\t{gene.start}\t"
                    f"{gene.end}\t.\t{gene.strand}\t0\t"
                    f"{';'.join(attrs_parts)}\n")


def write_tsv(genes: List[AnnotatedGene], outfile: str) -> None:
    """Write annotation matrix as TSV."""
    columns = [
        "gene_id", "seqid", "start", "end", "strand", "length",
        "product", "product_source", "merge_source", "region_type",
        "callers", "start_evidence",
        "kofamscan_ko", "kofamscan_score", "emapper_ko",
        "emapper_cog", "bakta_cog",
        "interproscan_go", "emapper_go",
        "interproscan_pfam", "bakta_pfam",
        "kofamscan_ec", "interproscan_ec", "emapper_ec",
        "rast_subsystem", "rast_role",
        "dram_merops", "dram_cazy", "dram_vogdb",
        "swissprot_hit", "swissprot_evalue", "swissprot_product",
        "bakta_product", "bakta_dbxref",
    ]

    with open(outfile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for gene in genes:
            product, product_src = resolve_product(gene)
            row = {c: getattr(gene, c, "") for c in columns}
            row["length"] = gene.length
            row["product"] = product
            row["product_source"] = product_src
            row["callers"] = ",".join(gene.callers)
            writer.writerow(row)


def write_genbank(genes: List[AnnotatedGene], assembly_file: str,
                  outfile: str, sample_name: str) -> None:
    """Write GenBank output from merged annotations."""
    records = []
    assembly_seqs = {r.id: r for r in SeqIO.parse(assembly_file, "fasta")}

    # Group genes by seqid
    genes_by_seq: Dict[str, List[AnnotatedGene]] = defaultdict(list)
    for gene in genes:
        genes_by_seq[gene.seqid].append(gene)

    for seqid, seq_record in assembly_seqs.items():
        record = SeqRecord(
            seq_record.seq,
            id=seqid,
            name=sample_name,
            description=f"{sample_name} ensemble annotation",
        )
        record.annotations["molecule_type"] = "DNA"

        for gene in genes_by_seq.get(seqid, []):
            strand = 1 if gene.strand == "+" else -1
            location = FeatureLocation(gene.start - 1, gene.end,
                                       strand=strand)
            product, _ = resolve_product(gene)
            qualifiers = {
                "locus_tag": [gene.gene_id],
                "product": [product],
            }
            if gene.kofamscan_ko:
                qualifiers["note"] = [f"KEGG:{gene.kofamscan_ko}"]
            if gene.interproscan_go:
                qualifiers["db_xref"] = gene.interproscan_go.split(",")

            feature = SeqFeature(location, type="CDS",
                                 qualifiers=qualifiers)
            record.features.append(feature)

        records.append(record)

    SeqIO.write(records, outfile, "genbank")


def write_summary_json(genes: List[AnnotatedGene], sample_name: str,
                       outfile: str) -> None:
    """Write per-isolate summary JSON."""
    total = len(genes)
    annotated = sum(1 for g in genes
                    if resolve_product(g)[0] != "hypothetical protein")

    summary = {
        "sample": sample_name,
        "total_cds": total,
        "annotated_cds": annotated,
        "hypothetical_cds": total - annotated,
        "annotation_rate": round(annotated / total, 4) if total > 0 else 0,
        "sources": {
            "consensus_genes": sum(1 for g in genes
                                   if g.merge_source == "consensus"),
            "bakta_only": sum(1 for g in genes
                              if g.merge_source == "bakta_only"),
            "prophage": sum(1 for g in genes
                            if g.region_type == "prophage"),
        },
        "kegg_coverage": sum(1 for g in genes if g.kofamscan_ko),
        "cog_coverage": sum(1 for g in genes
                            if g.emapper_cog or g.bakta_cog),
        "go_coverage": sum(1 for g in genes
                           if g.interproscan_go or g.emapper_go),
        "pfam_coverage": sum(1 for g in genes
                             if g.interproscan_pfam or g.bakta_pfam),
        "seed_subsystems": sum(1 for g in genes if g.rast_subsystem),
        "dram_cazy": sum(1 for g in genes if g.dram_cazy),
        "dram_merops": sum(1 for g in genes if g.dram_merops),
        "dram_vogdb": sum(1 for g in genes if g.dram_vogdb),
        "swissprot_hits": sum(1 for g in genes if g.swissprot_hit),
    }

    with open(outfile, "w") as f:
        json.dump(summary, f, indent=2)


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Reconcile multi-source annotations into unified output")
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--consensus-gff", required=True)
    parser.add_argument("--bakta-dir", required=True)
    parser.add_argument("--rast-gto")
    parser.add_argument("--dram-dir")
    parser.add_argument("--emapper-dir")
    parser.add_argument("--interproscan-dir")
    parser.add_argument("--kofamscan-dir")
    parser.add_argument("--diamond-dir")
    parser.add_argument("--phanotate-gff")
    parser.add_argument("--overlap-threshold", type=float, default=0.80)
    parser.add_argument("--sample-name", default="MyBacteria")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-file", default="reconciliation.log")

    args = parser.parse_args()

    # --- Logging --------------------------------------------------------------
    os.makedirs(args.output_dir, exist_ok=True)

    logger = logging.getLogger("reconcile")
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

    # =========================================================================
    # 1. Load consensus gene set (backbone)
    # =========================================================================
    logger.info("Loading consensus genes: %s", args.consensus_gff)
    genes = parse_consensus_gff(args.consensus_gff)
    logger.info("  Loaded %d consensus genes", len(genes))

    # =========================================================================
    # 2. Match and merge Bakta annotations (Track A)
    # =========================================================================
    logger.info("Loading Bakta annotations: %s", args.bakta_dir)
    bakta_genes = parse_bakta_gff(args.bakta_dir)
    logger.info("  Bakta: %d CDS features", len(bakta_genes))

    bakta_matches = match_coordinates(
        genes, bakta_genes, args.overlap_threshold, "bakta", logger)

    for ci, ai in bakta_matches.items():
        bg = bakta_genes[ai]
        genes[ci].bakta_product = bg.get("product", "")
        genes[ci].bakta_dbxref = bg.get("dbxref", "")

    # Find Bakta-only genes (not in consensus)
    matched_bakta = set(bakta_matches.values())
    bakta_only_count = 0
    for ai, bg in enumerate(bakta_genes):
        if ai not in matched_bakta:
            bakta_only_count += 1
            gene = AnnotatedGene(
                gene_id=f"bakta_only_{bakta_only_count:05d}",
                seqid=bg["seqid"],
                start=bg["start"],
                end=bg["end"],
                strand=bg["strand"],
                merge_source="bakta_only",
                bakta_product=bg.get("product", ""),
                bakta_dbxref=bg.get("dbxref", ""),
            )
            genes.append(gene)
            logger.info("ADDED bakta_only %s:%d-%d(%s) product=%s",
                        gene.seqid, gene.start, gene.end, gene.strand,
                        gene.bakta_product)

    logger.info("  Added %d bakta-only genes", bakta_only_count)

    # =========================================================================
    # 3. Merge eggNOG-mapper (Track B)
    # =========================================================================
    if args.emapper_dir:
        logger.info("Loading eggNOG-mapper: %s", args.emapper_dir)
        emapper = parse_emapper(args.emapper_dir)
        logger.info("  eggNOG: %d annotated proteins", len(emapper))
        for gene in genes:
            if gene.gene_id in emapper:
                em = emapper[gene.gene_id]
                gene.emapper_ko = em.get("ko", "")
                gene.emapper_cog = em.get("cog", "")
                gene.emapper_go = em.get("go", "")
                gene.emapper_ec = em.get("ec", "")
                # Store eggNOG product for resolution
                setattr(gene, "emapper_product", em.get("product", ""))

    # =========================================================================
    # 4. Merge InterProScan (Track B)
    # =========================================================================
    if args.interproscan_dir:
        logger.info("Loading InterProScan: %s", args.interproscan_dir)
        ips = parse_interproscan(args.interproscan_dir)
        logger.info("  InterProScan: %d proteins with hits", len(ips))
        for gene in genes:
            if gene.gene_id in ips:
                ip = ips[gene.gene_id]
                gene.interproscan_pfam = ip.get("pfam", "")
                gene.interproscan_go = ip.get("go", "")

    # =========================================================================
    # 5. Merge KofamScan (Track B — authoritative KEGG KOs)
    # =========================================================================
    if args.kofamscan_dir:
        logger.info("Loading KofamScan: %s", args.kofamscan_dir)
        kofam = parse_kofamscan(args.kofamscan_dir)
        logger.info("  KofamScan: %d significant KO hits", len(kofam))
        for gene in genes:
            if gene.gene_id in kofam:
                kf = kofam[gene.gene_id]
                gene.kofamscan_ko = kf.get("ko", "")
                gene.kofamscan_score = kf.get("score", "")

    # =========================================================================
    # 6. Merge DIAMOND/SwissProt (Track B — optional)
    # =========================================================================
    if args.diamond_dir:
        logger.info("Loading DIAMOND/SwissProt: %s", args.diamond_dir)
        diamond = parse_diamond_swissprot(args.diamond_dir)
        logger.info("  DIAMOND: %d proteins with SwissProt hits", len(diamond))
        for gene in genes:
            if gene.gene_id in diamond:
                dh = diamond[gene.gene_id]
                gene.swissprot_hit = dh.get("hit", "")
                gene.swissprot_evalue = dh.get("evalue", "")
                gene.swissprot_product = dh.get("product", "")

    # =========================================================================
    # 7. Merge DRAM (Track A — unique DBs)
    # =========================================================================
    if args.dram_dir:
        logger.info("Loading DRAM: %s", args.dram_dir)
        dram = parse_dram(args.dram_dir)
        logger.info("  DRAM: %d annotated genes", len(dram))
        # DRAM uses its own gene IDs → match by coordinates via Bakta
        # For now, direct gene ID matching against consensus
        for gene in genes:
            if gene.gene_id in dram:
                dr = dram[gene.gene_id]
                gene.dram_kegg = dr.get("kegg", "")
                gene.dram_pfam = dr.get("pfam", "")
                gene.dram_cazy = dr.get("cazy", "")
                gene.dram_merops = dr.get("merops", "")
                gene.dram_vogdb = dr.get("vogdb", "")

    # =========================================================================
    # 8. Merge Phanotate prophage genes (Track C)
    # =========================================================================
    if args.phanotate_gff:
        logger.info("Loading Phanotate: %s", args.phanotate_gff)
        phano_genes = parse_phanotate_gff(args.phanotate_gff)
        logger.info("  Phanotate: %d prophage CDS features", len(phano_genes))
        phano_count = 0
        for pg in phano_genes:
            phano_count += 1
            gene = AnnotatedGene(
                gene_id=f"prophage_{phano_count:05d}",
                seqid=pg["seqid"],
                start=pg["start"],
                end=pg["end"],
                strand=pg["strand"],
                merge_source="phanotate",
                region_type="prophage",
                callers=["phanotate"],
            )
            genes.append(gene)
        logger.info("  Added %d prophage genes", phano_count)

    # =========================================================================
    # 9. Sort by genomic position
    # =========================================================================
    genes.sort(key=lambda g: (g.seqid, g.start))

    # =========================================================================
    # 10. Write outputs
    # =========================================================================
    gff_out = os.path.join(args.output_dir, "final_annotation.gff3")
    tsv_out = os.path.join(args.output_dir, "annotation_matrix.tsv")
    gbk_out = os.path.join(args.output_dir, "final_annotation.gbk")
    json_out = os.path.join(args.output_dir, "summary.json")

    logger.info("Writing GFF3: %s", gff_out)
    write_gff3(genes, gff_out)

    logger.info("Writing TSV: %s", tsv_out)
    write_tsv(genes, tsv_out)

    logger.info("Writing GenBank: %s", gbk_out)
    write_genbank(genes, args.assembly, gbk_out, args.sample_name)

    logger.info("Writing summary JSON: %s", json_out)
    write_summary_json(genes, args.sample_name, json_out)

    # =========================================================================
    # Summary
    # =========================================================================
    logger.info("=== RECONCILIATION SUMMARY ===")
    logger.info("Total genes in final annotation: %d", len(genes))
    source_counts = defaultdict(int)
    for g in genes:
        source_counts[g.merge_source] += 1
    for src, count in sorted(source_counts.items()):
        logger.info("  %s: %d", src, count)

    annotated = sum(1 for g in genes
                    if resolve_product(g)[0] != "hypothetical protein")
    logger.info("Annotated (non-hypothetical): %d/%d (%.1f%%)",
                annotated, len(genes),
                100 * annotated / len(genes) if genes else 0)


if __name__ == "__main__":
    main()
