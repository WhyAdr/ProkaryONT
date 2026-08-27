#!/usr/bin/env python3
"""Remove conservatively contained contigs using a validated self-PAF.

Retained FASTA records are written to stdout. A mandatory per-record TSV report
records every keep/drop decision. Equal-length records are collapsed only when
their complete sequences are identical, including reverse-complement identity.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
REPORT_FIELDS = (
    "query_id",
    "query_length",
    "action",
    "reason",
    "target_id",
    "target_length",
    "strand",
    "covered_bp",
    "coverage",
    "identity",
    "min_cov",
    "min_identity",
    "min_len",
)


class InputError(ValueError):
    """Raised when FASTA, PAF, or arguments violate the input contract."""


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    header: str
    sequence: str


@dataclass(frozen=True)
class Alignment:
    query_start: int
    query_end: int
    matches: int
    block_length: int


@dataclass(frozen=True)
class Decision:
    action: str
    reason: str
    target_id: str = ""
    target_length: int | None = None
    strand: str = ""
    covered_bp: int = 0
    coverage: float = 0.0
    identity: float = 0.0


def fraction(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number in (0, 1]") from exc
    if not 0.0 < parsed <= 1.0:
        raise argparse.ArgumentTypeError("must be in (0, 1]")
    return parsed


def nonnegative_integer(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a non-negative integer") from exc
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be a non-negative integer")
    return parsed


def parse_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    seen_ids: set[str] = set()
    header: str | None = None
    identifier: str | None = None
    sequence_parts: list[str] = []

    def finish_record(line_number: int) -> None:
        nonlocal header, identifier, sequence_parts
        if header is None or identifier is None:
            return
        sequence = "".join(sequence_parts)
        if not sequence:
            raise InputError(
                f"{path}: FASTA record {identifier!r} has an empty sequence "
                f"(detected near line {line_number})"
            )
        invalid = sorted(set(sequence.upper()) - DNA_ALPHABET)
        if invalid:
            shown = "".join(invalid[:10])
            raise InputError(
                f"{path}: FASTA record {identifier!r} contains non-IUPAC DNA "
                f"character(s): {shown!r}"
            )
        records.append(FastaRecord(identifier, header, sequence))

    try:
        handle = path.open("r", encoding="utf-8", newline=None)
    except OSError as exc:
        raise InputError(f"cannot open FASTA {path}: {exc}") from exc

    with handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                raise InputError(
                    f"{path}:{line_number}: blank lines are not allowed in FASTA"
                )
            if line.startswith(">"):
                finish_record(line_number)
                header = line
                header_text = line[1:].strip()
                if not header_text:
                    raise InputError(
                        f"{path}:{line_number}: FASTA header has no identifier"
                    )
                identifier = header_text.split()[0]
                if identifier in seen_ids:
                    raise InputError(
                        f"{path}:{line_number}: duplicate FASTA identifier {identifier!r}"
                    )
                seen_ids.add(identifier)
                sequence_parts = []
                continue
            if header is None:
                raise InputError(
                    f"{path}:{line_number}: sequence data occurs before the first FASTA header"
                )
            if any(character.isspace() for character in line):
                raise InputError(
                    f"{path}:{line_number}: whitespace is not allowed inside a FASTA sequence"
                )
            sequence_parts.append(line)

    finish_record(line_number + 1 if "line_number" in locals() else 1)
    if not records:
        raise InputError(f"{path}: no FASTA records found")
    return records


def parse_paf(
    path: Path, records_by_id: dict[str, FastaRecord]
) -> dict[tuple[str, str, str], list[Alignment]]:
    groups: dict[tuple[str, str, str], list[Alignment]] = defaultdict(list)
    try:
        handle = path.open("r", encoding="utf-8", newline=None)
    except OSError as exc:
        raise InputError(f"cannot open PAF {path}: {exc}") from exc

    with handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                raise InputError(
                    f"{path}:{line_number}: PAF row has {len(fields)} fields; at least 12 required"
                )
            query_id, target_id, strand = fields[0], fields[5], fields[4]
            if query_id not in records_by_id:
                raise InputError(
                    f"{path}:{line_number}: unknown PAF query identifier {query_id!r}"
                )
            if target_id not in records_by_id:
                raise InputError(
                    f"{path}:{line_number}: unknown PAF target identifier {target_id!r}"
                )
            if strand not in {"+", "-"}:
                raise InputError(
                    f"{path}:{line_number}: PAF strand must be '+' or '-', got {strand!r}"
                )
            try:
                query_length = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                target_length = int(fields[6])
                target_start = int(fields[7])
                target_end = int(fields[8])
                matches = int(fields[9])
                block_length = int(fields[10])
                mapping_quality = int(fields[11])
            except ValueError as exc:
                raise InputError(
                    f"{path}:{line_number}: PAF mandatory numeric field is not an integer"
                ) from exc

            expected_query_length = len(records_by_id[query_id].sequence)
            expected_target_length = len(records_by_id[target_id].sequence)
            if query_length != expected_query_length:
                raise InputError(
                    f"{path}:{line_number}: query length for {query_id!r} is {query_length}, "
                    f"but FASTA length is {expected_query_length}"
                )
            if target_length != expected_target_length:
                raise InputError(
                    f"{path}:{line_number}: target length for {target_id!r} is {target_length}, "
                    f"but FASTA length is {expected_target_length}"
                )
            if not 0 <= query_start < query_end <= query_length:
                raise InputError(f"{path}:{line_number}: invalid PAF query coordinates")
            if not 0 <= target_start < target_end <= target_length:
                raise InputError(
                    f"{path}:{line_number}: invalid PAF target coordinates"
                )
            if block_length <= 0 or not 0 <= matches <= block_length:
                raise InputError(
                    f"{path}:{line_number}: invalid PAF match/block lengths"
                )
            if not 0 <= mapping_quality <= 255:
                raise InputError(
                    f"{path}:{line_number}: PAF mapping quality must be 0..255"
                )

            if query_id == target_id:
                continue
            if target_length <= query_length:
                # Equal-length handling uses exact sequence comparison below.
                continue
            groups[(query_id, target_id, strand)].append(
                Alignment(query_start, query_end, matches, block_length)
            )
    return groups


def merged_coverage(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    sorted_intervals = sorted(intervals)
    start, end = sorted_intervals[0]
    covered = 0
    for next_start, next_end in sorted_intervals[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            covered += end - start
            start, end = next_start, next_end
    return covered + (end - start)


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(COMPLEMENT)[::-1]


def exact_duplicate_target(
    record: FastaRecord, records: list[FastaRecord]
) -> tuple[FastaRecord, str] | None:
    query = record.sequence.upper()
    reverse_query = reverse_complement(record.sequence)
    candidates: list[tuple[str, FastaRecord, str]] = []
    for target in records:
        if target.identifier >= record.identifier or len(target.sequence) != len(
            record.sequence
        ):
            continue
        target_sequence = target.sequence.upper()
        if query == target_sequence:
            candidates.append((target.identifier, target, "+"))
        elif reverse_query == target_sequence:
            candidates.append((target.identifier, target, "-"))
    if not candidates:
        return None
    _, target, strand = min(candidates, key=lambda item: item[0])
    return target, strand


def decide(
    records: list[FastaRecord],
    groups: dict[tuple[str, str, str], list[Alignment]],
    min_cov: float,
    min_identity: float,
    min_len: int,
) -> dict[str, Decision]:
    by_id = {record.identifier: record for record in records}
    decisions: dict[str, Decision] = {}

    for record in records:
        length = len(record.sequence)
        if min_len > 0 and length < min_len:
            decisions[record.identifier] = Decision(action="drop", reason="short")
            continue

        duplicate = exact_duplicate_target(record, records)
        if duplicate is not None:
            target, strand = duplicate
            decisions[record.identifier] = Decision(
                action="drop",
                reason="exact_duplicate",
                target_id=target.identifier,
                target_length=len(target.sequence),
                strand=strand,
                covered_bp=length,
                coverage=1.0,
                identity=1.0,
            )
            continue

        qualifying: list[tuple[float, float, int, str, str, int]] = []
        for (query_id, target_id, strand), alignments in groups.items():
            if query_id != record.identifier:
                continue
            covered_bp = merged_coverage(
                [
                    (alignment.query_start, alignment.query_end)
                    for alignment in alignments
                ]
            )
            coverage = covered_bp / length
            total_block_length = sum(alignment.block_length for alignment in alignments)
            identity = (
                sum(alignment.matches for alignment in alignments) / total_block_length
            )
            if coverage >= min_cov and identity >= min_identity:
                target_length = len(by_id[target_id].sequence)
                qualifying.append(
                    (coverage, identity, target_length, target_id, strand, covered_bp)
                )

        if qualifying:
            coverage, identity, target_length, target_id, strand, covered_bp = min(
                qualifying,
                key=lambda item: (-item[0], -item[1], -item[2], item[3], item[4]),
            )
            decisions[record.identifier] = Decision(
                action="drop",
                reason="contained",
                target_id=target_id,
                target_length=target_length,
                strand=strand,
                covered_bp=covered_bp,
                coverage=coverage,
                identity=identity,
            )
        else:
            decisions[record.identifier] = Decision(
                action="keep", reason="not_contained"
            )
    return decisions


def report_row(
    record: FastaRecord,
    decision: Decision,
    min_cov: float,
    min_identity: float,
    min_len: int,
) -> dict[str, str | int]:
    return {
        "query_id": record.identifier,
        "query_length": len(record.sequence),
        "action": decision.action,
        "reason": decision.reason,
        "target_id": decision.target_id,
        "target_length": ""
        if decision.target_length is None
        else decision.target_length,
        "strand": decision.strand,
        "covered_bp": decision.covered_bp,
        "coverage": f"{decision.coverage:.6f}",
        "identity": f"{decision.identity:.6f}",
        "min_cov": f"{min_cov:.6f}",
        "min_identity": f"{min_identity:.6f}",
        "min_len": min_len,
    }


def write_report_atomic(
    path: Path,
    records: list[FastaRecord],
    decisions: dict[str, Decision],
    min_cov: float,
    min_identity: float,
    min_len: int,
) -> None:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            newline="",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=REPORT_FIELDS, delimiter="\t")
            writer.writeheader()
            for record in records:
                writer.writerow(
                    report_row(
                        record,
                        decisions[record.identifier],
                        min_cov,
                        min_identity,
                        min_len,
                    )
                )
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    except OSError as exc:
        try:
            temporary_path.unlink(missing_ok=True)
        except (OSError, UnboundLocalError):
            pass
        raise InputError(f"cannot write report {path}: {exc}") from exc


def write_retained_fasta(
    handle: TextIO, records: list[FastaRecord], decisions: dict[str, Decision]
) -> int:
    retained = 0
    for record in records:
        if decisions[record.identifier].action == "drop":
            continue
        handle.write(f"{record.header}\n{record.sequence}\n")
        retained += 1
    return retained


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Conservatively remove contained contigs using a self-alignment PAF."
    )
    parser.add_argument("fasta", type=Path, help="input FASTA")
    parser.add_argument("paf", type=Path, help="minimap2 self-alignment PAF")
    parser.add_argument(
        "--report-tsv", type=Path, required=True, help="per-record decision report"
    )
    parser.add_argument(
        "--min-cov", type=fraction, default=0.99, help="minimum single-target coverage"
    )
    parser.add_argument(
        "--min-identity",
        type=fraction,
        default=0.99,
        help="minimum aggregate alignment identity",
    )
    parser.add_argument(
        "--min-len",
        type=nonnegative_integer,
        default=0,
        help="length-only removal threshold; 0 disables it",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        resolved_inputs = {args.fasta.resolve(), args.paf.resolve()}
        if args.report_tsv.resolve() in resolved_inputs:
            raise InputError("--report-tsv must not overwrite the input FASTA or PAF")
        records = parse_fasta(args.fasta)
        records_by_id = {record.identifier: record for record in records}
        groups = parse_paf(args.paf, records_by_id)
        decisions = decide(
            records,
            groups,
            min_cov=args.min_cov,
            min_identity=args.min_identity,
            min_len=args.min_len,
        )
        retained = sum(
            decisions[record.identifier].action == "keep" for record in records
        )
        if retained == 0:
            raise InputError("containment policy would remove every FASTA record")
        write_report_atomic(
            args.report_tsv,
            records,
            decisions,
            min_cov=args.min_cov,
            min_identity=args.min_identity,
            min_len=args.min_len,
        )
        written = write_retained_fasta(sys.stdout, records, decisions)
    except (InputError, OSError, UnicodeError) as exc:
        parser.exit(2, f"error: {exc}\n")

    dropped_contained = sum(
        decision.reason in {"contained", "exact_duplicate"}
        for decision in decisions.values()
    )
    dropped_short = sum(decision.reason == "short" for decision in decisions.values())
    sys.stderr.write(
        f"[{args.fasta}] Total: {len(records)} | Retained: {written} | "
        f"Dropped: {len(records) - written} "
        f"({dropped_contained} contained/exact, {dropped_short} short) | "
        f"Report: {args.report_tsv}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
