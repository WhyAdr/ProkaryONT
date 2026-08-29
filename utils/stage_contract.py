#!/usr/bin/env python3
"""Small, dependency-free helpers for ProkaryONT stage contracts.

The numbered Bash stages use this module for portable path/identity handling,
stable signatures, and strict FASTQ validation.  Keeping these operations in
one place avoids platform-specific ``stat``/``realpath`` branches in Bash.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import sys
from pathlib import Path
from typing import BinaryIO


def _reject_control_characters(value: str, label: str) -> None:
    if any(ord(char) < 32 or ord(char) == 127 for char in value):
        raise ValueError(f"{label} contains a control character")


def canonical_path(raw_path: str) -> str:
    _reject_control_characters(raw_path, "path")
    return str(Path(raw_path).expanduser().resolve(strict=False))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def input_identity(raw_path: str, include_sha256: bool) -> str:
    path = Path(raw_path).expanduser().resolve(strict=True)
    if not path.is_file():
        raise ValueError(f"not a regular file: {path}")
    stat_result = path.stat()
    digest = sha256_file(path) if include_sha256 else "metadata-only"
    fields = (str(path), str(stat_result.st_size), str(stat_result.st_mtime_ns), digest)
    for field in fields:
        _reject_control_characters(field, "identity field")
    return "\t".join(fields)


def signature_from_stdin() -> str:
    digest = hashlib.sha256()
    for chunk in iter(lambda: sys.stdin.buffer.read(1024 * 1024), b""):
        digest.update(chunk)
    return digest.hexdigest()


def _open_fastq(path: Path) -> BinaryIO:
    with path.open("rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rb")
    return path.open("rb")


def validate_fastq(raw_path: str) -> tuple[int, int]:
    path = Path(raw_path)
    if not path.is_file():
        raise ValueError(f"FASTQ does not exist: {path}")

    reads = 0
    bases = 0
    try:
        with _open_fastq(path) as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline()
                plus = handle.readline()
                quality = handle.readline()
                record_number = reads + 1
                if not sequence or not plus or not quality:
                    raise ValueError(
                        f"truncated FASTQ record {record_number} in {path}"
                    )
                header = header.rstrip(b"\r\n")
                sequence = sequence.rstrip(b"\r\n")
                plus = plus.rstrip(b"\r\n")
                quality = quality.rstrip(b"\r\n")
                if not header.startswith(b"@"):
                    raise ValueError(
                        f"FASTQ record {record_number} does not start with @ in {path}"
                    )
                if not plus.startswith(b"+"):
                    raise ValueError(
                        f"FASTQ record {record_number} has no + separator in {path}"
                    )
                if not sequence:
                    raise ValueError(
                        f"FASTQ record {record_number} has an empty sequence in {path}"
                    )
                if len(sequence) != len(quality):
                    raise ValueError(
                        f"FASTQ record {record_number} sequence/quality lengths differ in {path}"
                    )
                reads += 1
                bases += len(sequence)
    except (OSError, EOFError, gzip.BadGzipFile) as exc:
        raise ValueError(f"could not read FASTQ {path}: {exc}") from exc

    if reads == 0:
        raise ValueError(f"FASTQ contains zero reads: {path}")
    return reads, bases


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    canonical = subparsers.add_parser("canonical")
    canonical.add_argument("path")

    identity = subparsers.add_parser("identity")
    identity.add_argument("path")
    identity.add_argument("--sha256", action="store_true")

    subparsers.add_parser("signature")

    validate = subparsers.add_parser("validate-fastq")
    validate.add_argument("path")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "canonical":
            print(canonical_path(args.path))
        elif args.command == "identity":
            print(input_identity(args.path, args.sha256))
        elif args.command == "signature":
            print(signature_from_stdin())
        elif args.command == "validate-fastq":
            reads, bases = validate_fastq(args.path)
            print(f"{reads}\t{bases}")
        else:  # pragma: no cover - argparse enforces the command set.
            raise AssertionError(args.command)
    except (OSError, ValueError) as exc:
        print(f"stage_contract: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
