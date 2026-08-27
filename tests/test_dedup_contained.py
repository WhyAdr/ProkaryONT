from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "dedup_contained.py"


def paf_row(
    query: str,
    query_length: int,
    query_start: int,
    query_end: int,
    strand: str,
    target: str,
    target_length: int,
    target_start: int,
    target_end: int,
    matches: int,
    block_length: int,
) -> str:
    return "\t".join(
        map(
            str,
            (
                query,
                query_length,
                query_start,
                query_end,
                strand,
                target,
                target_length,
                target_start,
                target_end,
                matches,
                block_length,
                60,
            ),
        )
    )


class DedupContainedTests(unittest.TestCase):
    def run_case(
        self,
        fasta: str,
        paf: str = "",
        *extra_args: str,
    ) -> tuple[subprocess.CompletedProcess[str], list[dict[str, str]]]:
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary = Path(temporary_directory)
            fasta_path = temporary / "input.fasta"
            paf_path = temporary / "input.paf"
            report_path = temporary / "events.tsv"
            fasta_path.write_text(fasta, encoding="utf-8")
            paf_path.write_text(paf, encoding="utf-8")
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    str(fasta_path),
                    str(paf_path),
                    "--report-tsv",
                    str(report_path),
                    *extra_args,
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            rows: list[dict[str, str]] = []
            if report_path.exists():
                with report_path.open(encoding="utf-8", newline="") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
            return result, rows

    def test_contained_shorter_query_is_dropped(self) -> None:
        fasta = ">long\n" + "A" * 120 + "\n>short\n" + "A" * 100 + "\n"
        paf = paf_row("short", 100, 0, 100, "+", "long", 120, 0, 100, 100, 100)
        result, rows = self.run_case(fasta, paf)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, ">long\n" + "A" * 120 + "\n")
        self.assertEqual(rows[1]["reason"], "contained")
        self.assertEqual(rows[1]["target_id"], "long")

    def test_noncontained_query_is_retained(self) -> None:
        fasta = ">long\n" + "A" * 120 + "\n>short\n" + "C" * 100 + "\n"
        paf = paf_row("short", 100, 0, 80, "+", "long", 120, 0, 80, 80, 80)
        result, rows = self.run_case(fasta, paf)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(">short\n", result.stdout)
        self.assertEqual(rows[1]["action"], "keep")

    def test_reverse_strand_containment_is_reported(self) -> None:
        fasta = ">long\n" + "A" * 120 + "\n>short\n" + "C" * 100 + "\n"
        paf = paf_row("short", 100, 0, 100, "-", "long", 120, 20, 120, 100, 100)
        result, rows = self.run_case(fasta, paf)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows[1]["action"], "drop")
        self.assertEqual(rows[1]["strand"], "-")

    def test_equal_length_near_duplicates_are_kept(self) -> None:
        first = "A" * 99 + "C"
        second = "A" * 100
        fasta = f">a\n{first}\n>b\n{second}\n"
        paf = paf_row("b", 100, 0, 100, "+", "a", 100, 0, 100, 99, 100)
        result, rows = self.run_case(fasta, paf)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([row["action"] for row in rows], ["keep", "keep"])

    def test_equal_length_exact_duplicates_use_lexical_tie_break(self) -> None:
        fasta = ">zeta\nACGTTGCA\n>alpha\nACGTTGCA\n"
        result, rows = self.run_case(fasta)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, ">alpha\nACGTTGCA\n")
        zeta = next(row for row in rows if row["query_id"] == "zeta")
        self.assertEqual(zeta["reason"], "exact_duplicate")
        self.assertEqual(zeta["target_id"], "alpha")

    def test_reverse_complement_exact_duplicates_are_collapsed(self) -> None:
        fasta = ">b\nACGTTT\n>a\nAAACGT\n"
        result, rows = self.run_case(fasta)
        self.assertEqual(result.returncode, 0, result.stderr)
        b_row = next(row for row in rows if row["query_id"] == "b")
        self.assertEqual(b_row["action"], "drop")
        self.assertEqual(b_row["strand"], "-")

    def test_split_alignments_to_one_target_can_satisfy_coverage(self) -> None:
        fasta = ">target\n" + "A" * 150 + "\n>query\n" + "C" * 100 + "\n"
        paf = "\n".join(
            [
                paf_row("query", 100, 0, 50, "+", "target", 150, 0, 50, 50, 50),
                paf_row("query", 100, 50, 100, "+", "target", 150, 100, 150, 50, 50),
            ]
        )
        result, rows = self.run_case(fasta, paf)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows[1]["action"], "drop")
        self.assertEqual(rows[1]["covered_bp"], "100")

    def test_coverage_is_not_pooled_across_targets(self) -> None:
        fasta = (
            ">target_a\n"
            + "A" * 150
            + "\n>target_b\n"
            + "G" * 140
            + "\n>query\n"
            + "C" * 100
            + "\n"
        )
        paf = "\n".join(
            [
                paf_row("query", 100, 0, 50, "+", "target_a", 150, 0, 50, 50, 50),
                paf_row("query", 100, 50, 100, "+", "target_b", 140, 50, 100, 50, 50),
            ]
        )
        result, rows = self.run_case(fasta, paf)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows[2]["action"], "keep")

    def test_overlapping_intervals_are_not_double_counted(self) -> None:
        fasta = ">target\n" + "A" * 150 + "\n>query\n" + "C" * 100 + "\n"
        paf = "\n".join(
            [
                paf_row("query", 100, 0, 60, "+", "target", 150, 0, 60, 60, 60),
                paf_row("query", 100, 40, 80, "+", "target", 150, 80, 120, 40, 40),
            ]
        )
        result, rows = self.run_case(fasta, paf, "--min-cov", "0.9")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows[1]["action"], "keep")

    def test_threshold_boundary_is_inclusive(self) -> None:
        fasta = ">target\n" + "A" * 120 + "\n>query\n" + "C" * 100 + "\n"
        paf = paf_row("query", 100, 0, 99, "+", "target", 120, 0, 99, 98, 99)
        result, rows = self.run_case(
            fasta, paf, "--min-cov", "0.99", "--min-identity", str(98 / 99)
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows[1]["action"], "drop")

    def test_length_filter_is_disabled_by_default_and_explicit_when_positive(
        self,
    ) -> None:
        fasta = ">long\nAAAA\n>tiny\nC\n"
        default_result, default_rows = self.run_case(fasta)
        self.assertEqual(default_result.returncode, 0, default_result.stderr)
        self.assertEqual([row["action"] for row in default_rows], ["keep", "keep"])
        filtered_result, filtered_rows = self.run_case(fasta, "", "--min-len", "2")
        self.assertEqual(filtered_result.returncode, 0, filtered_result.stderr)
        self.assertEqual(filtered_rows[1]["reason"], "short")

    def test_duplicate_fasta_ids_fail_without_report(self) -> None:
        result, rows = self.run_case(">dup\nAAAA\n>dup extra\nCCCC\n")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("duplicate FASTA identifier", result.stderr)
        self.assertEqual(rows, [])

    def test_malformed_paf_and_stale_lengths_fail(self) -> None:
        fasta = ">a\nAAAA\n>b\nCCCCCC\n"
        malformed, _ = self.run_case(fasta, "a\t4\t0")
        self.assertNotEqual(malformed.returncode, 0)
        stale = paf_row("a", 5, 0, 4, "+", "b", 6, 0, 4, 4, 4)
        stale_result, _ = self.run_case(fasta, stale)
        self.assertNotEqual(stale_result.returncode, 0)
        self.assertIn("FASTA length", stale_result.stderr)

    def test_unknown_paf_identifier_and_invalid_threshold_fail(self) -> None:
        fasta = ">a\nAAAA\n>b\nCCCCCC\n"
        unknown = paf_row("missing", 4, 0, 4, "+", "b", 6, 0, 4, 4, 4)
        unknown_result, _ = self.run_case(fasta, unknown)
        self.assertNotEqual(unknown_result.returncode, 0)
        invalid_result, _ = self.run_case(fasta, "", "--min-cov", "0")
        self.assertNotEqual(invalid_result.returncode, 0)

    def test_report_and_output_counts_reconcile_and_order_is_preserved(self) -> None:
        fasta = ">first note\nAAAA\n>second note\nCCCC\n>third note\nACAC\n"
        result, rows = self.run_case(fasta)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            result.stdout,
            ">first note\nAAAA\n>second note\nCCCC\n>third note\nACAC\n",
        )
        kept = sum(row["action"] == "keep" for row in rows)
        self.assertEqual(kept, result.stdout.count(">"))


if __name__ == "__main__":
    unittest.main()
