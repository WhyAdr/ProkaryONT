from __future__ import annotations

import csv
import hashlib
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAGE3 = ROOT / "03_autocycler_assemble.sh"
DEDUP = ROOT / "dedup_contained.py"


def bash_path() -> str | None:
    configured = os.environ.get("BASH")
    if configured and Path(configured).exists():
        return configured
    discovered = shutil.which("bash")
    if discovered and "WindowsApps" not in discovered:
        return discovered
    git_bash = Path(r"C:\Program Files\Git\bin\bash.exe")
    return str(git_bash) if git_bash.exists() else None


def posix(path: Path) -> str:
    rendered = path.resolve().as_posix()
    if len(rendered) >= 3 and rendered[1:3] == ":/":
        return f"/{rendered[0].lower()}{rendered[2:]}"
    return rendered


def sequence_for(index: int, length: int = 110) -> str:
    digits = []
    value = index + 1
    alphabet = "ACGT"
    while value:
        digits.append(alphabet[value % 4])
        value //= 4
    suffix = "".join(digits) or "A"
    return ("AC" * length + suffix)[-length:]


def fragmented_fasta(count: int) -> str:
    if count < 2:
        raise ValueError("fixture needs at least two records")
    records = [("target", "A" * 120), ("contained", "C" * 100)]
    records.extend(
        (f"unique_{index:02d}", sequence_for(index, 121 + index))
        for index in range(count - 2)
    )
    return "".join(f">{identifier}\n{sequence}\n" for identifier, sequence in records)


class Stage3ContigPolicyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.bash = bash_path()
        if cls.bash is None:
            raise unittest.SkipTest("Bash is unavailable")
        source = STAGE3.read_text(encoding="utf-8")
        start = source.index("# --- Contig assessment/preparation helpers")
        end = source.index(
            "# ==============================================================================\n# STEP 9",
            start,
        )
        cls.function_source = source[start:end]

    def setUp(self) -> None:
        temporary_root = ROOT / ".test-tmp"
        temporary_root.mkdir(parents=True, exist_ok=True)
        self.temporary_directory = tempfile.TemporaryDirectory(dir=temporary_root)
        self.work = Path(self.temporary_directory.name)
        self.assemblies = self.work / "assemblies"
        self.prepared = self.work / "assemblies_prepared"
        self.assessment = self.work / "assembly_assessment"
        self.mock_bin = self.work / "mock_bin"
        self.marker = self.work / "minimap2.invocations"
        self.assemblies.mkdir()
        self.mock_bin.mkdir()
        mock_minimap2 = self.mock_bin / "minimap2"
        mock_minimap2.write_text(
            """#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == "--version" ]]; then
    echo "mock-minimap2 1.0"
    exit 0
fi
printf 'invoked\\n' >> "${MOCK_MINIMAP_MARKER}"
printf 'contained\\t100\\t0\\t100\\t+\\ttarget\\t120\\t0\\t100\\t100\\t100\\t60\\n'
""",
            encoding="utf-8",
            newline="\n",
        )
        mock_minimap2.chmod(
            mock_minimap2.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH
        )
        python3_wrapper = self.mock_bin / "python3"
        python3_wrapper.write_text(
            f'#!/usr/bin/env bash\nexec {posix(Path(sys.executable))!r} "$@"\n',
            encoding="utf-8",
            newline="\n",
        )
        python3_wrapper.chmod(
            python3_wrapper.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH
        )

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def run_policy(
        self,
        *,
        trigger: int = 25,
        skip: bool = False,
        max_contigs: str = "",
        resume_target: str = "",
    ) -> subprocess.CompletedProcess[str]:
        harness = self.work / "harness.sh"
        harness.write_text(
            f"""#!/usr/bin/env bash
set -euo pipefail
log_info() {{ printf 'INFO: %s\\n' "$*"; }}
log_warn() {{ printf 'WARN: %s\\n' "$*" >&2; }}
log_error() {{ printf 'ERROR: %s\\n' "$*" >&2; exit 1; }}
assemblies_dir={posix(self.assemblies)!r}
prepared_assemblies_dir={posix(self.prepared)!r}
assembly_assessment_dir={posix(self.assessment)!r}
dedup_script={posix(DEDUP)!r}
dedup_trigger_contigs={trigger}
dedup_min_cov=0.99
dedup_min_identity=0.99
dedup_min_len=0
dedup_threads=4
skip_contained_dedup={"true" if skip else "false"}
autocycler_max_contigs={max_contigs!r}
dry_run=
_resume_target={resume_target!r}
export PATH={posix(self.mock_bin)!r}:"$PATH"
export MOCK_MINIMAP_MARKER={posix(self.marker)!r}
{self.function_source}
_prepare_autocycler_assemblies
printf 'decision=%s\\neffective=%s\\nerror=%s\\n' \\
    "${{contig_policy_decision}}" "${{effective_max_contigs}}" "${{contig_policy_error:-}}"
""",
            encoding="utf-8",
            newline="\n",
        )
        return subprocess.run(
            [self.bash, posix(harness)],
            cwd=self.work,
            check=False,
            capture_output=True,
            text=True,
            timeout=20,
        )

    def report_rows(self) -> list[dict[str, str]]:
        with (self.assessment / "assemblies.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_below_trigger_is_copied_without_minimap2(self) -> None:
        source = self.assemblies / "raven_01.fasta"
        source.write_text(">one\nACGT\n>two\nAACG\n", encoding="utf-8")
        source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
        result = self.run_policy()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(self.marker.exists())
        self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), source_hash)
        self.assertEqual(
            (self.prepared / source.name).read_bytes(), source.read_bytes()
        )
        self.assertEqual(self.report_rows()[0]["status"], "not_fragmented")

    def test_only_fragmented_file_is_cleaned_and_successful_paf_is_deleted(
        self,
    ) -> None:
        fragmented = self.assemblies / "flye_01.fasta"
        fragmented.write_text(fragmented_fasta(26), encoding="utf-8")
        normal = self.assemblies / "raven_01.fasta"
        normal.write_text(">one\nACGT\n>two\nAACG\n", encoding="utf-8")
        source_bytes = fragmented.read_bytes()
        result = self.run_policy()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.marker.read_text(encoding="utf-8").count("invoked"), 1)
        rows = {row["assembly"]: row for row in self.report_rows()}
        self.assertEqual(rows["flye_01.fasta"]["raw_contigs"], "26")
        self.assertEqual(rows["flye_01.fasta"]["post_contigs"], "25")
        self.assertEqual(rows["raven_01.fasta"]["status"], "not_fragmented")
        self.assertEqual(fragmented.read_bytes(), source_bytes)
        prepared_text = (self.prepared / "flye_01.fasta").read_text(encoding="utf-8")
        self.assertNotIn(">contained\n", prepared_text)
        self.assertIn("Autocycler_consensus_weight=2", prepared_text)
        self.assertFalse(any((self.assessment / "dedup").glob("*.paf")))
        self.assertTrue(
            (self.assessment / "dedup" / "flye_01.fasta.events.tsv").is_file()
        )
        self.assertTrue(
            (self.assessment / "dedup" / "flye_01.fasta.minimap2.log").is_file()
        )
        self.assertTrue(
            (self.assessment / "dedup" / "flye_01.fasta.dedup.log").is_file()
        )

    def test_skip_contained_dedup_still_assesses_without_invoking_minimap2(
        self,
    ) -> None:
        source = self.assemblies / "flye_01.fasta"
        source.write_text(fragmented_fasta(26), encoding="utf-8")
        result = self.run_policy(skip=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(self.marker.exists())
        row = self.report_rows()[0]
        self.assertEqual(row["status"], "fragmented_skipped")
        self.assertEqual(row["post_contigs"], "26")

    def test_cleanup_can_bring_mean_to_default_guard(self) -> None:
        (self.assemblies / "flye_01.fasta").write_text(
            fragmented_fasta(26), encoding="utf-8"
        )
        result = self.run_policy()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("decision=default_25_accepted", result.stdout)
        self.assertIn("effective=25", result.stdout)

    def test_post_mean_above_default_stops_with_exact_override(self) -> None:
        (self.assemblies / "flye_01.fasta").write_text(
            fragmented_fasta(27), encoding="utf-8"
        )
        result = self.run_policy()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("decision=stop_requires_override", result.stdout)
        self.assertIn("--max-contigs 26 or higher", result.stdout)

    def test_sufficient_and_insufficient_explicit_overrides(self) -> None:
        (self.assemblies / "flye_01.fasta").write_text(
            fragmented_fasta(27), encoding="utf-8"
        )
        sufficient = self.run_policy(max_contigs="26")
        self.assertEqual(sufficient.returncode, 0, sufficient.stderr)
        self.assertIn("decision=explicit_override_accepted", sufficient.stdout)
        self.assertIn("effective=26", sufficient.stdout)

        insufficient = self.run_policy(max_contigs="25")
        self.assertEqual(insufficient.returncode, 0, insufficient.stderr)
        self.assertIn("decision=stop_override_insufficient", insufficient.stdout)
        self.assertIn("requires at least 26", insufficient.stdout)

    def test_matching_manifest_reuses_prepared_cohort(self) -> None:
        (self.assemblies / "flye_01.fasta").write_text(
            fragmented_fasta(26), encoding="utf-8"
        )
        first = self.run_policy()
        self.assertEqual(first.returncode, 0, first.stderr)
        prepared_hash = hashlib.sha256(
            (self.prepared / "flye_01.fasta").read_bytes()
        ).hexdigest()
        second = self.run_policy()
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertIn("reusing validated assemblies_prepared", second.stdout)
        self.assertEqual(self.marker.read_text(encoding="utf-8").count("invoked"), 1)
        self.assertEqual(
            hashlib.sha256((self.prepared / "flye_01.fasta").read_bytes()).hexdigest(),
            prepared_hash,
        )

    def test_empty_output_is_excluded_and_reported(self) -> None:
        valid = self.assemblies / "raven_01.fasta"
        valid.write_text(">one\nACGT\n", encoding="utf-8")
        (self.assemblies / "failed_01.fasta").write_text("", encoding="utf-8")
        result = self.run_policy()
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = {row["assembly"]: row for row in self.report_rows()}
        self.assertEqual(rows["failed_01.fasta"]["status"], "excluded_empty")
        self.assertTrue((self.prepared / "raven_01.fasta").is_file())
        self.assertFalse((self.prepared / "failed_01.fasta").exists())

    def test_missing_expected_output_is_excluded_and_reported(self) -> None:
        valid = self.assemblies / "raven_01.fasta"
        valid.write_text(">one\nACGT\n", encoding="utf-8")
        (self.assemblies / "jobs.txt").write_text(
            f"autocycler helper flye --out_prefix {self.assemblies / 'flye_01'}\n",
            encoding="utf-8",
        )
        result = self.run_policy()
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = {row["assembly"]: row for row in self.report_rows()}
        self.assertEqual(rows["flye_01.fasta"]["status"], "excluded_missing")
        self.assertFalse((self.prepared / "flye_01.fasta").exists())

    def test_stale_manifest_rebuilds_prepared_cohort(self) -> None:
        source = self.assemblies / "flye_01.fasta"
        source.write_text(fragmented_fasta(26), encoding="utf-8")
        first = self.run_policy()
        self.assertEqual(first.returncode, 0, first.stderr)
        source.write_text(
            fragmented_fasta(26) + ">new_record\nACGACGACG\n", encoding="utf-8"
        )
        second = self.run_policy(max_contigs="26")
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertNotIn("reusing validated assemblies_prepared", second.stdout)
        self.assertEqual(self.marker.read_text(encoding="utf-8").count("invoked"), 2)
        self.assertIn(">new_record ", (self.prepared / source.name).read_text())

    def test_modified_prepared_fasta_is_not_reused(self) -> None:
        source = self.assemblies / "raven_01.fasta"
        source.write_text(">one\nACGT\n>two\nAACG\n", encoding="utf-8")
        first = self.run_policy()
        self.assertEqual(first.returncode, 0, first.stderr)
        prepared = self.prepared / source.name
        prepared.write_text(">one\nTTTT\n>two\nAACG\n", encoding="utf-8")
        second = self.run_policy()
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertNotIn("reusing validated assemblies_prepared", second.stdout)
        self.assertEqual(prepared.read_bytes(), source.read_bytes())

    def test_late_resume_does_not_rebuild_changed_inputs(self) -> None:
        source = self.assemblies / "flye_01.fasta"
        source.write_text(fragmented_fasta(26), encoding="utf-8")
        first = self.run_policy()
        self.assertEqual(first.returncode, 0, first.stderr)
        prepared_hash = hashlib.sha256(
            (self.prepared / source.name).read_bytes()
        ).hexdigest()
        source.write_text(
            fragmented_fasta(26) + ">changed_after_graph\nACGACGACG\n",
            encoding="utf-8",
        )
        resumed = self.run_policy(resume_target="trim")
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        self.assertIn("will not be changed underneath trim/resolve", resumed.stderr)
        self.assertEqual(self.marker.read_text(encoding="utf-8").count("invoked"), 1)
        self.assertEqual(
            hashlib.sha256((self.prepared / source.name).read_bytes()).hexdigest(),
            prepared_hash,
        )

    def test_stage3_passes_same_explicit_guard_to_both_commands(self) -> None:
        source = STAGE3.read_text(encoding="utf-8")
        self.assertIn('autocycler compress -i "${autocycler_input_dir}"', source)
        self.assertGreaterEqual(
            source.count('--max_contigs "${effective_max_contigs}"'), 2
        )


if __name__ == "__main__":
    unittest.main()
