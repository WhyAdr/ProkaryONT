from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAGE1 = ROOT / "01_qc_estimate.sh"
STAGE2 = ROOT / "02_preprocess_filter.sh"
MOCK_GENERATOR = ROOT / "utils" / "generate_mock_environment.py"
CONTRACT = ROOT / "utils" / "stage_contract.py"


def bash_executable() -> str | None:
    if os.name == "nt":
        git_bash = Path(r"C:\Program Files\Git\bin\bash.exe")
        if git_bash.is_file():
            return str(git_bash)
    return shutil.which("bash")


class StageContractTests(unittest.TestCase):
    def test_fastq_validator_accepts_plain_and_gzip_and_rejects_truncation(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            content = b"@r1\nACGT\n+\nIIII\n@r2\nAA\n+\nII\n"
            plain = root / "reads.fastq"
            plain.write_bytes(content)
            compressed = root / "reads.fastq.gz"
            with gzip.open(compressed, "wb") as handle:
                handle.write(content)
            truncated = root / "truncated.fastq"
            truncated.write_bytes(b"@r1\nACGT\n+\n")

            for path in (plain, compressed):
                result = subprocess.run(
                    [sys.executable, str(CONTRACT), "validate-fastq", str(path)],
                    check=False,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout.strip(), "2\t6")

            result = subprocess.run(
                [sys.executable, str(CONTRACT), "validate-fastq", str(truncated)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("truncated FASTQ", result.stderr)

    def test_metadata_identity_is_default_and_sha256_is_opt_in(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "reads.fastq"
            path.write_text("@r1\nA\n+\nI\n", encoding="ascii")
            metadata = subprocess.check_output(
                [sys.executable, str(CONTRACT), "identity", str(path)], text=True
            ).strip()
            strict = subprocess.check_output(
                [
                    sys.executable,
                    str(CONTRACT),
                    "identity",
                    "--sha256",
                    str(path),
                ],
                text=True,
            ).strip()
            self.assertEqual(metadata.split("\t")[3], "metadata-only")
            self.assertRegex(strict.split("\t")[3], r"^[0-9a-f]{64}$")


@unittest.skipUnless(bash_executable(), "Bash is required for stage contract tests")
class StageBashIntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.work = Path(self.tempdir.name)
        subprocess.run([sys.executable, str(MOCK_GENERATOR)], cwd=self.work, check=True)
        self.environment = os.environ.copy()
        self.environment["PATH"] = (
            str(self.work / "mock_bin") + os.pathsep + self.environment.get("PATH", "")
        )
        self.environment["PROKARYONT_PYTHON"] = sys.executable

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def run_stage(
        self, script: Path, *arguments: str, expected: int | None = 0
    ) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            [bash_executable(), str(script), *arguments],
            cwd=self.work,
            env=self.environment,
            check=False,
            capture_output=True,
            text=True,
        )
        if expected is not None:
            self.assertEqual(
                result.returncode,
                expected,
                f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}",
            )
        return result

    def test_strict_dry_run_with_every_optional_stage2_branch_is_mutation_free(
        self,
    ) -> None:
        output = self.work / "dry-output"
        before = sorted(path.relative_to(self.work) for path in self.work.rglob("*"))
        result = self.run_stage(
            STAGE2,
            "--input-fastq",
            "input.fastq.gz",
            "--output-dir",
            "dry-output",
            "--enable-porechopabi-trim",
            "--enable-chopper-trim",
            "--chopper-trim-approach",
            "split-by-low-quality",
            "--chopper-trim-cutoff",
            "10",
            "--enable-fastcat-lint",
            "--keep-percent",
            "90",
            "--dry-run",
        )
        after = sorted(path.relative_to(self.work) for path in self.work.rglob("*"))
        self.assertFalse(output.exists())
        self.assertEqual(before, after)
        self.assertIn("No directories, logs, temporary files", result.stdout)
        self.assertIn("--keep_percent 90", result.stdout)

        stage1_output = self.work / "stage1-dry-output"
        before = sorted(path.relative_to(self.work) for path in self.work.rglob("*"))
        self.run_stage(
            STAGE1,
            "--input-fastq",
            "input.fastq.gz",
            "--output-dir",
            "stage1-dry-output",
            "--memory",
            "1",
            "--dry-run",
        )
        after = sorted(path.relative_to(self.work) for path in self.work.rglob("*"))
        self.assertFalse(stage1_output.exists())
        self.assertEqual(before, after)

    def test_stage1_default_diagnostics_manifest_and_fail_closed_resume(self) -> None:
        output = self.work / "stage1-output"
        self.run_stage(
            STAGE1,
            "--input-fastq",
            "input.fastq.gz",
            "--output-dir",
            "stage1-output",
            "--memory",
            "1",
        )
        self.assertTrue((output / ".prokaryont_stage1.complete.tsv").is_file())
        self.assertTrue(
            (output / "02_genome_size" / "genome_size_estimates.tsv").is_file()
        )
        self.assertFalse((output / "02_genome_size" / "mean_genome_size.txt").exists())
        histogram = (
            output / "01_qc" / "fastcat_sdust" / "sdust_fraction_hist.tsv"
        ).read_text(encoding="utf-8")
        self.assertIn("exact_zero\t0.00", histogram)
        self.assertIn("positive_reported\t0.99", histogram)

        resumed = self.run_stage(
            STAGE1,
            "--input-fastq",
            "input.fastq.gz",
            "--output-dir",
            "stage1-output",
            "--memory",
            "1",
        )
        self.assertIn("nothing to do", resumed.stdout)

        mismatch = self.run_stage(
            STAGE1,
            "--input-fastq",
            "input.fastq.gz",
            "--output-dir",
            "stage1-output",
            "--memory",
            "1",
            "--lint-window",
            "65",
            expected=None,
        )
        self.assertNotEqual(mismatch.returncode, 0)
        self.assertIn("--force", mismatch.stderr)

    def test_stage2_atomic_output_metrics_resume_and_force(self) -> None:
        output = self.work / "stage2-output"
        self.run_stage(
            STAGE2,
            "--input-fastq",
            "input.fastq.gz",
            "--output-dir",
            "stage2-output",
        )
        filtered = output / "filtered_input.fastq.gz"
        self.assertTrue(filtered.is_file())
        subprocess.run(
            [sys.executable, str(CONTRACT), "validate-fastq", str(filtered)],
            check=True,
        )
        metrics = (output / "02_filtering_metrics.tsv").read_text(encoding="utf-8")
        self.assertIn("base_retention_pct", metrics)
        self.assertIn("probability_mean_q", metrics)
        source = STAGE2.read_text(encoding="utf-8")
        self.assertNotIn("65.0", source)

        mismatch = self.run_stage(
            STAGE2,
            "--input-fastq",
            "input.fastq.gz",
            "--output-dir",
            "stage2-output",
            "--min-qscore",
            "8",
            expected=None,
        )
        self.assertNotEqual(mismatch.returncode, 0)
        self.run_stage(
            STAGE2,
            "--input-fastq",
            "input.fastq.gz",
            "--output-dir",
            "stage2-output",
            "--min-qscore",
            "8",
            "--force",
        )

    def test_invalid_values_and_input_output_alias_fail_before_execution(self) -> None:
        invalid_cases = (
            ("--threads", "0"),
            ("--min-length", "0"),
            ("--min-qscore", "not-a-number"),
            ("--keep-percent", "0"),
            ("--keep-percent", "100"),
            ("--keep-percent", "not-a-number"),
            ("--lint-max-proportion", "1.5"),
        )
        for option, value in invalid_cases:
            with self.subTest(option=option, value=value):
                result = self.run_stage(
                    STAGE2,
                    "--input-fastq",
                    "input.fastq.gz",
                    option,
                    value,
                    expected=None,
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(option, result.stderr)

        missing = self.run_stage(
            STAGE2,
            "--input-fastq",
            "input.fastq.gz",
            "--threads",
            expected=None,
        )
        self.assertNotEqual(missing.returncode, 0)
        self.assertIn("requires a value", missing.stderr)

        missing_cutoff = self.run_stage(
            STAGE2,
            "--input-fastq",
            "input.fastq.gz",
            "--enable-chopper-trim",
            "--chopper-trim-approach",
            "trim-by-quality",
            expected=None,
        )
        self.assertNotEqual(missing_cutoff.returncode, 0)
        self.assertIn("--chopper-trim-cutoff", missing_cutoff.stderr)

        invalid_config = self.work / "invalid.conf"
        invalid_config.write_text("enable_fastcat_lint=maybe\n", encoding="utf-8")
        invalid_bool = self.run_stage(
            STAGE2,
            "--input-fastq",
            "input.fastq.gz",
            "--config",
            "invalid.conf",
            expected=None,
        )
        self.assertNotEqual(invalid_bool.returncode, 0)
        self.assertIn("must be true or false", invalid_bool.stderr)

        valid_fraction = self.run_stage(
            STAGE2,
            "--input-fastq",
            "input.fastq.gz",
            "--keep-percent",
            "99.9",
            "--dry-run",
        )
        self.assertIn("--keep_percent 99.9", valid_fraction.stdout)

        alias = self.run_stage(
            STAGE2,
            "--input-fastq",
            "input.fastq.gz",
            "--output-dir",
            ".",
            "--output-fastq",
            "input.fastq.gz",
            expected=None,
        )
        self.assertNotEqual(alias.returncode, 0)
        self.assertIn("same file", alias.stderr)

    def test_failed_filter_never_publishes_final_output_or_manifest(self) -> None:
        filtlong = self.work / "mock_bin" / "filtlong"
        filtlong.write_text("#!/bin/sh\nexit 9\n", encoding="ascii", newline="\n")
        result = self.run_stage(
            STAGE2,
            "--input-fastq",
            "input.fastq.gz",
            "--output-dir",
            "failed-output",
            "--keep-percent",
            "90",
            expected=None,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(
            (self.work / "failed-output" / "filtered_input.fastq.gz").exists()
        )
        self.assertFalse(
            (self.work / "failed-output" / ".prokaryont_stage2.complete.tsv").exists()
        )

    def test_failed_post_filter_diagnostic_preserves_previous_final(self) -> None:
        output = self.work / "diagnostic-failure"
        output.mkdir()
        previous = output / "filtered_input.fastq.gz"
        previous_bytes = (self.work / "input.fastq.gz").read_bytes()
        previous.write_bytes(previous_bytes)
        nanoplot = self.work / "mock_bin" / "NanoPlot"
        nanoplot.write_text(
            "#!/bin/sh\n"
            'if [ "${1:-}" = "--version" ]; then echo "NanoPlot mock"; exit 0; fi\n'
            "exit 9\n",
            encoding="ascii",
            newline="\n",
        )

        result = self.run_stage(
            STAGE2,
            "--input-fastq",
            "input.fastq.gz",
            "--output-dir",
            "diagnostic-failure",
            "--force",
            expected=None,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(previous.read_bytes(), previous_bytes)
        self.assertFalse((output / ".prokaryont_stage2.complete.tsv").exists())

    def test_plain_fastq_and_spaces_in_output_path(self) -> None:
        with gzip.open(self.work / "input.fastq.gz", "rb") as source:
            (self.work / "input.fastq").write_bytes(source.read())
        self.run_stage(
            STAGE2,
            "--input-fastq",
            "input.fastq",
            "--output-dir",
            "output with spaces",
            "--skip-nanoplot",
            "--skip-fastplong",
        )
        filtered = self.work / "output with spaces" / "filtered_input.fastq.gz"
        self.assertTrue(filtered.is_file())
        subprocess.run(
            [sys.executable, str(CONTRACT), "validate-fastq", str(filtered)],
            check=True,
        )

    def test_mock_fastcat_rejects_the_obsolete_interface(self) -> None:
        result = subprocess.run(
            [
                bash_executable(),
                str(self.work / "mock_bin" / "fastcat"),
                "--histograms",
                "legacy",
            ],
            cwd=self.work,
            env=self.environment,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)


class Stage3GenomeSizeContractTests(unittest.TestCase):
    def test_stage3_does_not_consume_weighted_stage1_estimate(self) -> None:
        source = (ROOT / "03_autocycler_assemble.sh").read_text(encoding="utf-8")
        self.assertNotIn("mean_genome_size.txt", source)
        self.assertIn("Estimating genome size from the selected filtered reads", source)


if __name__ == "__main__":
    unittest.main()
