import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class DocumentationCliContractTests(unittest.TestCase):
    def test_stage1_stage2_layout_and_identity_flags_are_documented(self) -> None:
        for filename in ("01_qc_estimate.sh", "02_preprocess_filter.sh"):
            with self.subTest(filename=filename):
                source = (ROOT / filename).read_text(encoding="utf-8")
                for option in (
                    "--output-dir",
                    "--sample-name",
                    "--tmp-dir",
                    "--input-sha256",
                    "--force",
                    "--dry-run",
                ):
                    self.assertIn(option, source)
        self.assertIn(
            "--output-fastq",
            (ROOT / "02_preprocess_filter.sh").read_text(encoding="utf-8"),
        )

    def test_stage1_stage2_tool_interfaces_are_pinned(self) -> None:
        stage1 = (ROOT / "envs" / "env_stage1_qc.yml").read_text()
        stage2 = (ROOT / "envs" / "env_stage2_preproc.yml").read_text()
        self.assertIn("fastcat>=1.0.1,<2", stage1)
        self.assertIn("fastcat>=1.0.1,<2", stage2)
        self.assertIn("chopper>=0.13.0,<0.14", stage2)

    def test_stage3_default_and_optional_synopsis(self) -> None:
        source = (ROOT / "03_autocycler_assemble.sh").read_text(encoding="utf-8")
        self.assertIn('assembler_list="flye,canu,hifiasm,miniasm,myloasm"', source)
        self.assertIn('echo "Usage: $(basename "$0") [OPTIONS]"', source)

    def test_stage4_read_group_is_validated_and_forwarded_twice(self) -> None:
        source = (ROOT / "04_polish_orient.sh").read_text(encoding="utf-8")
        self.assertIn("--read-group)", source)
        self.assertIn("validate_read_group_in_bam aligned.sorted.bam", source)
        self.assertIn(
            "validate_read_group_in_bam aligned_reoriented.sorted.bam", source
        )
        self.assertIn('dorado_rg_flag=(--RG "${dorado_read_group}")', source)
        self.assertIn('"${field#ID:}" == "${requested_id}"', source)
        self.assertEqual(source.count('"${dorado_rg_flag[@]}"'), 2)

    def test_stage_environment_files_are_committed_artifacts(self) -> None:
        expected = {
            "env_stage1_qc.yml": "prokaryont-stage1-qc",
            "env_stage2_preproc.yml": "prokaryont-stage2-preproc",
            "env_stage3_assemble.yml": "prokaryont-stage3-assemble",
            "env_stage4_polish.yml": "prokaryont-stage4-polish",
            "env_stage5_taxonomy.yml": "prokaryont-stage5-taxonomy",
        }
        for filename, environment_name in expected.items():
            with self.subTest(filename=filename):
                path = ROOT / "envs" / filename
                self.assertTrue(path.is_file())
                self.assertIn(f"name: {environment_name}", path.read_text())

    def test_shipped_docs_do_not_contain_machine_local_links(self) -> None:
        for filename in (
            "README.md",
            "conda-environment-setup-for-each-prokaryont-script.md",
            "testing_guidelines.md",
        ):
            with self.subTest(filename=filename):
                source = (ROOT / filename).read_text(encoding="utf-8")
                self.assertNotIn("file:///", source)


if __name__ == "__main__":
    unittest.main()
