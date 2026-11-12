#!/usr/bin/env python3
"""Real Nextflow command/cache tests; invented empty inputs, no Docker or caller.

Runs the production UPSTREAM workflow and DEEPVARIANT script, replacing only
the external goleft/DeepVariant executables with local argument recorders.
This deliberately does not use -stub: the real calling command must be tested.
"""
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class SmallModelTest(unittest.TestCase):
    def test_commands_and_resume(self):
        with tempfile.TemporaryDirectory(prefix="og-dv-test-") as tmp:
            work = Path(tmp)
            ref = work / "reference.fasta"
            ref.write_text(">chr20\nA\n")
            Path(f"{ref}.fai").write_text("chr20\t1\t7\t1\t2\n")
            bam = work / "example.bam"
            bai = work / "example.bam.bai"
            bam.write_text("invented placeholder, not a BAM\n")
            bai.write_text("invented placeholder, not a BAI\n")
            sheet = work / "bam-only.csv"
            sheet.write_text(f"sample,bam,bam_index,sex\nexample,{bam},{bai},male\n")
            fake = work / "commands.py"
            fake.write_text('''import gzip, json, pathlib, sys
mode, *args = sys.argv[1:]
if mode == 'goleft':
    assert args[:3] == ['indexcov', '--directory', 'indexcov'], args
    pathlib.Path('indexcov').mkdir()
    pathlib.Path('indexcov/indexcov-indexcov.ped').write_text(
        '#family_id sample_id sex CNchrX CNchrY\\nexample example 1 1 1\\n')
else:
    pathlib.Path('argv.json').write_text(json.dumps(args))
    for key in ('--output_vcf=', '--output_gvcf='):
        output, = [a[len(key):] for a in args if a.startswith(key)]
        with gzip.open(output, 'wt') as f:
            f.write('placeholder output, no variant calls\\n')
        pathlib.Path(output + '.tbi').touch()
''')
            harness = work / "main.nf"
            harness.write_text(f'''include {{ UPSTREAM }} from '{ROOT}/workflows/upstream'
workflow {{
    rows = file(params.input).splitCsv(header: true)
    bams = Channel.fromList(rows).map {{ row ->
        [[id: row.sample, sex: row.sex], file(row.bam), file(row.bam_index)]
    }}
    UPSTREAM(Channel.empty(), bams, Channel.empty(),
        Channel.value(file(params.reference)),
        Channel.value(file(params.reference + '.fai')),
        Channel.value(file('{ROOT}/assets/par_grch38.bed')))
}}
''')
            shell_env = work / "command-env.sh"
            shell_env.write_text(f'''
function /opt/deepvariant/bin/run_deepvariant() {{ python3 '{fake}' dv "$@"; }}
function goleft() {{ python3 '{fake}' goleft "$@"; }}
''')
            override = work / "local.config"
            override.write_text(f'''docker.enabled = false
singularity.enabled = false
timeline.enabled = false
report.enabled = false
dag.enabled = false
trace.enabled = false
process {{
    executor = 'local'
    withName: '.*' {{
        cpus = 2
        memory = '256.MB'
        beforeScript = "export BASH_ENV='{shell_env}'"
    }}
}}
''')
            env = os.environ.copy()
            env["NXF_ANSI_LOG"] = "false"

            def run(name, flags=(), resume=False):
                trace = work / f"{name}.tsv"
                cmd = ["nextflow", "-C", f"{ROOT}/nextflow.config,{override}",
                       "run", str(harness), "-ansi-log", "false", "-with-trace", str(trace),
                       "--input", str(sheet), "--reference", str(ref),
                       "--outdir", str(work / "results"),
                       "--intervals", "chr20:1-10 chr22:1-10", *flags]
                if resume:
                    cmd.append("-resume")
                result = subprocess.run(cmd, cwd=work, env=env, text=True,
                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                        timeout=180)
                self.assertEqual(result.returncode, 0, result.stdout)
                with trace.open() as f:
                    rows = {r["name"].split()[0].split(":")[-1]: r
                            for r in csv.DictReader(f, delimiter="\t")}
                self.assertEqual(set(rows), {"INDEXCOV", "DEEPVARIANT"}, result.stdout)
                return rows

            def argv_for(row):
                prefix, suffix = row["hash"].split("/")
                task, = (work / "work" / prefix).glob(suffix + "*")
                return json.loads((task / "argv.json").read_text()), task

            default = run("default")
            args, task = argv_for(default["DEEPVARIANT"])
            expected = ["--model_type=WGS", "--ref=reference.fasta", "--reads=example.bam",
                        "--output_vcf=example.vcf.gz", "--output_gvcf=example.g.vcf.gz",
                        "--intermediate_results_dir=deepvariant_tmp", "--sample_name=example",
                        "--num_shards=2", "--haploid_contigs=chrX,chrY",
                        "--par_regions_bed=par_grch38.bed", "--regions", "chr20:1-10 chr22:1-10"]
            self.assertEqual(args, expected)
            self.assertIn("disable_small_model=false", (task / ".command.out").read_text())

            resumed = run("default-again", resume=True)
            for process in default:
                self.assertEqual(resumed[process]["hash"], default[process]["hash"], process)
                self.assertEqual(resumed[process]["status"], "CACHED")

            # Nextflow can hash a config Boolean and a CLI string differently,
            # even though their effective command is identical.
            explicit_false = run("false", ["--deepvariant_disable_small_model", "false"], True)
            self.assertEqual(argv_for(explicit_false["DEEPVARIANT"])[0], expected)
            self.assertEqual(explicit_false["INDEXCOV"]["status"], "CACHED")

            disabled = run("disabled", ["--deepvariant_disable_small_model", "true"], True)
            args, task = argv_for(disabled["DEEPVARIANT"])
            self.assertEqual(args.count("--disable_small_model=true"), 1)
            self.assertEqual([a for a in args if a != "--disable_small_model=true"], expected)
            self.assertIn("disable_small_model=true", (task / ".command.out").read_text())
            self.assertNotEqual(disabled["DEEPVARIANT"]["hash"], default["DEEPVARIANT"]["hash"])
            self.assertEqual(disabled["DEEPVARIANT"]["status"], "COMPLETED")
            self.assertEqual(disabled["INDEXCOV"]["hash"], default["INDEXCOV"]["hash"])
            self.assertEqual(disabled["INDEXCOV"]["status"], "CACHED")
            repeated = run("disabled-again", ["--deepvariant_disable_small_model", "true"], True)
            for process in repeated:
                self.assertEqual(repeated[process]["hash"], disabled[process]["hash"])
                self.assertEqual(repeated[process]["status"], "CACHED")

            bare = run("bare", ["--deepvariant_disable_small_model"], True)
            self.assertEqual(argv_for(bare["DEEPVARIANT"])[0], args)

            # Also exercise the public entry point and its samplesheet routing.
            trace = work / "pipeline.tsv"
            result = subprocess.run(
                ["nextflow", "-C", f"{ROOT}/nextflow.config,{override}",
                 "run", str(ROOT / "main.nf"), "-stub", "-ansi-log", "false",
                 "-with-trace", str(trace), "--input", str(sheet),
                 "--reference", str(ref), "--outdir", str(work / "pipeline-results"),
                 "--tools", " ", "--deepvariant_disable_small_model", "true"],
                cwd=work, env=env, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, timeout=180)
            self.assertEqual(result.returncode, 0, result.stdout)
            with trace.open() as f:
                processes = {r["name"].split()[0].split(":")[-1]
                             for r in csv.DictReader(f, delimiter="\t")}
            self.assertTrue({"INDEXCOV", "DEEPVARIANT"} <= processes, processes)
            self.assertFalse({"FASTP", "MINIMAP2_INDEX", "ALIGN_MINIMAP2"} & processes,
                             processes)

    def test_invalid_value_rejected_at_pipeline_launch(self):
        with tempfile.TemporaryDirectory(prefix="og-dv-invalid-") as tmp:
            result = subprocess.run(
                ["nextflow", "run", str(ROOT / "main.nf"), "-ansi-log", "false",
                 "--input", "unused.csv", "--reference", "unused.fasta",
                 "--deepvariant_disable_small_model", "typo"], cwd=tmp, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--deepvariant_disable_small_model takes no value (or true / false)",
                          result.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
