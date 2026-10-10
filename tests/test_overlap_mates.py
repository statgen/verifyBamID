#!/usr/bin/env python3
"""Regression for overlap filtering, using the actual binary and indexed BAMs.

Requires Python 3 and samtools. Run after building verifyBamID:
  python3 tests/test_overlap_mates.py [--verifybamid bin/verifyBamID]
"""
import argparse
import csv
from pathlib import Path
import subprocess
import tempfile
import unittest


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--verifybamid", type=Path,
                    default=Path(__file__).resolve().parents[1] / "bin/verifyBamID")
parser.add_argument("--samtools", default="samtools")
args = parser.parse_args()


class OverlappingMatesTest(unittest.TestCase):
    def count_reads(self, observations, ignore_overlap=False):
        with tempfile.TemporaryDirectory(prefix="verifybamid-overlap-") as directory:
            directory = Path(directory)
            sam, bam = directory / "reads.sam", directory / "reads.bam"
            vcf, output = directory / "sites.vcf", directory / "result"
            rows = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:1\tLN:1000",
                    "@RG\tID:rg1\tSM:sample"]
            for i, (kind, name) in enumerate(observations):
                # All reads start at 90, so query offset 10 is marker 100.
                sequence, quality, cigar = "A" * 21, "I" * 21, "21M"
                if kind == "low_quality":
                    quality = "I" * 10 + "!" + "I" * 10
                elif kind == "unknown_base":
                    sequence = "A" * 10 + "N" + "A" * 10
                elif kind == "deletion":
                    sequence, quality, cigar = "A" * 20, "I" * 20, "10M1D10M"
                flag = 99 if i == 0 else 147
                rows.append("\t".join(map(str, [name, flag, "1", 90, 60, cigar,
                                                "=", 90, 21 if i == 0 else -21,
                                                sequence, quality, "RG:Z:rg1"])))
            sam.write_text("\n".join(rows) + "\n")
            subprocess.run([args.samtools, "view", "-b", "-o", str(bam), str(sam)],
                           check=True, capture_output=True)
            subprocess.run([args.samtools, "index", str(bam)], check=True,
                           capture_output=True)
            vcf.write_text(
                '##fileformat=VCFv4.2\n'
                '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n'
                '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                '1\t100\t.\tA\tC\t.\tPASS\tAF=0.5\n')
            command = [str(args.verifybamid.resolve()), "--bam", str(bam), "--vcf", str(vcf),
                       "--out", str(output), "--site", "--free-none", "--chip-none",
                       "--noPhoneHome"]
            if ignore_overlap:
                command.append("--ignoreOverlapPair")
            process = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            with output.with_suffix(".selfSM").open() as result:
                return int(next(csv.DictReader(result, delimiter="\t"))["#READS"])

    def test_rejected_mate_does_not_consume_name(self):
        for rejected in ("low_quality", "unknown_base", "deletion"):
            for observations in ([(rejected, "pair"), ("good", "pair")],
                                 [("good", "pair"), (rejected, "pair")]):
                with self.subTest(observations=observations):
                    self.assertEqual(self.count_reads(observations), 1)

    def test_two_valid_mates_still_count_once(self):
        self.assertEqual(self.count_reads([("good", "pair"), ("good", "pair")]), 1)

    def test_overlap_opt_out_counts_both_mates(self):
        self.assertEqual(self.count_reads([("good", "pair"), ("good", "pair")], True), 2)

    def test_independent_reads_are_retained(self):
        self.assertEqual(self.count_reads([("good", "read1"), ("good", "read2")]), 2)


if __name__ == "__main__":
    unittest.main(argv=[__file__])
