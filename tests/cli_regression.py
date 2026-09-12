"""End-to-end regression checks using temporary inputs and outputs only."""

import gzip
import itertools
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


EXECUTABLE = Path(sys.argv.pop(1) if len(sys.argv) > 1 else "./LCA").resolve()
LINEAGE = "k__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__Genus {}"
HEADER = "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tqlen"


def hit(query="q", subject="a", identity=99, length=100, qlen=100):
    return "\t".join(map(str, [query, subject, identity, length, 0, 0, 1, length, 1, length, qlen])) + "\n"


class CLIRegression(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="lca-regression-")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.mapping = self.directory / "input.m8"
        self.taxonomy = self.directory / "reference.tax"
        self.output = self.directory / "output.tsv"
        self.mapping.write_text(hit())
        self.taxonomy.write_text("a\t" + LINEAGE.format("alpha") + "\nb\t" + LINEAGE.format("beta") + "\n")

    def run_lca(self, *arguments, expected=0, mapping=None, taxonomy=None, output=None):
        result = subprocess.run(
            [str(EXECUTABLE), "-i", str(mapping or self.mapping), "-r", str(taxonomy or self.taxonomy),
             "-o", str(output or self.output), *map(str, arguments)],
            capture_output=True, text=True, timeout=20,
        )
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        return result

    def rows(self):
        return [line.split("\t") for line in self.output.read_text().splitlines()]

    def test_tab_header_numeric_padding_and_crlf(self):
        fields = hit().rstrip("\n").split("\t")
        fields[2:] = [" " + field + " " for field in fields[2:]]
        self.mapping.write_bytes(("  # comment\r\n\t \r\n" + HEADER + "\r\n" +
                                  "\t".join(fields) + "\r\n").encode())
        self.run_lca()
        self.assertEqual(self.rows()[1], ["q", "Bacteria", "P", "C", "O", "F", "Genus", "Genus alpha"])

    def test_legacy_coverage_is_skipped(self):
        header = HEADER.replace("qlen", "evalue\tbitscore")
        self.mapping.write_text(header + "\nq\ta\t99\t100\t0\t0\t101\t200\t1\t100\t1e-35\t152\n")
        self.run_lca("-cover", "1")
        self.assertEqual(self.rows()[1][-1], "Genus alpha")

    def test_malformed_and_mixed_layouts_report_line(self):
        for invalid in [hit().replace("99", "nan", 1), "q a 99\n",
                        "q\ta\t99\t100\t0\t0\t1\t100\t1\t100\t1e-20\t80\n",
                        "q\ta\t99\t100\t0\t0\t0\t100\t1\t100\t100\n"]:
            with self.subTest(row=invalid):
                self.mapping.write_text(hit("first") + invalid)
                result = self.run_lca(expected=25)
                self.assertIn("line 2", result.stderr)

    def test_reported_illumina_header_preserves_query_identity(self):
        spot = "LH00409:413:22WTJGLT4:6:2481:9713:19847"
        annotated = spot + " 1:N:0:CGTATCTC+CTCGAACA"
        other_annotation = spot + " 2:N:0:CGTATCTC+CTCGAACA"
        subject = "OBEP010038528"
        self.taxonomy.write_text(subject + "\t" + LINEAGE.format("alpha") + "\n")
        reported = annotated + "\tOBEP010038528\t 95.3 \t150\t7\t0\t1\t150\t795\t944\t150\n"
        self.mapping.write_text(hit(spot, subject, 95.3, 150, 150) + reported +
                                hit(other_annotation, subject, 95.3, 150, 150))
        self.run_lca("-showHitRead", "-reportID")
        rows = self.rows()[1:]
        self.assertEqual([row[0] for row in rows], [spot, annotated, other_annotation])
        self.assertTrue(all(row[-2:] == [subject, "95.300000"] for row in rows))

    def test_subject_descriptions_are_preserved(self):
        first, second = "reference description A", "reference description B"
        self.taxonomy.write_text(first + "\t" + LINEAGE.format("alpha") + "\n" +
                                 second + "\t" + LINEAGE.format("beta") + "\n")
        self.mapping.write_text(hit(subject=first) + hit(subject=second))
        self.run_lca("-showHitRead")
        self.assertEqual(self.rows()[1][-2:], ["?", "?"])
        self.run_lca("-reportBestHit", "-showHitRead")
        self.assertEqual(self.rows()[1][-2:], ["Genus alpha", first])

    def test_legacy_identifiers_with_spaces(self):
        query, subject = "query with barcode", "reference with description"
        self.taxonomy.write_text(subject + "\t" + LINEAGE.format("alpha") + "\n")
        self.mapping.write_text(HEADER.replace("qlen", "evalue\tbitscore") + "\n" +
                                query + "\t" + subject + "\t99\t100\t0\t0\t1\t100\t1\t100\t 1e-35 \t 152 \n")
        self.run_lca("-showHitRead")
        self.assertEqual(self.rows()[1][0], query)
        self.assertEqual(self.rows()[1][-2:], ["Genus alpha", subject])

    def test_empty_tab_fields_fail(self):
        for column in range(11):
            fields = hit().rstrip("\n").split("\t")
            fields[column] = ""
            self.mapping.write_text("\t".join(fields) + "\n")
            result = self.run_lca(expected=25)
            self.assertIn("line 1", result.stderr)

    def test_noncontiguous_queries_fail(self):
        self.mapping.write_text(hit("q1") + hit("q2") + hit("q1"))
        self.run_lca(expected=27)

    def test_best_hit_taxonomy_matches_reported_subject(self):
        for order, filtering in itertools.product(itertools.permutations(["a", "b"]), [[], ["-no_bl_filter"]]):
            with self.subTest(order=order, filtering=filtering):
                self.mapping.write_text("".join(hit(subject=subject) for subject in order))
                self.run_lca("-reportBestHit", "-showHitRead", *filtering)
                self.assertEqual(self.rows()[1][-2:], ["Genus alpha", "a"])

    def test_ineligible_duplicates_do_not_hide_valid_hits(self):
        for ineligible in [hit(identity=100, length=50), hit(identity=100, qlen=300)]:
            for order in itertools.permutations([ineligible, hit()]):
                with self.subTest(order=order):
                    self.mapping.write_text("".join(order))
                    self.run_lca("-reportID")
                    self.assertEqual(self.rows()[1][-2:], ["Genus alpha", "99.000000"])

    def test_relative_length_filter_precedes_deduplication(self):
        alignments = [hit(identity=97.5, length=100), hit(subject="b", identity=97.5, length=80),
                      hit(subject="b", identity=97, length=100)]
        for order in itertools.permutations(alignments):
            self.mapping.write_text("".join(order))
            self.run_lca("-reportID")
            self.assertEqual(self.rows()[1][-2:], ["?", "97.250000"])

    def test_duplicate_subjects_do_not_gain_extra_votes(self):
        self.mapping.write_text(hit() * 20 + hit(subject="b"))
        self.run_lca("-showHitRead")
        self.assertEqual(self.rows()[1][-2:], ["?", "?"])

    def test_exact_length_threshold_is_inclusive(self):
        self.mapping.write_text(hit(identity=97) + hit(subject="b", identity=97, length=85))
        self.run_lca()
        self.assertEqual(self.rows()[1][-1], "?")

    def test_multi_database_selection_uses_eligible_ranks(self):
        second_mapping = self.directory / "second.m8"
        second_taxonomy = self.directory / "second.tax"
        self.mapping.write_text(hit(identity=94))
        second_mapping.write_text(hit(subject="c", identity=96))
        second_taxonomy.write_text("c\tk__Bacteria;p__P;c__C;o__O;f__F;g__BetterGenus\n")
        for mappings, taxonomies in [([self.mapping, second_mapping], [self.taxonomy, second_taxonomy]),
                                    ([second_mapping, self.mapping], [second_taxonomy, self.taxonomy])]:
            self.run_lca("-reportID", mapping=",".join(map(str, mappings)), taxonomy=",".join(map(str, taxonomies)))
            self.assertEqual(self.rows()[1][-3:], ["BetterGenus", "?", "96.000000"])

    def test_identity_mask_matches_matrix_and_hit_pattern(self):
        pattern = self.directory / "pattern.tsv"
        self.mapping.write_text(hit(identity=97) + hit(subject="b", identity=96))
        self.run_lca("-matHigh", "-reportHitPattern", pattern)
        self.assertEqual(self.rows()[1][-1], "?")
        self.assertEqual(Path(str(self.output) + "_Species").read_text(), "Bacteria;P;C;O;F;Genus;?\t1\n")
        self.assertEqual(pattern.read_text().splitlines()[1].split("\t")[-1], "6")

    def test_decimal_identity_threshold_is_inclusive(self):
        self.taxonomy.write_text("".join(f"s{i}\t" + LINEAGE.format("alpha") + "\n" for i in range(30)))
        for count in [1, 9, 30]:
            self.mapping.write_text("".join(hit(subject=f"s{i}", identity=97.2) for i in range(count)))
            self.run_lca("-id", "97.2,95,93,91,88,78,0", "-reportID")
            self.assertEqual(self.rows()[1][-2:], ["Genus alpha", "97.200000"])

    def test_missing_species_is_not_evidence_of_a_certain_species(self):
        self.taxonomy.write_text("a\t" + LINEAGE.format("sp.") + ";t__strainA\n" +
                                 "b\tk__Bacteria;p__P;c__C;o__O;f__F;g__Genus;s__;t__strainB\n")
        self.mapping.write_text(hit() + hit(subject="b"))
        self.run_lca("-tdep", "8")
        self.assertEqual(self.rows()[1][-2], "Genus sp.")

    def test_empty_and_single_final_line_inputs(self):
        self.mapping.write_text("")
        self.run_lca()
        self.assertEqual(len(self.rows()), 1)
        self.mapping.write_text(hit().rstrip("\n"))
        self.run_lca()
        self.assertEqual(len(self.rows()), 2)

    @unittest.skipIf(os.environ.get("LCA_TEST_NO_GZIP"), "gzip disabled in this build")
    def test_gzip_success_and_corruption(self):
        compressed = self.directory / "input.m8.gz"
        good = gzip.compress(hit().encode())
        compressed.write_bytes(good)
        self.run_lca(mapping=compressed)
        self.assertEqual(len(self.rows()), 2)
        for broken in [good[:-8], good[:-8] + b"\0" * 8]:
            with self.subTest(compressed=broken):
                compressed.write_bytes(broken)
                self.run_lca(mapping=compressed, expected=28)

    @unittest.skipUnless(os.name == "posix", "directory read errors are platform dependent")
    def test_read_errors_fail(self):
        self.run_lca(mapping=self.directory, expected=28)
        self.run_lca(taxonomy=self.directory, expected=13)

    def test_missing_subject_fails(self):
        self.mapping.write_text(hit(subject="absent"))
        self.run_lca(expected=74)

    def test_overlapping_outputs_preserve_files(self):
        matrix = Path(str(self.output) + "_Domain")
        matrix.write_text(hit())
        self.output.write_text("existing output\n")
        original = {path: path.read_bytes() for path in [self.mapping, self.taxonomy, self.output, matrix]}
        cases = [([], {"output": self.mapping}), ([], {"output": self.taxonomy}),
                 (["-reportHitPattern", self.mapping], {}), (["-reportHitPattern", self.output], {}),
                 (["-matHigh"], {"mapping": matrix}), (["-matHigh", "-reportHitPattern", matrix], {})]
        for args, kwargs in cases:
            with self.subTest(args=args, kwargs=kwargs):
                self.run_lca(*args, **kwargs, expected=29)
                for path, content in original.items():
                    self.assertEqual(path.read_bytes(), content)

    @unittest.skipUnless(os.name == "posix", "link creation needs platform support")
    def test_input_aliases_are_protected(self):
        original = self.mapping.read_bytes()
        symlink = self.directory / "alias.m8"
        symlink.symlink_to(self.mapping)
        hardlink = self.directory / "hardlink.m8"
        os.link(self.mapping, hardlink)
        for alias in [symlink, hardlink, str(self.directory) + "/./input.m8"]:
            self.run_lca(output=alias, expected=29)
            self.assertEqual(self.mapping.read_bytes(), original)

    @unittest.skipUnless(Path("/dev/full").exists(), "requires /dev/full")
    def test_write_errors_fail(self):
        self.run_lca(output="/dev/full", expected=32)
        self.run_lca("-reportHitPattern", "/dev/full", expected=33)
        Path(str(self.output) + "_Domain").symlink_to("/dev/full")
        self.run_lca("-matHigh", expected=34)


if __name__ == "__main__":
    unittest.main()
