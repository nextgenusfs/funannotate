"""#1240: minimap2 refuses a genome larger than one index part unless -I covers it.

--split-prefix is not an option: minimap2 rejects it together with --cs, and the BAM
parsers need the cs tag.
"""
import os
import tempfile
import unittest
from unittest import mock

from funannotate import library as lib


class Minimap2LargeGenomeTests(unittest.TestCase):
    def test_small_genome_unchanged(self):
        with tempfile.NamedTemporaryFile(suffix=".fa") as f:
            f.write(b">c\nACGT\n")
            f.flush()
            mm2, st = lib.minimap2_align_cmds(f.name, "t.fa", 8, 3000, "x.bam")
        self.assertNotIn("-I", mm2)
        self.assertNotIn("--split-prefix", mm2)
        self.assertIn("--cs", mm2)
        self.assertEqual(st[st.index("-@") + 1], "4")

    def test_large_genome_gets_index_part_covering_it(self):
        with mock.patch.object(os.path, "getsize", return_value=30_400_000_000):
            self.assertEqual(lib.minimap2_index_part_size("g.fa"), "32G")
            mm2, _ = lib.minimap2_align_cmds("g.fa", "t.fa", 48, 20000, "x.bam")
        self.assertEqual(mm2[mm2.index("-I") + 1], "32G")
        self.assertIn("--cs", mm2)
        self.assertEqual(mm2[-2:], ["g.fa", "t.fa"])

    def test_missing_genome_is_not_an_error(self):
        self.assertIsNone(lib.minimap2_index_part_size("/no/such/genome.fa"))


if __name__ == "__main__":
    unittest.main()
