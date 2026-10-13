import os
import tempfile
import unittest

from funannotate import library as lib

DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "deeptmhmm")
TMRS_GFF3 = os.path.join(DATA_DIR, "DeepTMHMM.TMRs.gff3")
TOPOLOGIES_3LINE = os.path.join(DATA_DIR, "DeepTMHMM.predicted_topologies.3line")


class ParseTMHMMTests(unittest.TestCase):
    def test_real_deeptmhmm_gff3_is_parsed_without_error(self):
        # DeepTMHMM pads each data line with trailing empty tab-separated
        # fields (id, type, start, end, "", "", "", ""); parseTMHMM must not
        # mistake that for 9-column GFF3 (id, source, type, start, end).
        with tempfile.TemporaryDirectory() as tmpdir:
            membrane_annot = os.path.join(tmpdir, "annotations.transmembrane.txt")
            lib.parseTMHMM(TMRS_GFF3, membrane_annot)

            with open(membrane_annot) as handle:
                lines = [line.rstrip("\n") for line in handle if line.strip()]

        self.assertEqual(len(lines), 1)
        pid, note_tag, note = lines[0].split("\t")
        self.assertEqual(pid, "5H2A_CRIGR")
        self.assertEqual(note_tag, "note")
        self.assertTrue(note.startswith("TransMembrane:7 ("))
        # topology must alternate orientation/segment and use 1-based
        # residue coordinates copied straight from the gff3 TMhelix rows
        self.assertEqual(
            note,
            "TransMembrane:7 (o77-101i108-133o150-171i192-212o234-256i324-348o360-382i)",
        )

    def test_protein_with_no_tmhelix_segments_is_skipped(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gff3 = os.path.join(tmpdir, "TMRs.gff3")
            with open(gff3, "w") as handle:
                handle.write("##gff-version 3\n")
                handle.write("# SOLUBLE_PROT Length: 100\n")
                handle.write("# SOLUBLE_PROT Number of predicted TMRs: 0\n")
                handle.write("SOLUBLE_PROT\toutside\t1\t100\t\t\t\t\n")

            membrane_annot = os.path.join(tmpdir, "annotations.transmembrane.txt")
            lib.parseTMHMM(gff3, membrane_annot)

            with open(membrane_annot) as handle:
                content = handle.read()

        self.assertEqual(content, "")

    def test_predicted_topologies_3line_reference_matches_gff3_call_count(self):
        # Cross-check: the human-readable topology string ("3line" format)
        # for the same protein reports the same number of TM segments (M
        # runs) that the gff3-derived note reports.
        with open(TOPOLOGIES_3LINE) as handle:
            record = [next(handle).rstrip("\n") for _ in range(3)]

        header, _sequence, topology_string = record
        self.assertTrue(header.startswith(">5H2A_CRIGR"))

        # count contiguous runs of "M" (membrane) characters
        tm_runs = 0
        in_run = False
        for char in topology_string:
            if char == "M":
                if not in_run:
                    tm_runs += 1
                    in_run = True
            else:
                in_run = False

        with tempfile.TemporaryDirectory() as tmpdir:
            membrane_annot = os.path.join(tmpdir, "annotations.transmembrane.txt")
            lib.parseTMHMM(TMRS_GFF3, membrane_annot)
            with open(membrane_annot) as handle:
                note = handle.read().strip()

        self.assertEqual(tm_runs, 7)
        self.assertIn("TransMembrane:7 ", note)


if __name__ == "__main__":
    unittest.main()
