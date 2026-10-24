"""NCBI translation tables: the code tables, translate(), and gene models.

Covers what `--table` changes without running any external tool:

* funannotate.genetic_codes: which table ids are valid, stop and start codons;
* library.translate(): codon -> amino acid for alternative tables, both strands, phases;
* library.convertgff2tbl(): a GFF3 gene model becomes a protein and a .tbl record whose
  `transl_table` qualifier is written only for non-standard tables, and an in-frame TGA is
  read as Trp (table 4) or kept as a stop (table 1).

Tables used and why: 12 (alternative yeast nuclear, CUG = Ser) for CUG-clade yeasts;
4 (TGA = Trp) as the stop-reassigning table (mold mitochondrial / Mycoplasma); 6 (TAA, TAG = Gln,
ciliate nuclear), the other stop-reassigning table GeneMark accepts.
"""
import logging
import os
import random
import shutil
import tempfile
import unittest

from funannotate import library as lib
from funannotate.genetic_codes import (
    get_codon_table,
    get_start_codons,
    get_stop_codons,
    is_valid_table,
    valid_tables,
)


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class IsValidTableTests(unittest.TestCase):
    def test_common_tables_are_valid(self):
        for table in (1, 2, 3, 4, 5, 6, 11, 12, 26):
            self.assertTrue(is_valid_table(table), table)

    def test_numeric_strings_are_accepted(self):
        self.assertTrue(is_valid_table("12"))

    def test_unassigned_and_malformed_ids_are_rejected(self):
        # NCBI has no tables 7, 8 and 17-20; 0 and negatives are not ids
        for table in (0, 7, 8, 17, 18, 19, 20, 99, -1, "abc", "", None, 1.5):
            self.assertFalse(is_valid_table(table), repr(table))

    def test_valid_tables_are_sorted_and_include_standard(self):
        tables = valid_tables()
        self.assertEqual(tables, sorted(tables))
        self.assertIn(1, tables)
        self.assertNotIn(7, tables)


class StopAndStartCodonTests(unittest.TestCase):
    def test_standard_stops(self):
        self.assertEqual(get_stop_codons(1), {"TAA", "TAG", "TGA"})

    def test_table_4_reads_tga_as_tryptophan(self):
        self.assertEqual(get_stop_codons(4), {"TAA", "TAG"})

    def test_table_6_reads_taa_and_tag_as_glutamine(self):
        self.assertEqual(get_stop_codons(6), {"TGA"})

    def test_table_12_keeps_standard_stops(self):
        # CUG-Ser changes a sense codon only, so no stop-codon warning applies
        self.assertEqual(get_stop_codons(12), get_stop_codons(1))

    def test_atg_starts_every_common_table(self):
        for table in (1, 4, 6, 11, 12):
            self.assertIn("ATG", get_start_codons(table), table)

    def test_table_11_allows_gtg_start(self):
        self.assertIn("GTG", get_start_codons(11))

    def test_string_ids_work(self):
        self.assertEqual(get_stop_codons("4"), get_stop_codons(4))


class CodonTableTests(unittest.TestCase):
    def test_every_table_covers_all_64_codons(self):
        for table in valid_tables():
            self.assertEqual(len(get_codon_table(table)), 64, table)

    def test_stops_map_to_asterisk(self):
        codons = get_codon_table(1)
        for stop in ("TAA", "TAG", "TGA"):
            self.assertEqual(codons[stop], "*")

    def test_cug_is_serine_only_in_table_12(self):
        self.assertEqual(get_codon_table(1)["CTG"], "L")
        self.assertEqual(get_codon_table(12)["CTG"], "S")

    def test_tga_is_tryptophan_in_table_4(self):
        self.assertEqual(get_codon_table(1)["TGA"], "*")
        self.assertEqual(get_codon_table(4)["TGA"], "W")

    def test_taa_is_glutamine_in_table_6(self):
        self.assertEqual(get_codon_table(1)["TAA"], "*")
        self.assertEqual(get_codon_table(6)["TAA"], "Q")

    def test_a_table_does_not_change_after_another_was_loaded(self):
        # tables are cached; loading one must not alter another
        get_codon_table(12)
        get_codon_table(4)
        self.assertEqual(get_codon_table(1)["CTG"], "L")
        self.assertEqual(get_codon_table(1)["TGA"], "*")


class TranslateTests(unittest.TestCase):
    def test_standard_code(self):
        self.assertEqual(lib.translate("ATGCTGTAA", "+", 0), "ML*")

    def test_alternative_yeast_code(self):
        self.assertEqual(lib.translate("ATGCTGTAA", "+", 0, transl_table=12), "MS*")

    def test_tga_is_read_through_in_table_4(self):
        self.assertEqual(lib.translate("ATGTGATAA", "+", 0, transl_table=1), "M**")
        self.assertEqual(lib.translate("ATGTGATAA", "+", 0, transl_table=4), "MW*")

    def test_taa_is_glutamine_in_table_6(self):
        self.assertEqual(lib.translate("ATGTAATGA", "+", 0, transl_table=6), "MQ*")

    def test_phase_skips_leading_bases(self):
        self.assertEqual(lib.translate("AATGCTGTAA", "+", 1, transl_table=12), "MS*")
        self.assertEqual(lib.translate("AAATGCTGTAA", "+", 2, transl_table=12), "MS*")

    def test_minus_strand_is_reverse_complemented(self):
        gene = "ATGCTGTAA"
        self.assertEqual(lib.translate(revcomp(gene), "-", 0, transl_table=12), "MS*")
        self.assertEqual(lib.translate(revcomp(gene), -1, 0, transl_table=12), "MS*")

    def test_lowercase_is_translated(self):
        self.assertEqual(lib.translate("atgctgtaa", "+", 0, transl_table=12), "MS*")

    def test_incomplete_last_codon_is_dropped(self):
        self.assertEqual(lib.translate("ATGCT", "+", 0), "M")

    def test_ambiguous_codon_is_x(self):
        self.assertEqual(lib.translate("ATGNNNTAA", "+", 0), "MX*")

    def test_table_given_as_string(self):
        self.assertEqual(lib.translate("ATGCTGTAA", "+", 0, transl_table="12"), "MS*")


class GeneModelFixture(unittest.TestCase):
    """One contig with two genes whose CDS read differently under different tables.

    g1 (+ strand):  ATG CTG AAA TGA GGG TAA   table 1: M L K * G *   table 12: M S K * G *   table 4: M L K W G *
    g2 (- strand):  ATG CTG CCC TGA TAA       table 1: M L P * *     table 4: M L P W *
    Both CDS rows carry phase 0, as GFF3 from a gene finder does.
    """

    CDS1 = "ATGCTGAAATGAGGGTAA"
    CDS2 = "ATGCTGCCCTGATAA"

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        if not hasattr(lib, "log"):
            lib.log = logging.getLogger(__name__)
        rng = random.Random(3)
        flank = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
        f1, f2, f3 = flank(60), flank(60), flank(60)
        genome = f1 + self.CDS1 + f2 + revcomp(self.CDS2) + f3
        self.fasta = os.path.join(self.tmp, "genome.fa")
        with open(self.fasta, "w") as f:
            f.write(">ctg1\n" + genome + "\n")
        s1 = len(f1) + 1
        e1 = s1 + len(self.CDS1) - 1
        s2 = e1 + len(f2) + 1
        e2 = s2 + len(self.CDS2) - 1
        self.gff3 = os.path.join(self.tmp, "genes.gff3")
        with open(self.gff3, "w") as f:
            f.write("##gff-version 3\n")
            for gid, strand, a, b in (("g1", "+", s1, e1), ("g2", "-", s2, e2)):
                f.write("ctg1\tt\tgene\t{}\t{}\t.\t{}\t.\tID={}\n".format(a, b, strand, gid))
                f.write("ctg1\tt\tmRNA\t{}\t{}\t.\t{}\t.\tID={}.t1;Parent={}\n".format(a, b, strand, gid, gid))
                f.write("ctg1\tt\texon\t{}\t{}\t.\t{}\t.\tID={}.t1.e1;Parent={}.t1\n".format(a, b, strand, gid, gid))
                f.write("ctg1\tt\tCDS\t{}\t{}\t.\t{}\t0\tID={}.t1.c1;Parent={}.t1\n".format(a, b, strand, gid, gid))

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def convert(self, table):
        prots = os.path.join(self.tmp, "p{}.fa".format(table))
        tbl = os.path.join(self.tmp, "o{}.tbl".format(table))
        lib.convertgff2tbl(
            self.gff3, "", self.fasta, prots,
            os.path.join(self.tmp, "t{}.fa".format(table)), tbl, transl_table=table,
        )
        proteins = {}
        with open(prots) as f:
            name = None
            for line in f:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    proteins[name] = ""
                elif name:
                    proteins[name] += line
        with open(tbl) as f:
            return proteins, f.read()

    def genes(self, table):
        genes = {}
        lib.gff2dict(self.gff3, self.fasta, genes, transl_table=table)
        return genes


class GeneModelTranslationTests(GeneModelFixture):
    def test_table_4_reads_in_frame_tga_as_tryptophan(self):
        genes = self.genes(4)
        self.assertEqual(genes["g1"]["protein"][0], "MLKWG*")
        self.assertEqual(genes["g2"]["protein"][0], "MLPW*")

    def test_table_4_models_are_complete(self):
        for gene in self.genes(4).values():
            self.assertFalse(gene["partialStart"][0])
            self.assertFalse(gene["partialStop"][0])
            self.assertEqual(gene["codon_start"][0], 1)

    def test_table_12_translates_cug_as_serine_on_both_strands(self):
        genes = self.genes(12)
        for gene in genes.values():
            self.assertEqual(gene["protein"][0][:2], "MS")

    def test_table_1_keeps_cug_as_leucine(self):
        for gene in self.genes(1).values():
            self.assertEqual(gene["protein"][0][:2], "ML")

    def test_table_1_does_not_read_through_in_frame_tga(self):
        # same models under the standard code: the TGA is a stop, so the protein is cut
        for gene in self.genes(1).values():
            self.assertNotIn("W", gene["protein"][0])
            self.assertIn("*", gene["protein"][0][:-1])

    def test_alternative_table_gives_different_protein(self):
        self.assertNotEqual(self.genes(1)["g1"]["protein"][0], self.genes(4)["g1"]["protein"][0])


class TblQualifierTests(GeneModelFixture):
    def qualifiers(self, tbl):
        return [line.strip() for line in tbl.splitlines() if line.strip().startswith("transl_table")]

    def test_standard_table_writes_no_transl_table(self):
        _, tbl = self.convert(1)
        self.assertEqual(self.qualifiers(tbl), [])

    def test_alternative_table_writes_one_qualifier_per_cds(self):
        _, tbl = self.convert(4)
        self.assertEqual(self.qualifiers(tbl), ["transl_table\t4", "transl_table\t4"])

    def test_table_12_qualifier(self):
        _, tbl = self.convert(12)
        self.assertEqual(set(self.qualifiers(tbl)), {"transl_table\t12"})

    def test_protein_file_matches_the_table(self):
        proteins, _ = self.convert(4)
        self.assertEqual(proteins["g1.t1"], "MLKWG")
        self.assertEqual(proteins["g2.t1"], "MLPW")

    def test_complete_models_have_no_partial_markers_under_table_4(self):
        _, tbl = self.convert(4)
        gene_lines = [l for l in tbl.splitlines() if l.endswith("\tgene")]
        self.assertEqual(len(gene_lines), 2)
        for line in gene_lines:
            self.assertNotIn("<", line)
            self.assertNotIn(">", line)


if __name__ == "__main__":
    unittest.main()
