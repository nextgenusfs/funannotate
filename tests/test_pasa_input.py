import os
import tempfile
import unittest

from funannotate.pasa_input import (
    contained_transcripts,
    read_alignment_exons,
    top_isoforms,
    trinity_gene,
    write_fasta_without,
    write_table_without,
    write_gff3_without,
)

G = "Trinity_GG_1_c0_g1"


def aln(*exons, ch="chr1", st="+"):
    return [(ch, st, sorted(exons))]


class ContainedTests(unittest.TestCase):
    def test_gene_key(self):
        self.assertEqual(trinity_gene(G + "_i12"), G)
        self.assertEqual(trinity_gene("asmbl_7"), "asmbl_7")

    def test_fragment_with_shared_introns_is_dropped(self):
        a = {G + "_i1": aln((100, 200), (300, 400), (500, 600)),
             G + "_i2": aln((150, 200), (300, 400))}
        self.assertEqual(contained_transcripts(a), {G + "_i2"})

    def test_alternative_splice_site_is_kept(self):
        # i2 uses a different acceptor (310), so it adds a junction
        a = {G + "_i1": aln((100, 200), (300, 400), (500, 600)),
             G + "_i2": aln((150, 200), (310, 400))}
        self.assertEqual(contained_transcripts(a), set())

    def test_exon_skipping_isoform_is_kept(self):
        a = {G + "_i1": aln((100, 200), (300, 400), (500, 600)),
             G + "_i2": aln((100, 200), (500, 600))}
        self.assertEqual(contained_transcripts(a), set())

    def test_alternative_last_exon_extending_beyond_is_kept(self):
        # i2 shares the first intron, but its last exon runs past exon 2 of i1
        a = {G + "_i1": aln((100, 200), (300, 400), (500, 600)),
             G + "_i2": aln((150, 200), (300, 450))}
        self.assertEqual(contained_transcripts(a), set())

    def test_single_exon_inside_an_exon_is_dropped(self):
        a = {G + "_i1": aln((100, 200), (300, 400)),
             G + "_i2": aln((310, 390))}
        self.assertEqual(contained_transcripts(a), {G + "_i2"})

    def test_single_exon_spanning_an_intron_is_kept(self):
        # retained intron: new information
        a = {G + "_i1": aln((100, 200), (300, 400)),
             G + "_i2": aln((150, 350))}
        self.assertEqual(contained_transcripts(a), set())

    def test_identical_chains_keep_one(self):
        a = {G + "_i1": aln((100, 200), (300, 400)),
             G + "_i2": aln((100, 200), (300, 400))}
        self.assertEqual(contained_transcripts(a), {G + "_i2"})

    def test_other_gene_strand_or_chrom_not_compared(self):
        a = {G + "_i1": aln((100, 200), (300, 400), (500, 600)),
             "Trinity_GG_2_c0_g1_i1": aln((150, 200), (300, 400)),
             G + "_i3": aln((150, 200), (300, 400), st="-"),
             G + "_i4": aln((150, 200), (300, 400), ch="chr2")}
        self.assertEqual(contained_transcripts(a), set())

    def test_multi_locus_transcript_never_dropped(self):
        a = {G + "_i1": aln((100, 200), (300, 400), (500, 600)),
             G + "_i2": [("chr1", "+", [(150, 200), (300, 400)]),
                         ("chr3", "+", [(10, 90)])]}
        self.assertEqual(contained_transcripts(a), set())


class IntronsModeTests(unittest.TestCase):
    def test_fragment_with_longer_terminal_exon_dropped(self):
        # strict keeps it (end runs past exon 2 of i1); introns mode drops it
        a = {G + "_i1": aln((100, 200), (300, 400), (500, 600)),
             G + "_i2": aln((50, 200), (300, 450))}
        self.assertEqual(contained_transcripts(a, "strict"), set())
        self.assertEqual(contained_transcripts(a, "introns"), {G + "_i2"})

    def test_new_junction_still_kept(self):
        a = {G + "_i1": aln((100, 200), (300, 400), (500, 600)),
             G + "_i2": aln((150, 200), (310, 400))}
        self.assertEqual(contained_transcripts(a, "introns"), set())

    def test_identical_chain_keeps_widest(self):
        a = {G + "_i1": aln((100, 200), (300, 400)),
             G + "_i2": aln((20, 200), (300, 480))}
        self.assertEqual(contained_transcripts(a, "introns"), {G + "_i1"})

    def test_unknown_mode(self):
        with self.assertRaises(ValueError):
            contained_transcripts({}, "loose")


class TopIsoformTests(unittest.TestCase):
    def test_zero_keeps_all(self):
        ids = [G + "_i%d" % i for i in range(5)]
        self.assertEqual(top_isoforms(ids, {}, {}, 0), set())

    def test_keeps_most_abundant_then_longest(self):
        ids = [G + "_i1", G + "_i2", G + "_i3", "Trinity_GG_2_c0_g1_i1"]
        tpm = {G + "_i1": 5.0, G + "_i2": 50.0, G + "_i3": 5.0}
        lens = {G + "_i1": 900, G + "_i2": 100, G + "_i3": 1200}
        self.assertEqual(top_isoforms(ids, tpm, lens, 2), {G + "_i1"})


class FileFilterTests(unittest.TestCase):
    def test_read_and_filter_files(self):
        with tempfile.TemporaryDirectory() as d:
            gff = os.path.join(d, "a.gff3")
            with open(gff, "w") as fh:
                fh.write("##gff-version 3\n")
                for tid, s, e in ((G + "_i1", 100, 200), (G + "_i1", 300, 400),
                                  (G + "_i2", 310, 390)):
                    fh.write("chr1\tgenome\tcDNA_match\t%d\t%d\t99\t+\t.\tID=%s;Target=%s 1 10\n"
                             % (s, e, tid, tid))
            a = read_alignment_exons(gff)
            self.assertEqual(a[G + "_i1"], [("chr1", "+", [(100, 200), (300, 400)])])
            drop = contained_transcripts(a)
            out = os.path.join(d, "b.gff3")
            write_gff3_without(gff, drop, out)
            with open(out) as fh:
                body = fh.read()
            self.assertNotIn(G + "_i2", body)
            self.assertIn("##gff-version 3", body)
            fa = os.path.join(d, "t.fa")
            with open(fa, "w") as fh:
                fh.write(">%s len=5\nACGTA\n>%s\nAC\nGT\n" % (G + "_i1", G + "_i2"))
            fo = os.path.join(d, "o.fa")
            self.assertEqual(write_fasta_without(fa, drop, fo), 1)
            with open(fo) as fh:
                self.assertEqual(fh.read(), ">%s len=5\nACGTA\n" % (G + "_i1"))
            cln = os.path.join(d, "t.fa.cln")
            with open(cln, "w") as fh:
                fh.write("%s\t 0.00\t    1\t 5\t 5\t\t\n%s\t 0.00\t    1\t 4\t 4\t\t\n" % (G + "_i1", G + "_i2"))
            co = os.path.join(d, "o.fa.cln")
            write_table_without(cln, drop, co)
            with open(co) as fh:
                self.assertEqual(fh.read().split("\t")[0], G + "_i1")


if __name__ == "__main__":
    unittest.main()
