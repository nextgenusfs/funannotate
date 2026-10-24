"""User-supplied FASTA with a blank or comment line before the first record.

Biopython 1.85+ `"fasta"` raises ValueError on such files (1.79 accepted them); funannotate reads user FASTA with
`"fasta-pearson"`, which accepts them. Genomes that went through `funannotate clean` or `sort` are unaffected.
"""
import random

from test_genetic_code_tables import GeneModelFixture

from funannotate import library as lib


class LeadingLinesTests(GeneModelFixture):
    def with_prefix(self, prefix):
        with open(self.fasta) as f:
            text = f.read()
        with open(self.fasta, "w") as f:
            f.write(prefix + text)

    def test_gff2dict_reads_a_genome_with_a_leading_blank_line(self):
        self.with_prefix("\n")
        self.assertEqual(sorted(self.genes(1)), ["g1", "g2"])

    def test_gff2dict_reads_a_genome_with_a_leading_comment(self):
        self.with_prefix("; made by a tool\n#another\n")
        genes = self.genes(4)
        self.assertEqual(genes["g1"]["protein"][0], "MLKWG*")

    def test_gff2dict_matches_the_clean_genome(self):
        clean = self.genes(1)
        self.with_prefix("\n\n")
        self.assertEqual(self.genes(1)["g1"]["protein"], clean["g1"]["protein"])

    def test_clean_proteins_reads_a_file_with_a_leading_blank_line(self):
        rng = random.Random(5)
        prot = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(80))
        src = self.tmp + "/prot.fa"
        with open(src, "w") as f:
            f.write("\n>sp|P1|X desc\n" + prot + "\n")
        out = self.tmp + "/clean.fa"
        lib.cleanProteins([src], out)
        with open(out) as f:
            text = f.read()
        self.assertIn(prot, text)
