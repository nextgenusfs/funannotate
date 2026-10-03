"""Regression test: annotate --genbank without --table crashed.

annotate.py set the default --table only for predict-folder input. With
--genbank (or --gff) input args.table stayed None and the tbl2asn step
crashed at int(args.table). For --genbank input the table now comes from the
GenBank CDS transl_table qualifier (the code gb2parts writes to the tbl).

Found running `funannotate test -t annotate` in the v1.9.0-rc.5 container.
"""
import os
import tempfile
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from funannotate import library as lib


def write_gbk(path, tables):
    """one record; one CDS per entry in tables (None = no transl_table)"""
    rec = SeqRecord(Seq("ATG" + "GCT" * 20 + "TAA" + "N" * 30), id="contig1",
                    name="contig1", annotations={"molecule_type": "DNA"})
    for i, table in enumerate(tables):
        quals = {"locus_tag": ["G{}".format(i)]}
        if table is not None:
            quals["transl_table"] = [str(table)]
        rec.features.append(SeqFeature(FeatureLocation(0, 66, strand=1),
                                       type="CDS", qualifiers=quals))
    SeqIO.write([rec], path, "genbank")


class GbTranslTableTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.gbk = os.path.join(self.tmp.name, "in.gbk")

    def tearDown(self):
        self.tmp.cleanup()

    def test_declared_table(self):
        write_gbk(self.gbk, [12])
        self.assertEqual(lib.gb_transl_table(self.gbk), 12)

    def test_first_declared_table_wins(self):
        write_gbk(self.gbk, [None, 4, 12])
        self.assertEqual(lib.gb_transl_table(self.gbk), 4)

    def test_no_table_is_none(self):
        write_gbk(self.gbk, [None])
        self.assertIsNone(lib.gb_transl_table(self.gbk))


if __name__ == "__main__":
    unittest.main()
