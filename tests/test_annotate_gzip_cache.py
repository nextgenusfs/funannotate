import importlib
import importlib.util
import os
import tempfile
import unittest
from unittest import mock


if not hasattr(importlib.util, "module_for_loader"):
    def module_for_loader(func):
        return func


    importlib.util.module_for_loader = module_for_loader


annotate = importlib.import_module("funannotate.annotate")
annotate.FUNDB = "/tmp/funannotate-db"


class _FakeHSP(object):
    aln_span = 100
    ident_num = 95


class _FakeHit(object):
    hsps = [_FakeHSP()]
    description = "Valid protein OS=Example species GN=ABC1"
    id = "sp|P12345|ABC1"


class _FakeResult(object):
    hits = [_FakeHit()]
    seq_len = 100
    id = "gene1"


class SwissProtBlastGzipCacheTests(unittest.TestCase):
    def test_zero_byte_gzip_cache_is_ignored(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gzip_cache = os.path.join(tmpdir, "uniprot.xml.gz")
            open(gzip_cache, "w").close()

            def _write_plain_xml(*args, **kwargs):
                with open(os.path.join(tmpdir, "uniprot.xml"), "w") as handle:
                    handle.write("<xml />\n")

            with mock.patch.object(annotate.lib, "log", mock.Mock(), create=True), \
                mock.patch.object(annotate.lib, "runSubprocess", side_effect=_write_plain_xml) as run_subprocess, \
                mock.patch.object(annotate.SearchIO, "parse", return_value=[]):
                annotate.SwissProtBlast("proteins.fa", 1, 1e-5, tmpdir, {}, diamond=True)

            self.assertEqual(run_subprocess.call_count, 1)

    def test_invalid_gzip_cache_is_regenerated(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gzip_cache = os.path.join(tmpdir, "uniprot.xml.gz")
            with open(gzip_cache, "wb") as handle:
                handle.write(b"not-a-valid-gzip")

            def _write_plain_xml(*args, **kwargs):
                with open(os.path.join(tmpdir, "uniprot.xml"), "w") as handle:
                    handle.write("<xml />\n")

            def _parse_results(handle, fmt):
                self.assertEqual(fmt, "blast-xml")
                if handle.name.endswith(".gz"):
                    raise OSError("truncated gzip cache")
                return [_FakeResult()]

            gene_dict = {}
            with mock.patch.object(annotate.lib, "log", mock.Mock(), create=True), \
                mock.patch.object(annotate.lib, "runSubprocess", side_effect=_write_plain_xml) as run_subprocess, \
                mock.patch.object(annotate.SearchIO, "parse", side_effect=_parse_results):
                annotate.SwissProtBlast("proteins.fa", 1, 1e-5, tmpdir, gene_dict, diamond=True)

            self.assertEqual(run_subprocess.call_count, 1)
            self.assertIn("gene1", gene_dict)
            self.assertEqual(gene_dict["gene1"][0]["name"], "ABC1")
            self.assertTrue(os.path.isfile(gzip_cache))
            self.assertFalse(os.path.exists(os.path.join(tmpdir, "uniprot.xml")))


if __name__ == "__main__":
    unittest.main()
