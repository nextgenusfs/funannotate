"""Regression test for GitHub issue #1171.

`funannotate setup` crashed with `TypeError: cannot unpack non-iterable NoneType
object` when a database data file was present on disk but its metadata was missing
from the in-memory `info` dict (e.g. an interrupted earlier setup left the data
file behind but `funannotate.txt` was absent or incomplete). The per-database
helpers only populated `info[<db>]` inside the download branch, then unconditionally
unpacked `info.get(<db>)` -- which is `None` when that branch is skipped.

The fix makes the download branch also trigger when the metadata entry is missing,
so the entry is always (re)populated before it is unpacked.
"""
import importlib
import os
import tempfile
import types
import unittest
from unittest import mock


setupDB = importlib.import_module("funannotate.setupDB")


class MeropsMissingInfoTests(unittest.TestCase):
    def _run_with_existing_fasta_but_empty_info(self, info):
        """Simulate: data file already on disk, but `info` has no 'merops' entry."""
        with tempfile.TemporaryDirectory() as tmp:
            fasta = os.path.join(tmp, "meropsscan.lib")
            with open(fasta, "w") as fh:                # pre-existing data file
                fh.write(">MER0001 #A01.001\nMKQ\n")

            def fake_download(url, name, wget=False):
                with open(name, "w") as fh:
                    fh.write(">MER0001 #A01.001\nMKQ\n")

            args = types.SimpleNamespace(update=False, wget=False, force=False)
            with mock.patch.object(setupDB, "FUNDB", tmp, create=True), \
                    mock.patch.object(setupDB, "DBURL", {"merops": "http://example/merops"}, create=True), \
                    mock.patch.object(setupDB.lib, "log", mock.MagicMock(), create=True), \
                    mock.patch.object(setupDB, "download", fake_download), \
                    mock.patch.object(setupDB, "calcmd5", lambda *_a, **_k: "deadbeef"), \
                    mock.patch.object(setupDB.lib, "runSubprocess", lambda *_a, **_k: None), \
                    mock.patch.object(setupDB.lib, "countfasta", lambda *_a, **_k: 1):
                setupDB.meropsDB(info, force=False, args=args)
            return info

    def test_missing_metadata_does_not_crash_and_is_repopulated(self):
        info = {}                                       # no 'merops' key -> used to crash
        result = self._run_with_existing_fasta_but_empty_info(info)
        self.assertIn("merops", result)
        self.assertEqual(len(result["merops"]), 6)      # full metadata tuple recovered

    def test_present_metadata_is_left_untouched(self):
        # When metadata is already known and the data file exists, no re-download
        # should happen and the existing record must be preserved.
        existing = ("diamond", "/x/merops.dmnd", "12.5", "2023-01-19", 42, "cafef00d")
        info = {"merops": existing}
        with tempfile.TemporaryDirectory() as tmp:
            fasta = os.path.join(tmp, "meropsscan.lib")
            open(fasta, "w").close()
            args = types.SimpleNamespace(update=False, wget=False, force=False)

            def fail_download(*_a, **_k):               # must NOT be called
                raise AssertionError("re-download triggered for an up-to-date DB")

            with mock.patch.object(setupDB, "FUNDB", tmp, create=True), \
                    mock.patch.object(setupDB.lib, "log", mock.MagicMock(), create=True), \
                    mock.patch.object(setupDB, "download", fail_download):
                setupDB.meropsDB(info, force=False, args=args)
        self.assertEqual(info["merops"], existing)


if __name__ == "__main__":
    unittest.main()
