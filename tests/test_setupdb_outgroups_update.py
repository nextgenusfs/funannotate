"""Regression test: `funannotate setup -u` always re-downloaded the BUSCO outgroups.

outgroupsDB stores its metadata as info['busco_outgroups'] but called
check4newDB('outgroups', info), which looks up info['outgroups']. That key never
exists, so check4newDB reported the database missing and forced a download on
every update, without comparing md5s.
"""
import importlib
import os
import tempfile
import types
import unittest
from unittest import mock


setupDB = importlib.import_module("funannotate.setupDB")


class OutgroupsUpdateTests(unittest.TestCase):
    def _run_update(self, remote_md5):
        existing = ("outgroups", "/x/outgroups", "1.0", "2020-01-01", 6, "cafef00d")
        info = {"busco_outgroups": existing}
        downloads = []

        def fake_download(url, name, wget=False):
            downloads.append(url)
            open(name, "w").close()

        with tempfile.TemporaryDirectory() as tmp:
            os.mkdir(os.path.join(tmp, "outgroups"))

            def fake_untar(*_a, **_k):                  # tar -zxf recreates outgroups/
                os.makedirs(os.path.join(tmp, "outgroups"), exist_ok=True)
                return 0

            args = types.SimpleNamespace(update=True, wget=False, force=False)
            with mock.patch.object(setupDB, "FUNDB", tmp, create=True), \
                    mock.patch.object(setupDB, "DBURL", {"outgroups": "http://example/outgroups"}, create=True), \
                    mock.patch.object(setupDB, "today", "2026-10-01", create=True), \
                    mock.patch.object(setupDB.lib, "log", mock.MagicMock(), create=True), \
                    mock.patch.object(setupDB, "calcmd5remote", lambda *_a, **_k: remote_md5), \
                    mock.patch.object(setupDB, "download", fake_download), \
                    mock.patch.object(setupDB, "calcmd5", lambda *_a, **_k: remote_md5), \
                    mock.patch.object(setupDB.subprocess, "call", fake_untar):
                setupDB.outgroupsDB(info, force=False, args=args)
        return info, downloads

    def test_current_outgroups_not_redownloaded(self):
        info, downloads = self._run_update("cafef00d")    # remote md5 == stored md5
        self.assertEqual(downloads, [])
        self.assertEqual(info["busco_outgroups"][5], "cafef00d")

    def test_changed_outgroups_redownloaded(self):
        info, downloads = self._run_update("0ddba11")     # remote md5 differs
        self.assertEqual(downloads, ["http://example/outgroups"])
        self.assertEqual(info["busco_outgroups"][5], "0ddba11")


if __name__ == "__main__":
    unittest.main()
