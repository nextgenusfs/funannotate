"""downloads.json: dbCAN and MIBiG links, and the versions `funannotate setup` records for them.

funannotate 1.9.x reads downloads.json from the master branch at run time, so a
wrong link here changes every new `funannotate setup`.

- dbCAN: the HMM link must be the V14 file (dbCAN HMMdb release 14.0, 875 HMMs).
  The readme link (Databases/V11/readme.txt) describes release 14.0, so with the
  V14 HMMs the recorded version is 14.0 with the readme date. With the V11 HMMs
  (699 HMMs) the old setup recorded "14.0" for V11 models.
- MIBiG: 4.0. setup takes the version from the file name.
"""
import json
import os

from funannotate import setupDB

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
URLS = json.load(open(os.path.join(ROOT, "funannotate", "downloads.json")))["downloads"]
README_V14_HEAD = ["# dbCAN HMMdb release 14.0\n", "# 08/26/2025\n"]


def test_dbcan_link_is_v14():
    assert "/V14/dbCAN-HMMdb-V14.txt" in URLS["dbCAN"]


def test_dbcan_v14_with_readme_records_14_and_readme_date():
    assert setupDB.dbcan_version_info(URLS["dbCAN"], README_V14_HEAD, "2026-11-03") == ("14.0", "2025-08-26")


def test_mibig_link_is_4_0_and_version_parses():
    url = URLS["mibig"]
    assert url.endswith("/mibig_prot_seqs_4.0.fasta")
    # same expression as setupDB.mibigDB
    assert os.path.basename(url).split("_")[-1].replace(".fasta", "") == "4.0"
