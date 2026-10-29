"""dbCAN: the readme next to the V11 files says "release 14.0"; the recorded version follows the V<n> folder."""
from funannotate import setupDB

URL = "https://pro.unl.edu/dbCAN2/download_file.php?file=Databases/V11/dbCAN-HMMdb-V11.txt"


def test_readme_of_a_newer_release_does_not_relabel_the_files():
    head = ["# dbCAN HMMdb release 14.0\n", "# 08/26/2025\n"]
    assert setupDB.dbcan_version_info(URL, head, "2026-10-28") == ("11.0", "2026-10-28")


def test_matching_readme_keeps_its_date():
    head = ["# dbCAN HMMdb release 11.0\n", "# 08/09/2022\n"]
    assert setupDB.dbcan_version_info(URL, head, "2026-10-28") == ("11.0", "2022-08-09")


def test_link_without_version_folder_uses_readme():
    head = ["# dbCAN HMMdb release 12.0\n", "# 08/09/2023\n"]
    assert setupDB.dbcan_version_info("http://x/dbCAN.txt", head, "t") == ("12.0", "2023-08-09")
