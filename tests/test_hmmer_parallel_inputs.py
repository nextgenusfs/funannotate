"""aux_scripts/hmmer_parallel.py searches only the protein chunks (*.fa).

annotate runs hmmer_parallel.py twice on the same split_prots folder: Pfam
first, then dbCAN. The Pfam step leaves chunk_<n>.pfam.txt there, and the
script passed every file in the folder to hmmscan, so the dbCAN step printed
"error: hmmscan (dbCAN) failed (exit 1)" once per chunk (16 times on a 9,326-gene
genome). A re-run Pfam step would also search the dbCAN outputs.
"""
import importlib.util
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, "funannotate", "aux_scripts", "hmmer_parallel.py")


def load_script():
    spec = importlib.util.spec_from_file_location("hmmer_parallel", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # must not parse sys.argv on import
    return mod


def test_only_fa_chunks_in_natural_order(tmp_path):
    names = [
        "chunk_1.fa", "chunk_2.fa", "chunk_10.fa",
        "chunk_1.pfam.txt", "chunk_2.pfam.txt", "chunk_1.dbcan.txt",
        "dbCAN.txt", "dbCAN.filtered.txt", "pfam.txt",
    ]
    for n in names:
        (tmp_path / n).write_text(">x\nM\n" if n.endswith(".fa") else "# text\n")
    (tmp_path / "sub.fa").mkdir()  # a directory is not a chunk
    got = [os.path.basename(p) for p in load_script().split_prot_files(str(tmp_path))]
    assert got == ["chunk_1.fa", "chunk_2.fa", "chunk_10.fa"]


def test_empty_folder_gives_no_chunks(tmp_path):
    (tmp_path / "chunk_1.pfam.txt").write_text("# text\n")
    assert load_script().split_prot_files(str(tmp_path)) == []
