"""A resumed `funannotate train` must give PASA the same transcripts as the first run.

mapTranscripts() makes the combined Trinity + long-read FASTA and is skipped when the
alignment BAM exists. Without resumedTranscripts(), PASA then got the Trinity-only FASTA
(1,442 sequences in one test) while the imported alignments also held the long reads (the
first run used 6,180 sequences).
"""
import os

from funannotate.train import resumedTranscripts


def touch(path, text=">seq\nACGT\n"):
    with open(path, "w") as fh:
        fh.write(text)


def paths(tmp_path):
    trinity = str(tmp_path / "trinity.fasta")
    return trinity, trinity + ".clean"


def test_trinity_and_long_reads_use_the_combined_files(tmp_path):
    trinity, clean = paths(tmp_path)
    touch(trinity)
    touch(clean)
    combined = str(tmp_path / "trinity.long-reads.fasta")
    touch(combined)
    touch(combined + ".clean")
    assert resumedTranscripts(str(tmp_path), trinity, clean) == (combined, combined + ".clean")


def test_trinity_only_keeps_the_trinity_files(tmp_path):
    trinity, clean = paths(tmp_path)
    touch(trinity)
    touch(clean)
    assert resumedTranscripts(str(tmp_path), trinity, clean) == (trinity, clean)


def test_long_reads_without_a_trinity_assembly_use_the_long_read_files(tmp_path):
    trinity, clean = paths(tmp_path)  # neither file exists
    long_only = str(tmp_path / "long-reads.mapped.fasta")
    touch(long_only)
    touch(long_only + ".clean")
    assert resumedTranscripts(str(tmp_path), trinity, clean) == (long_only, long_only + ".clean")


def test_an_unfinished_combined_set_is_not_used(tmp_path):
    # the combined FASTA exists but its cleaned copy does not: do not hand PASA half a result
    trinity, clean = paths(tmp_path)
    touch(trinity)
    touch(clean)
    touch(str(tmp_path / "trinity.long-reads.fasta"))
    assert resumedTranscripts(str(tmp_path), trinity, clean) == (trinity, clean)


def test_long_reads_do_not_replace_an_existing_trinity_assembly(tmp_path):
    trinity, clean = paths(tmp_path)
    touch(trinity)
    touch(clean)
    long_only = str(tmp_path / "long-reads.mapped.fasta")
    touch(long_only)
    touch(long_only + ".clean")
    assert resumedTranscripts(str(tmp_path), trinity, clean) == (trinity, clean)
