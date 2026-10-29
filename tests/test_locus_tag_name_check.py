"""#1074: --name with curly quotes ended up in the BED file and bedtools failed at the end of predict."""
from funannotate import library as lib


def test_valid_names():
    for name in ("FUN", "Al-Tr-26", "ABC_", "Sneo58285", "a.b"):
        assert lib.invalid_locus_tag_chars(name) == []


def test_curly_quotes_spaces_and_punctuation_are_reported_once():
    assert lib.invalid_locus_tag_chars("“Sneo58285”") == ["“", "”"]
    assert lib.invalid_locus_tag_chars("AB CD") == [" "]
    assert lib.invalid_locus_tag_chars('a"b"c') == ['"']
    assert lib.invalid_locus_tag_chars("café") == ["é"]
