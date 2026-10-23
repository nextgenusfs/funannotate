import io

from funannotate.check import ErrorCounter


def test_error_counter_counts_error_lines_and_passes_text_through():
    out = io.StringIO()
    counter = ErrorCounter(out)
    counter.write("\tblat found\n")
    counter.write("\tERROR: pslDnaFiler not installed\n")
    counter.write("   ERROR: PASAHOME not set\n")
    counter.flush()
    assert counter.count == 2
    assert out.getvalue().count("\n") == 3


def test_error_counter_is_zero_without_errors():
    counter = ErrorCounter(io.StringIO())
    counter.write("NOTE: Rust PASA not found\n")
    assert counter.count == 0
