import os
import tempfile
import unittest
import warnings

from funannotate import library as lib


class Gb2outputDeprecationTests(unittest.TestCase):
    def test_call_warns_that_it_is_deprecated(self):
        with tempfile.TemporaryDirectory() as tmp:
            empty = os.path.join(tmp, "empty.gbk")
            open(empty, "w").close()
            outs = [os.path.join(tmp, n) for n in ("p.fa", "t.fa", "s.fa")]
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                lib.gb2output(empty, *outs)
        messages = [str(w.message) for w in caught if issubclass(w.category, DeprecationWarning)]
        self.assertTrue(any("gb2output is deprecated" in m for m in messages), messages)


if __name__ == "__main__":
    unittest.main()
