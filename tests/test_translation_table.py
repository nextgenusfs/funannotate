import unittest
from unittest import mock

from funannotate import funannotate as cli
from funannotate.genetic_codes import get_codon_table


class TranslationTableTests(unittest.TestCase):
    @staticmethod
    def _translate(seq, table):
        codon_table = get_codon_table(table)
        protein = []
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            if len(codon) == 3:
                protein.append(codon_table.get(codon, "X"))
        return "".join(protein)

    def test_standard_table_translation(self):
        self.assertEqual(self._translate("ATGCTGTAA", table=1), "ML*")

    def test_alternative_yeast_translation(self):
        self.assertEqual(self._translate("ATGCTGTAA", table=12), "MS*")


class CliExitCodeTests(unittest.TestCase):
    def test_unknown_subcommand_exits_nonzero(self):
        with mock.patch("sys.argv", ["funannotate", "not-a-command"]):
            with self.assertRaises(SystemExit) as cm:
                cli.main()
        self.assertEqual(cm.exception.code, 1)

    def test_version_exits_zero(self):
        with mock.patch("sys.argv", ["funannotate", "--version"]):
            with self.assertRaises(SystemExit) as cm:
                cli.main()
        self.assertEqual(cm.exception.code, 0)


if __name__ == "__main__":
    unittest.main()
