# -*- coding: utf-8 -*-
"""NCBI translation table lookups.

Provides codon -> amino-acid mappings, stop-codon sets, and start-codon
sets keyed by NCBI genetic code id
(https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).

Built on top of Bio.Data.CodonTable and cached per table id so the dict
is constructed once per table and reused.
"""

from Bio.Data import CodonTable

_TABLE_CACHE = {}
_STOP_CACHE = {}
_START_CACHE = {}


def get_codon_table(table=1):
    """Return dict mapping uppercase codon -> amino acid (stop codons -> '*')."""
    table = int(table)
    if table in _TABLE_CACHE:
        return _TABLE_CACHE[table]
    ncbi = CodonTable.unambiguous_dna_by_id[table]
    d = dict(ncbi.forward_table)
    for stop in ncbi.stop_codons:
        d[stop] = "*"
    _TABLE_CACHE[table] = d
    return d


def get_stop_codons(table=1):
    """Return set of stop codons (uppercase) for the given NCBI table."""
    table = int(table)
    if table in _STOP_CACHE:
        return _STOP_CACHE[table]
    ncbi = CodonTable.unambiguous_dna_by_id[table]
    s = set(ncbi.stop_codons)
    _STOP_CACHE[table] = s
    return s


def get_start_codons(table=1):
    """Return set of start codons (uppercase) for the given NCBI table."""
    table = int(table)
    if table in _START_CACHE:
        return _START_CACHE[table]
    ncbi = CodonTable.unambiguous_dna_by_id[table]
    s = set(ncbi.start_codons)
    _START_CACHE[table] = s
    return s


def is_valid_table(table):
    """return True if table is a valid NCBI translation table id."""
    if isinstance(table, float) and not table.is_integer():
        return False  # int(1.5) would silently pass as table 1
    try:
        return int(table) in CodonTable.unambiguous_dna_by_id
    except (TypeError, ValueError):
        return False


def valid_tables():
    """Sorted list of supported NCBI translation table ids."""
    return sorted(CodonTable.unambiguous_dna_by_id.keys())


def check_tables(table, mtable=None):
    """Validate --table / --mtable and warn about stop-reassigning tables.

    Returns (errors, warnings), lists of message strings. Any error means the run
    should stop; warnings are printed and the run continues. mtable may be None.
    """
    errors, warnings = [], []
    if not is_valid_table(table):
        errors.append("--table {} is not a valid NCBI translation table id".format(table))
    if mtable is not None and not is_valid_table(mtable):
        errors.append("--mtable {} is not a valid NCBI translation table id".format(mtable))
    if not errors and int(table) != 1:
        standard_stops = get_stop_codons(1)
        table_stops = get_stop_codons(int(table))
        if standard_stops != table_stops:
            warnings.append(
                "--table {} reassigns stop codons (table stops {} vs standard {}). "
                "SNAP/GlimmerHMM/CodingQuarry cannot use this table; their gene calls "
                "may be truncated. Augustus and GeneMark accept the table.".format(
                    table, sorted(table_stops), sorted(standard_stops)
                )
            )
    return errors, warnings
