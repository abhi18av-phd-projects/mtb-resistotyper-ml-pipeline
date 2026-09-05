"""Classify a GARC mutation feature by its consequence.

Extracted so that the class FILTER and the class EXPERIMENT share one definition.
They previously would have had two, and a filter that means something slightly
different by "synonymous" than the experiment that justified dropping it is a
quiet way to invalidate both.

Import-light on purpose: a filter should not pull in a gradient-boosting library
to parse a string.
"""

from __future__ import annotations

import re

_SINGLE = re.compile(r"^([A-Za-z])([0-9-]+)([A-Za-z*!])$")

CLASSES = ("synonymous", "nonsynonymous", "promoter_nt", "indel",
           "minor_indel", "nonsense_Z", "other")


def mut_class(col: str) -> str:
    """raw__<gene>_<garc> -> one of CLASSES."""
    mut = col.removeprefix("raw__").partition("_")[2]
    if "minorindel" in mut:
        return "minor_indel"
    if "_ins_" in mut or "_del_" in mut or "indel" in mut:
        return "indel"
    m = _SINGLE.match(mut)
    if m:
        a, _, b = m.groups()
        # A lower-case reference or alternate base is a nucleotide call, which in
        # GARC means a promoter/non-coding position rather than an amino acid.
        if a.islower() or b.islower():
            return "promoter_nt"
        if b in "Zz*!":
            return "nonsense_Z"
        return "synonymous" if a == b else "nonsynonymous"
    return "other"
