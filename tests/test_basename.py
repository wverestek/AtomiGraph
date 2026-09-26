import pytest

from atomigraph import AtomiGraph


@pytest.mark.parametrize("infile, expected", [
    ("bonds.reaxff.dump", "bonds.reaxff"),
    ("/some/dir/bonds.reaxff.dump.gz", "bonds.reaxff"),
    ("pe_chain.*.data", "pe_chain"),                            # quoted glob pattern
    (["pe_chain.0.data", "pe_chain.1000.data"], "pe_chain"),     # glob expanded by the shell
    (["bonds.1000.dump", "bonds.1500.dump"], "bonds"),           # cut back to a separator
    (["bonds.reaxff.dump"], "bonds.reaxff"),
    ("*.dump", "AtomiGraph"),
    ("", "AtomiGraph"),
])
def test_derived_basename(infile, expected):
    assert AtomiGraph(infile=infile).basename == expected


def test_explicit_basename_wins():
    assert AtomiGraph(infile="pe_chain.*.data", basename="out/run1").basename == "out/run1"
