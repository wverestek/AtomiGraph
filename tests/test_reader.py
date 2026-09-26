from pathlib import Path

import pytest

from atomigraph.reader import read_lammps_data

PE_DATA = Path(__file__).resolve().parents[1] / "examples" / "01_PE_chain_basic" / "pe_chain.0.data"


def write_variant(tmp_path, old, new):
    f = tmp_path / "variant.data"
    text = PE_DATA.read_text()
    assert old in text
    f.write_text(text.replace(old, new))
    return str(f)


def test_read_pe_chain():
    [ts], [g] = read_lammps_data(str(PE_DATA))
    assert ts == 0
    assert (g.number_of_nodes(), g.number_of_edges()) == (40, 39)


def test_unknown_atom_style_raises(tmp_path):
    infile = write_variant(tmp_path, "Atoms # full", "Atoms # atomic")
    with pytest.raises(ValueError, match=r"variant\.data.*atom style.*atomic"):
        read_lammps_data(infile)


def test_unknown_section_raises(tmp_path):
    infile = write_variant(tmp_path, "Masses", "Ellipsoids")
    with pytest.raises(ValueError, match=r"variant\.data.*Ellipsoids"):
        read_lammps_data(infile)
