from pathlib import Path

import pytest

from atomigraph.reader import read_lammps_data
from common import PE_DIR

PE_DATA = PE_DIR / "pe_chain.0.data"


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


def test_unmapped_atom_types_warn_once(caplog):
    from atomigraph import AtomiGraph
    infile = str(PE_DATA.parent / "pe_chain.*.data")
    topo = AtomiGraph(infile=infile, informat="lammps_data", atom_type_map="1:C")

    with caplog.at_level("WARNING", logger="atomigraph"):
        topo.read()

    warnings = [r.getMessage() for r in caplog.records if r.levelname == "WARNING"]
    assert warnings == ["atom type(s) [2] not in atom_type_map, element set to 'X'"]
    elements = {g.nodes[n]["element"] for g in topo.frames["graph"] for n in g}
    assert elements == {"C", "X"}


def test_timestep_followed_by_units(tmp_path):
    # newer LAMMPS versions write "timestep = N, units = real"
    old = "LAMMPS data file, united-atom PE chain, timestep = 0"
    infile = write_variant(tmp_path, old,
                           "LAMMPS data file via write_data, version 4 Jul 2026, timestep = 250000, units = real")
    [ts], _ = read_lammps_data(infile)
    assert ts == 250000


def test_mixed_file_names_sort_and_warn(tmp_path, caplog):
    from atomigraph.reader import _expand_infiles
    for name in ("10.data", "2.data", "equi.data", "0.data"):
        (tmp_path / name).write_text("")

    with caplog.at_level("WARNING", logger="atomigraph"):
        files = _expand_infiles(str(tmp_path / "*.data"))

    assert [Path(f).name for f in files] == ["0.data", "2.data", "10.data", "equi.data"]
    assert any("equi.data" in r.getMessage() for r in caplog.records if r.levelname == "WARNING")


def test_numbered_series_does_not_warn(tmp_path, caplog):
    from atomigraph.reader import _expand_infiles
    for name in ("bonds.1000.dump", "bonds.0.dump", "bonds.200.dump"):
        (tmp_path / name).write_text("")

    with caplog.at_level("WARNING", logger="atomigraph"):
        files = _expand_infiles(str(tmp_path / "*.dump"))

    assert [Path(f).name for f in files] == ["bonds.0.dump", "bonds.200.dump", "bonds.1000.dump"]
    assert not [r for r in caplog.records if r.levelname == "WARNING"]
