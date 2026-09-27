import pytest

import atomigraph as ag
from common import PE_DIR, PEEK_DUMP, PEEK_TYPE_MAP


@pytest.fixture
def rxns(tmp_path, monkeypatch):
    # the plots are written relative to cwd; keep them out of the repo
    monkeypatch.chdir(tmp_path)
    topo = ag.AtomiGraph(infile=str(PEEK_DUMP), atom_type_map=PEEK_TYPE_MAP)
    topo.read()
    topo.find_reactions()
    return topo.rxns


def test_plot_reactions_writes_one_file_per_reaction(tmp_path, rxns):
    ag.plot_reactions(rxns, basename="rxn_plots", outformat="png")

    assert [p.name for p in (tmp_path / "rxn_plots").iterdir()] == \
        ["Reaction0000_timestep2000_Type0_Count1.png"]


def test_plot_rxns_is_deprecated_alias(tmp_path, rxns):
    with pytest.deprecated_call():
        ag.plot_rxns(rxns, basename="rxn_plots", outformat="png")

    assert len(list((tmp_path / "rxn_plots").iterdir())) == 1


def test_color_table_is_valid():
    from matplotlib.colors import is_color_like
    assert all(is_color_like(c) for c in ag.ON2HEX.values())
    assert all(is_color_like(c) for c in ag.ELEM2HEX.values())


def test_plot_unmapped_atom_types(tmp_path, monkeypatch):
    # without atom_type_map every atom gets element "X" and the fallback color
    monkeypatch.chdir(tmp_path)
    infile = PE_DIR / "pe_chain.*.data"
    topo = ag.AtomiGraph(infile=str(infile), informat="lammps_data")
    topo.read()
    topo.find_reactions()

    ag.plot_reactions(topo.rxns, basename="rxn_plots", outformat="png")

    assert len(list((tmp_path / "rxn_plots").iterdir())) == len(topo.rxns)


def test_write_and_plot_methods_use_basename(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    topo = ag.AtomiGraph(infile=str(PEEK_DUMP), atom_type_map=PEEK_TYPE_MAP)
    topo.read()
    topo.find_reactions()

    topo.write_reactions()
    topo.plot_reactions(outformat="png")

    assert sorted(p.name for p in tmp_path.iterdir()) == ["bonds.reaxff", "bonds.reaxff_rxnIDs.dat"]
    assert len(list((tmp_path / "bonds.reaxff").iterdir())) == 1


def test_functions_default_basename(tmp_path, rxns):
    ag.write_reactions(rxns)
    ag.plot_reactions(rxns, outformat="png")

    assert sorted(p.name for p in tmp_path.iterdir()) == ["AtomiGraph", "AtomiGraph_rxnIDs.dat"]


def test_write_reactions(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    topo = ag.AtomiGraph(infile=str(PE_DIR / "pe_chain.*.data"), informat="lammps_data")
    topo.read()
    topo.find_reactions()
    assert list(tmp_path.iterdir()) == []       # find_reactions itself writes nothing

    ag.write_reactions(topo.rxns, "pe")

    lines = (tmp_path / "pe_rxnIDs.dat").read_text().splitlines()
    assert lines[0].startswith("# Timestep\tRxnID\tRxnCount")
    assert len(lines) == 1 + len(topo.rxns)
    # first scission: one molecule before, two fragments after
    assert lines[1].split("\t")[:4] == ["1000", "0", "1", "[[7, 8, 9, 10]]:[[7, 8], [9, 10]]"]
