# The examples double as validation cases; expected values are documented in each example README.

import pytest

from atomigraph import AtomiGraph
from common import EXAMPLES, PEEK_TYPE_MAP


def run_example(tmp_path, monkeypatch, infile, **kwargs):
    # run in tmp_path so that no output can end up in the repo
    monkeypatch.chdir(tmp_path)
    topo = AtomiGraph(infile=str(EXAMPLES / infile), **kwargs)
    topo.read()
    topo.find_reactions()
    return topo


@pytest.mark.parametrize("cutoff, expected", [
    (0, {1000: 2, 2000: 1, 3000: 1, 4000: 2}),
    (1, {1000: 2, 2000: 1, 3000: 1, 4000: 1}),
    (2, {1000: 2, 2000: 1, 3000: 1, 4000: 1}),
])
def test_01_pe_chain_basic(tmp_path, monkeypatch, cutoff, expected):
    topo = run_example(tmp_path, monkeypatch, "01_PE_chain_basic/pe_chain.*.data",
                       informat="lammps_data", rxn_bond_cutoff=cutoff)

    assert list(topo.frames["timestep"]) == [0, 1000, 2000, 3000, 4000]
    assert topo.rxns.groupby("timestep").size().to_dict() == expected

    by_ts = topo.rxns.groupby("timestep")
    assert sorted(map(sorted, by_ts.get_group(1000)["atoms_rxn"])) == [[8, 9], [21, 22]]
    assert list(by_ts.get_group(2000)["edges_after"]) == [[{8, 22}]]
    rxn_3000 = by_ts.get_group(3000).iloc[0]
    assert (rxn_3000["edges_before"], rxn_3000["edges_after"]) == ([{30, 31}], [{15, 31}])
    assert sorted(a for atoms in by_ts.get_group(4000)["atoms_rxn"] for a in atoms) == [12, 33, 34, 36, 37]


@pytest.mark.parametrize("infile", ["02_PEEK_one_reaction/bonds.reaxff.dump",
                                    "03_PEEK_multiple_files/bonds.reaxff.*.dump"])
@pytest.mark.parametrize("cutoff", [0, 1, 2])
def test_02_03_peek_one_reaction(tmp_path, monkeypatch, infile, cutoff):
    topo = run_example(tmp_path, monkeypatch, infile,
                       atom_type_map=PEEK_TYPE_MAP, rxn_bond_cutoff=cutoff)

    assert list(topo.frames["timestep"]) == [0, 1000, 2000]
    assert len(topo.rxns) == 1
    rxn = topo.rxns.iloc[0]
    assert rxn["timestep"] == 2000
    assert (rxn["edges_before"], rxn["edges_after"]) == ([], [{703, 719}])


@pytest.mark.parametrize("cutoff, n_rxns, n_unique", [(0, 212, 41), (1, 206, 78)])
def test_04_peek_many_reactions(tmp_path, monkeypatch, cutoff, n_rxns, n_unique):
    # regression values recorded with the current code, not independently validated
    topo = run_example(tmp_path, monkeypatch, "04_PEEK_many_reactions/bonds.reaxff.dump",
                       atom_type_map=PEEK_TYPE_MAP, rxn_bond_cutoff=cutoff)

    assert len(topo.frames) == 101
    assert len(topo.rxns) == n_rxns
    assert topo.rxns["rxnID"].nunique() == n_unique


@pytest.mark.parametrize("cutoff", [0, 1])
def test_05_epoxy_network_reactions(tmp_path, monkeypatch, cutoff):
    topo = run_example(tmp_path, monkeypatch, "05_epoxy_network/*.data", informat="lammps_data",
                       rxn_bond_cutoff=cutoff)

    assert len(topo.frames) == 19
    assert topo.frames["timestep"].is_monotonic_increasing
    assert topo.rxns.groupby("timestep").size().tolist() == [1] * 18


def test_05_epoxy_network_script(tmp_path, monkeypatch):
    # run the analysis script on a copy, so that its output stays out of the repo
    import runpy, shutil
    src = EXAMPLES / "05_epoxy_network"
    for f in [*src.glob("*.data"), src / "run_ag.py"]:
        shutil.copy(f, tmp_path)
    monkeypatch.chdir(tmp_path)

    ns = runpy.run_path("run_ag.py")

    assert [g.number_of_nodes() for g in ns["topo"].frames["graph"]] == [195] * 19
    fractions = ns["fractions"]
    parts = ["cores", "connecting", "dangling", "sol"]
    assert (fractions[parts].sum(axis=1).round(6) == 100).all()
    assert fractions.iloc[-1][parts].round(1).tolist() == [71.3, 0.0, 8.7, 20.0]
    assert (tmp_path / "fractions.png").exists()
    assert (tmp_path / "fractions.csv").read_text() == (src / "fractions_reference.csv").read_text()
