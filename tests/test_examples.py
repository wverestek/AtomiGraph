# The examples double as validation cases; expected values are documented in each example README.
from pathlib import Path

import pytest

from atomigraph import AtomiGraph

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
PEEK_TYPE_MAP = "1:C,2:H,3:H,4:O,5:O,6:O,7:O,8:O"


def run_example(tmp_path, monkeypatch, infile, **kwargs):
    # find_reactions writes <basename>_rxnIDs.dat relative to cwd; keep it out of the repo
    monkeypatch.chdir(tmp_path)
    ag = AtomiGraph(infile=str(EXAMPLES / infile), **kwargs)
    ag.read()
    ag.find_reactions()
    return ag


@pytest.mark.parametrize("cutoff, expected", [
    (0, {1000: 2, 2000: 1, 3000: 1, 4000: 2}),
    (1, {1000: 2, 2000: 1, 3000: 1, 4000: 1}),
    (2, {1000: 2, 2000: 1, 3000: 1, 4000: 1}),
])
def test_01_pe_chain_basic(tmp_path, monkeypatch, cutoff, expected):
    ag = run_example(tmp_path, monkeypatch, "01_PE_chain_basic/pe_chain.*.data",
                     informat="lammps_data", rxn_bond_cutoff=cutoff)

    assert list(ag.frames["timestep"]) == [0, 1000, 2000, 3000, 4000]
    assert ag.rxns.groupby("timestep").size().to_dict() == expected

    by_ts = ag.rxns.groupby("timestep")
    assert sorted(map(sorted, by_ts.get_group(1000)["atoms_rxn"])) == [[8, 9], [21, 22]]
    assert list(by_ts.get_group(2000)["edges_after"]) == [[{8, 22}]]
    rxn_3000 = by_ts.get_group(3000).iloc[0]
    assert (rxn_3000["edges_before"], rxn_3000["edges_after"]) == ([{30, 31}], [{15, 31}])
    assert sorted(a for atoms in by_ts.get_group(4000)["atoms_rxn"] for a in atoms) == [12, 33, 34, 36, 37]


@pytest.mark.parametrize("infile", ["02_PEEK_one_reaction/bonds.reaxff.dump",
                                    "03_PEEK_multiple_files/bonds.reaxff.*.dump"])
@pytest.mark.parametrize("cutoff", [0, 1, 2])
def test_02_03_peek_one_reaction(tmp_path, monkeypatch, infile, cutoff):
    ag = run_example(tmp_path, monkeypatch, infile,
                     atom_type_map=PEEK_TYPE_MAP, rxn_bond_cutoff=cutoff)

    assert list(ag.frames["timestep"]) == [0, 1000, 2000]
    assert len(ag.rxns) == 1
    rxn = ag.rxns.iloc[0]
    assert rxn["timestep"] == 2000
    assert (rxn["edges_before"], rxn["edges_after"]) == ([], [{703, 719}])


@pytest.mark.slow
@pytest.mark.parametrize("cutoff, n_rxns, n_unique", [(0, 859, 38), (1, 856, 123)])
def test_04_peek_many_reactions(tmp_path, monkeypatch, cutoff, n_rxns, n_unique):
    # regression values recorded with the current code, not independently validated
    ag = run_example(tmp_path, monkeypatch, "04_PEEK_many_reactions/bonds.reaxff.dump",
                     atom_type_map=PEEK_TYPE_MAP, rxn_bond_cutoff=cutoff)

    assert len(ag.frames) == 1001
    assert len(ag.rxns) == n_rxns
    assert ag.rxns["rxnID"].nunique() == n_unique
