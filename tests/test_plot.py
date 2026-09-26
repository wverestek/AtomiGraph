from pathlib import Path

import pytest

import atomigraph as ag

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "02_PEEK_one_reaction" / "bonds.reaxff.dump"


@pytest.fixture
def rxns(tmp_path, monkeypatch):
    # find_reactions and the plots write relative to cwd; keep them out of the repo
    monkeypatch.chdir(tmp_path)
    net = ag.AtomiGraph(infile=str(EXAMPLE), atom_type_map="1:C,2:H,3:H,4:O,5:O,6:O,7:O,8:O")
    net.read()
    net.find_reactions()
    return net.rxns


def test_plot_reactions_writes_one_file_per_reaction(tmp_path, rxns):
    ag.plot_reactions(rxns, basename="rxn_plots", outformat="png")

    assert [p.name for p in (tmp_path / "rxn_plots").iterdir()] == \
        ["Reaction0000_timestep2000_Type0_Count1.png"]


def test_plot_rxns_is_deprecated_alias(tmp_path, rxns):
    with pytest.deprecated_call():
        ag.plot_rxns(rxns, basename="rxn_plots", outformat="png")

    assert len(list((tmp_path / "rxn_plots").iterdir())) == 1
