import networkx as nx
import pandas as pd
import pytest

from atomigraph import AtomiGraph
from common import PE_DIR


def chain_graph(n=30, broken=()):
    # united-atom PE chain 1..n, optionally without the given bonds
    g = nx.Graph()
    for i in range(1, n + 1):
        g.add_node(i, type=2 if i in (1, n) else 1, element="X")
    g.add_edges_from((i, i + 1) for i in range(1, n) if (i, i + 1) not in broken)
    return g


def run_on_graphs(tmp_path, monkeypatch, graphs, cutoff=1, **kwargs):
    monkeypatch.chdir(tmp_path)
    topo = AtomiGraph(infile="synthetic", rxn_bond_cutoff=cutoff, **kwargs)
    topo.frames = pd.DataFrame({"frame": list(range(len(graphs))),
                                "timestep": [1000 * i for i in range(len(graphs))],
                                "graph": graphs})
    topo.find_reactions()
    return topo


def test_close_scissions_are_merged(tmp_path, monkeypatch):
    # scissions 8-9 and 11-12 (2 bonds apart): separate cores, but the
    # environments {7..10} and {10..13} overlap for cutoff 1
    topo = run_on_graphs(tmp_path, monkeypatch,
                         [chain_graph(), chain_graph(broken={(8, 9), (11, 12)})], cutoff=1)
    assert len(topo.rxns) == 1
    assert topo.rxns["atoms_env"].iloc[0] == list(range(7, 14))
    assert len(topo.rxns["edges_before"].iloc[0]) == 2


def test_merge_is_transitive(tmp_path, monkeypatch):
    # scissions 5-6, 8-9, 11-12 with cutoff 2: environments A={3..8}, B={6..11}, C={9..14};
    # A overlaps B, B overlaps C, A and C are disjoint -> still a single reaction
    topo = run_on_graphs(tmp_path, monkeypatch,
                         [chain_graph(), chain_graph(broken={(5, 6), (8, 9), (11, 12)})], cutoff=2)
    assert len(topo.rxns) == 1
    assert topo.rxns["atoms_env"].iloc[0] == list(range(3, 15))
    assert len(topo.rxns["edges_before"].iloc[0]) == 3


def test_no_reactions_gives_empty_dataframe(tmp_path, monkeypatch):
    topo = run_on_graphs(tmp_path, monkeypatch, [chain_graph(), chain_graph()], cutoff=1)
    assert topo.rxns.empty
    assert {"timestep", "rxnID", "rxnCount", "edges_before"} <= set(topo.rxns.columns)

    topo.find_reactions()     # a second search must not fail on the empty result
    assert topo.rxns.empty



def pe_reactions(topo):
    return sorted((ts, tuple(atoms)) for ts, atoms in zip(topo.rxns["timestep"], topo.rxns["atoms_rxn"]))


def test_find_reactions_twice_does_not_duplicate(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    topo = AtomiGraph(infile=str(PE_DIR / "pe_chain.*.data"), informat="lammps_data")
    topo.read()
    topo.find_reactions()
    topo.find_reactions()
    assert len(topo.rxns) == 5


def test_chunked_read_equals_single_read(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    files = [str(PE_DIR / f"pe_chain.{ts}.data") for ts in (0, 1000, 2000, 3000, 4000)]
    single = AtomiGraph(infile=files, informat="lammps_data")
    single.read()
    single.find_reactions()

    chunked = AtomiGraph(infile=files[:3], informat="lammps_data")
    chunked.read()
    chunked.find_reactions()
    chunked.read(infile=files[3:])
    chunked.find_reactions()

    assert list(chunked.frames["frame"]) == [0, 1, 2, 3, 4]
    assert pe_reactions(chunked) == pe_reactions(single)
    assert list(chunked.rxns["frame"]) == list(single.rxns["frame"])


def test_chunked_read_drops_old_frames(tmp_path, monkeypatch):
    # more than 100 frames: old frames are dropped on the next read, numbering continues
    monkeypatch.chdir(tmp_path)
    intact = (PE_DIR / "pe_chain.0.data").read_text()
    broken = (PE_DIR / "pe_chain.1000.data").read_text().replace("timestep = 1000", "timestep = 0")
    for i in range(105):
        text = intact if i < 103 else broken             # two scissions between frame 102 and 103
        (tmp_path / f"f.{i}.data").write_text(text.replace("timestep = 0", f"timestep = {10 * i}"))
    topo = AtomiGraph(infile=[str(tmp_path / f"f.{i}.data") for i in range(102)], informat="lammps_data")
    topo.read()
    topo.find_reactions()
    assert topo.rxns.empty

    topo.read(infile=[str(tmp_path / f"f.{i}.data") for i in range(102, 105)])
    topo.find_reactions()

    assert list(topo.frames["frame"]) == [101, 102, 103, 104]
    assert list(topo.rxns["frame"]) == [103, 103]
    assert list(topo.rxns["timestep"]) == [1030, 1030]


def scission_series(nframes=12):
    # frame i: 60-atom chain with i scissions at 5-6, 10-11, ..., one new scission per frame
    return [chain_graph(n=60, broken={(5 * k, 5 * k + 1) for k in range(1, i + 1)}) for i in range(nframes)]


@pytest.mark.parametrize("kwargs, frames, per_frame", [
    ({}, list(range(1, 12)), 1),                                  # 0 vs 1, 1 vs 2, ...
    ({"stepframe": 5}, [1, 6, 11], 1),                            # 0 vs 1, 5 vs 6, 10 vs 11
    ({"checkframe": 2}, list(range(2, 12)), 2),                   # 0 vs 2: two scissions
    ({"checkframe": 2, "stepframe": 4}, [2, 6, 10], 2),
    ({"stabiframes": 2}, list(range(1, 10)), 1),                  # last 2 frames not evaluated
])
def test_frame_sampling(tmp_path, monkeypatch, kwargs, frames, per_frame):
    topo = run_on_graphs(tmp_path, monkeypatch, scission_series(), **kwargs)
    assert topo.rxns.groupby("frame").size().to_dict() == {f: per_frame for f in frames}


def test_rxnID_and_rxnCount(tmp_path, monkeypatch):
    # all scissions have the same local environment -> one rxnID, counted up
    topo = run_on_graphs(tmp_path, monkeypatch, scission_series())
    assert list(topo.rxns["rxnID"]) == [0] * 11
    assert list(topo.rxns["rxnCount"]) == list(range(1, 12))

    # a different reaction (recombination 5-6) gets the next rxnID, counts continue across calls
    last = topo.frames["graph"].iloc[-1].copy()
    last.add_edge(5, 6)
    topo.frames.loc[len(topo.frames)] = {"frame": 12, "timestep": 12000, "graph": last}
    topo.find_reactions()
    assert list(topo.rxns["rxnID"]) == [0] * 11 + [1]
    assert list(topo.rxns["rxnCount"]) == list(range(1, 12)) + [1]
