import networkx as nx
import pandas as pd

from atomigraph import AtomiGraph


def chain_graph(n=30, broken=()):
    # united-atom PE chain 1..n, optionally without the given bonds
    g = nx.Graph()
    for i in range(1, n + 1):
        g.add_node(i, type=2 if i in (1, n) else 1, element="X")
    g.add_edges_from((i, i + 1) for i in range(1, n) if (i, i + 1) not in broken)
    return g


def run_on_graphs(tmp_path, monkeypatch, graphs, cutoff):
    monkeypatch.chdir(tmp_path)
    ag = AtomiGraph(infile="synthetic", rxn_bond_cutoff=cutoff)
    ag.frames = pd.DataFrame({"frame": list(range(len(graphs))),
                              "timestep": [1000 * i for i in range(len(graphs))],
                              "graph": graphs})
    ag.find_reactions()
    return ag


def test_close_scissions_are_merged(tmp_path, monkeypatch):
    # scissions 8-9 and 11-12 (2 bonds apart): separate cores, but the
    # environments {7..10} and {10..13} overlap for cutoff 1
    ag = run_on_graphs(tmp_path, monkeypatch,
                       [chain_graph(), chain_graph(broken={(8, 9), (11, 12)})], cutoff=1)
    assert len(ag.rxns) == 1
    assert ag.rxns["atoms_env"].iloc[0] == list(range(7, 14))
    assert len(ag.rxns["edges_before"].iloc[0]) == 2


def test_merge_is_transitive(tmp_path, monkeypatch):
    # scissions 5-6, 8-9, 11-12 with cutoff 2: environments A={3..8}, B={6..11}, C={9..14};
    # A overlaps B, B overlaps C, A and C are disjoint -> still a single reaction
    ag = run_on_graphs(tmp_path, monkeypatch,
                       [chain_graph(), chain_graph(broken={(5, 6), (8, 9), (11, 12)})], cutoff=2)
    assert len(ag.rxns) == 1
    assert ag.rxns["atoms_env"].iloc[0] == list(range(3, 15))
    assert len(ag.rxns["edges_before"].iloc[0]) == 3


def test_no_reactions_gives_empty_dataframe(tmp_path, monkeypatch):
    ag = run_on_graphs(tmp_path, monkeypatch, [chain_graph(), chain_graph()], cutoff=1)
    assert ag.rxns.empty
    assert {"timestep", "rxnID", "rxnCount", "edges_before"} <= set(ag.rxns.columns)

    ag.find_reactions()     # a second search must not fail on the empty result
    assert ag.rxns.empty
