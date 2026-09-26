import networkx as nx
import pandas as pd
import pytest

from atomigraph import remove_atoms_by_type, remove_atoms_by_pattern


def make_frames():
    # two frames of a chain CH3-CH2-CH2-CH2-CH3 (types 2-1-1-1-2)
    graphs = []
    for _ in range(2):
        g = nx.path_graph(range(1, 6))
        nx.set_node_attributes(g, {i: 2 if i in (1, 5) else 1 for i in g}, name="type")
        nx.set_node_attributes(g, "X", name="element")
        graphs.append(g)
    return pd.DataFrame({"frame": [0, 1], "timestep": [0, 1000], "graph": graphs})


def snapshot(df):
    return [(sorted(g.nodes()), sorted(g.edges())) for g in df["graph"]]


def test_remove_atoms_by_type_keeps_original():
    df = make_frames()
    before = snapshot(df)

    reduced = remove_atoms_by_type(df, target_atoms=(2,))

    assert snapshot(df) == before
    assert all(sorted(g.nodes()) == [2, 3, 4] for g in reduced["graph"])
    assert all(g is not h for g, h in zip(df["graph"], reduced["graph"]))


def test_remove_atoms_by_pattern_keeps_original():
    df = make_frames()
    before = snapshot(df)

    # template CH3-CH2 (atoms 1-2), delete the CH3 end -> both chain ends are removed
    reduced = remove_atoms_by_pattern(df, template_node_ids=[1, 2], delete_node_ids=[1])

    assert snapshot(df) == before
    assert all(sorted(g.nodes()) == [2, 3, 4] for g in reduced["graph"])
    assert all(g is not h for g, h in zip(df["graph"], reduced["graph"]))


def test_remove_atoms_by_pattern_unknown_ids_raise():
    df = make_frames()
    with pytest.raises(ValueError, match=r"\[99\]"):
        remove_atoms_by_pattern(df, template_node_ids=[1, 2, 99], delete_node_ids=[1])
