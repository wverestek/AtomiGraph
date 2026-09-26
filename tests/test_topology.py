import networkx as nx
import pandas as pd

from atomigraph import get_degrees, find_minimum_cycle_basis


def frames(*graphs):
    return pd.DataFrame({"frame": range(len(graphs)), "timestep": range(len(graphs)), "graph": graphs})


def test_get_degrees():
    g = nx.path_graph(range(1, 6))                                   # 1-2-3-4-5
    nx.set_node_attributes(g, {i: 2 if i in (1, 5) else 1 for i in g}, name="type")
    df = frames(g, g.subgraph([1, 2, 3]).copy())

    assert get_degrees(df) == [[1, 2, 2, 2, 1], [1, 2, 1]]
    assert get_degrees(df, target_atoms=(2,)) == [[1, 1], [1]]


def fused_and_single_rings():
    # two fused 6-rings (naphthalene-like, 10 atoms) + a separate 8-ring bridged to it + a tail
    g = nx.cycle_graph(6)
    g.add_edges_from([(0, 10), (10, 11), (11, 12), (12, 13), (13, 1)])
    ring8 = nx.relabel_nodes(nx.cycle_graph(8), lambda i: 20 + i)
    g.add_edges_from(ring8.edges())
    g.add_edges_from([(3, 20), (27, 30), (30, 31)])                  # bridge and tail
    return g


def test_minimum_cycle_basis():
    df = frames(fused_and_single_rings())

    # SSSR: the two 6-rings, not the 10-ring perimeter; min_size filters small rings
    cycles = find_minimum_cycle_basis(df, min_size=3)[0]
    assert sorted(len(c) for c in cycles) == [6, 6, 8]
    assert sorted(len(c) for c in find_minimum_cycle_basis(df, min_size=7)[0]) == [8]


def test_minimum_cycle_basis_skips_large_blocks():
    df = frames(fused_and_single_rings())
    assert sorted(len(c) for c in find_minimum_cycle_basis(df, min_size=3, max_block_size=9)[0]) == [8]
