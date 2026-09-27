import pandas as pd

from atomigraph import filter_transient_reactions


def make_rxns(rows, index=None):
    # rows: (frame, hash_before, hash_after, atoms_env)
    return pd.DataFrame([{"frame": f, "rxn_hash_before": hb, "rxn_hash_after": ha, "atoms_env": env}
                         for f, hb, ha, env in rows], index=index)


def test_net_reaction_survives():
    # A->B, B->A, A->B: the first pair cancels, the net reaction A->B remains
    df = make_rxns([(1, "A", "B", [1, 2]), (2, "B", "A", [1, 2]), (3, "A", "B", [1, 2])])

    kept, removed = filter_transient_reactions(df, nframes=5)

    assert list(kept["frame"]) == [3]
    assert list(removed["frame"]) == [1, 2]
    assert len(df) == 3


def test_label_index_and_other_atoms():
    # non-positional index labels; reverse on different atoms and outside the window are kept
    df = make_rxns([(1, "A", "B", [1, 2]),
                    (2, "B", "A", [3, 4]),     # other atoms
                    (3, "C", "D", [5, 6]),
                    (4, "B", "A", [1, 2]),     # reverses frame 1
                    (20, "D", "C", [5, 6])],   # outside nframes window of frame 3
                   index=[10, 11, 12, 13, 14])

    kept, removed = filter_transient_reactions(df, nframes=5)

    assert list(removed.index) == [10, 13]
    assert list(kept.index) == [11, 12, 14]


def test_nothing_removed():
    df = make_rxns([(1, "A", "B", [1, 2]), (2, "C", "D", [1, 2])])

    kept, removed = filter_transient_reactions(df, nframes=5)

    assert list(kept["frame"]) == [1, 2]
    assert removed.empty
