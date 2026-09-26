#!/usr/bin/env python3
# python -m atomigraph.cli  or  AtomiGraph (entry point)

import sys
import argparse

from atomigraph.core import AtomiGraph, write_reactions, plot_reactions


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="AtomiGraph",
        description="Find and plot changes in bond topology (reactions) in LAMMPS output.",
        epilog="Example: AtomiGraph -i 'bonds.reaxff.*.dump' -a 1:C,2:H,3:H,4:O,5:O,6:O,7:O,8:O",
    )
    parser.add_argument("-i", "--input", required=True, nargs="+",
                        help="input file(s) or quoted glob pattern, e.g. 'bonds.reaxff.*.dump'")
    parser.add_argument("-f", "--format", default="reaxff", choices=["reaxff", "lammps_data"],
                        help="input format (default: reaxff)")
    parser.add_argument("-a", "--atom-map", default="",
                        help="atom type to element mapping, e.g. '1:C,2:H,3:O'")
    parser.add_argument("-b", "--basename", default="",
                        help="output basename: summary <basename>_rxnIDs.dat and plot folder <basename> "
                             "(default: input file name without extension)")
    parser.add_argument("-c", "--cutoff", type=int, default=1,
                        help="rxn_bond_cutoff: bonds around the changed bonds that belong to a reaction (default: 1)")
    parser.add_argument("--checkframe", type=int, default=1, help="frame difference to check (default: 1)")
    parser.add_argument("--stepframe", type=int, default=1, help="frames between evaluations (default: 1)")
    parser.add_argument("--plot-format", default="pdf", choices=["pdf", "png"],
                        help="file format of the reaction plots (default: pdf)")
    parser.add_argument("--no-plot", action="store_true", help="only find reactions, do not plot")
    args = parser.parse_args(argv)

    infile = args.input[0] if len(args.input) == 1 else args.input
    try:
        topo = AtomiGraph(infile=infile, informat=args.format, atom_type_map=args.atom_map,
                          rxn_bond_cutoff=args.cutoff,
                          checkframe=args.checkframe, stepframe=args.stepframe)
        topo.read()
        topo.find_reactions()
    except (OSError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    if topo.rxns.empty:
        print("No reactions found.")
        return 0
    print(f"{len(topo.rxns)} reaction(s) found, {topo.rxns['rxnID'].nunique()} unique.")

    basename = args.basename or topo.basename
    write_reactions(topo.rxns, filename=f"{basename}_rxnIDs.dat")
    if not args.no_plot:
        plot_reactions(topo.rxns, basename=basename, outformat=args.plot_format)
    return 0


if __name__ == "__main__":
    sys.exit(main())
