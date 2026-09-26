# AtomiGraph

**AtomiGraph** (from **Atomi**stic **Graph**) is a Python library for transforming bonded topologies from molecular dynamics (MD) simulations into graph representations using NetworkX. Its primary focus is the graph-based analysis of molecular connectivity, particularly reaction detection and polymer network topology.

It can be used to analyze bonded network structures in detail or to identify changes in the covalent bond network (i.e., reactions) over time.

This started as a side project during my PhD at IMWF, University of Stuttgart. While it is still under development and some features are not yet fully implemented (e.g., automated reaction kinetics), it already provides a useful framework for topology and reaction analysis in molecular dynamics simulations.

**Disclaimer:**  
Provided *as is*, without any warranty. Use at your own risk — but feel free to use, modify, and build on it. Contributions are very welcome.

---

## What AtomiGraph does

- Parses **LAMMPS data files** and **ReaxFF bond topology dumps**
- Transforms bonded topology into **NetworkX graph representations**
- Enables direct application of **NetworkX algorithms and graph-theoretical analyses**
- Supports removal of nodes by **atom type or pattern** (e.g., for cleanup or coarse-graining)
- Identifies **reaction events** by tracking connectivity changes between timesteps
- Filters out **reversible reactions** within a defined time window  
  (e.g. A + B → C followed by C → A + B)
- Generates **before/after visualizations** for detected reactions

---

## What AtomiGraph does *not* (yet)

- Does **not** consider atomic positions or geometry
- Does **not** include non-covalent interactions (e.g. hydrogen bonds, ionic interactions)
- Does **not** generate SMILES / SMARTS or other cheminformatics outputs
- Does **not** write out coarse-grained MD configurations (yet)
- No automated extraction of reaction rates (yet)

---

## Prerequisites

Python 3 and the following modules:

- default packages: sys, os, random
- numpy
- pandas
- matplotlib
- networkx

---

## Installation

```bash
pip install /path/to/AtomiGraph             # package and the AtomiGraph command
pip install -e "/path/to/AtomiGraph[dev]"   # editable, for development (incl. pytest)
```

Without installing: add the base folder to `PYTHONPATH` (no `AtomiGraph` command then).

---

## Usage

```python3
import atomigraph as ag

topo = ag.AtomiGraph(
    infile="bonds.reaxff.dump",
    atom_type_map="1:C,2:H,3:H,4:O,5:O,6:O,7:O,8:O"
)

topo.read()
topo.find_reactions()           # reactions in topo.rxns (pandas DataFrame)
topo.write_reactions()          # <basename>_rxnIDs.dat, one line per reaction
topo.plot_reactions("png")      # plots in folder <basename>
```

The output basename is derived from the input (here `bonds.reaxff`) unless `basename=...` is given.
For another reaction DataFrame, e.g. after `filter_transient_reactions`, use the functions
`ag.write_reactions(df, basename)` and `ag.plot_reactions(df, basename, outformat)`.

---

## Examples and tests

The examples in `examples/` double as validation cases; each README lists the expected result:

- `01_PE_chain_basic`: synthetic united-atom PE chain (LAMMPS data files): bond scission, formation and flip
- `02_PEEK_one_reaction`: small ReaxFF trajectory with a single bond formation
- `03_PEEK_multiple_files`: same trajectory split into one file per frame (glob input)
- `04_PEEK_many_reactions`: large ReaxFF trajectory (1001 frames), regression values

```Bash
pip install -e ".[dev]"
pytest              # fast tests
pytest -m slow      # large example (04)
```

---

## Command line usage:

After `pip install .` the command `AtomiGraph` is available (or use `python3 -m atomigraph.cli`):

```Bash
AtomiGraph -i 'bonds.reaxff.*.dump' -a 1:C,2:H,3:H,4:O,5:O,6:O,7:O,8:O -c 1 -b rxn_plots
AtomiGraph -i 'pe_chain.*.data' -f lammps_data --plot-format png
AtomiGraph --help
```

---
## Frame sampling
Frame comparison is controlled by these `AtomiGraph` arguments:
- `checkframe`: compare frame i with frame i - checkframe (default 1)
- `stepframe`: frames between evaluations (default 1)
- `stabiframes`: the last n frames are not evaluated (default 0)

Example: `checkframe=1, stepframe=5` compares 0 vs 1, 5 vs 6, 10 vs 11, ...

---

## Visualization

Color coding can be specified for each atom type in the dictionary "type2color" in utils.py.

Default colors follow Jmol conventions:
[Jmol element color convention](https://jmol.sourceforge.net/jscolors/).

---

## Citing

No formal publication yet.

If you use AtomiGraph in academic work, please cite:

> Wolfgang Verestek, AtomiGraph (GitHub repository)

Once a peer-reviewed publication becomes available, please cite the paper instead.

---

## Contributing

Contributions are always welcome, whether it's bug reports, feature requests, documentation improvements, or code contributions.

If you have ideas for new analysis tools or workflows, feel free to open an issue or start a discussion.

---

## License

AtomiGraph is free for non-commercial and academic use.

Commercial or production use requires explicit permission from the author.

Contact: Wolfgang Verestek

