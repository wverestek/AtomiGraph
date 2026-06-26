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

Add the AtomiGraph base folder to your `PYTHONPATH`, e.g.:

```bash
export PYTHONPATH=$PYTHONPATH:/path/to/AtomiGraph
```

---

## Usage

```python3
import atomigraph as ag

net = ag.AtomiGraph(
    infile="bonds.reaxff.dump",
    atom_type_map="1:C,2:H,3:H,4:O,5:O,6:O,7:O,8:O"
)

net.read()
net.find_rxns()
net.plot_rxns()
```
---

## Command line usage:

```Bash
tbd
```

---
## Frame/Reaction sampling
Frame comparison can be controlled via:
- startstep: starting frame index
- checkstep: compare frame i up to i + checkstep
- framestep: increment between evaluations

Example:
startstep = 0
checkstep = 1
framestep = 5

-> compares:
0 vs 1, 5 vs 6, 10 vs 11, ...

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

AtomiGraph is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0).

For licensing questions, including potential commercial or alternative licensing arrangements, please contact the author.
