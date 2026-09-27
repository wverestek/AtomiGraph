#!/usr/bin/env python3

import numpy as np
import pandas as pd
import networkx as nx
import atomigraph as ag

topo = ag.AtomiGraph(infile="*.data",informat="lammps_data",
                    atom_type_map="1:C,2:C,3:C,4:C,5:H,6:H,7:H,8:N,9:O,10:O,11:O")
topo.read()

def myprint(df:pd.core.frame.DataFrame, mystr:str = "", header:str="graph"):
    for idx,(df_idx,frame) in enumerate(df.iterrows()):
        nxg = frame[header]
        n1 = len(nxg.nodes)
        print(f"{mystr} Index {idx}: {n1} nodes")

myprint(topo.frames)

#####################
# remove atom types #
#####################
# Hydrogens and Oxygen of epoxy group -> linearized
topo.frames = ag.remove_atoms_by_type(topo.frames, target_atoms=(5,6,7,9,11))
myprint(topo.frames,mystr='[H,O]')

############################
# remove isomorph patterns #
############################
# DETDA 3 Carbons from aromatic ring, 5 Carbons dangling -> lineraized
DETDApattern = [65, 67, 70, 72, 542,   63, 64, 66, 68, 69, 71, 204, 545]
DETDAdel = [63, 64, 66, 68, 69, 71, 204, 545]
topo.frames = ag.remove_atoms_by_pattern(topo.frames,DETDApattern,DETDAdel)
# EPON 2 Carbons from aromatic rings -> linearized
EPONpattern = [169, 174, 176, 124, 177, 56,   47, 49]
EPONdel = [47, 49]
topo.frames = ag.remove_atoms_by_pattern(topo.frames,EPONpattern,EPONdel)



#################################
# identify connceted components #
#################################
# initialize arrays
cc = [[] for _ in range(len(topo.frames["graph"]))]
gc_nnodes = [[] for _ in range(len(topo.frames["graph"]))]
gc_cyclo_number = [[] for _ in range(len(topo.frames["graph"]))]
gc_dia = [[] for _ in range(len(topo.frames["graph"]))]
gc_avPath = [[] for _ in range(len(topo.frames["graph"]))]
gc_avClus = [[] for _ in range(len(topo.frames["graph"]))]
gc_sigma = [[] for _ in range(len(topo.frames["graph"]))]

gc_active_nodes = [[] for _ in range(len(topo.frames["graph"]))]
gc_active_nnodes = [[] for _ in range(len(topo.frames["graph"]))]
gc_active_cyclo_number = [[] for _ in range(len(topo.frames["graph"]))]
gc_active_dia = [[] for _ in range(len(topo.frames["graph"]))]
gc_active_avPath = [[] for _ in range(len(topo.frames["graph"]))]
cycle_basis_active = [[] for _ in range(len(topo.frames["graph"]))]

gc_core_nodes = [[] for _ in range(len(topo.frames["graph"]))]
gc_core_nnodes = [[] for _ in range(len(topo.frames["graph"]))]
gc_core_cyclo_number = [[] for _ in range(len(topo.frames["graph"]))]
gc_core_dia = [[] for _ in range(len(topo.frames["graph"]))]
gc_core_avPath = [[] for _ in range(len(topo.frames["graph"]))]
gc_core_avClus = [[] for _ in range(len(topo.frames["graph"]))]
gc_core_sigma = [[] for _ in range(len(topo.frames["graph"]))]

connecting_nodes = [[] for _ in range(len(topo.frames["graph"]))]
connecting_nnodes = [[] for _ in range(len(topo.frames["graph"]))]
connecting_dia = [[] for _ in range(len(topo.frames["graph"]))]
connecting_avPath = [[] for _ in range(len(topo.frames["graph"]))]

dangling_nodes = [[] for _ in range(len(topo.frames["graph"]))]
dangling_nnodes = [[] for _ in range(len(topo.frames["graph"]))]
dangling_dia = [[] for _ in range(len(topo.frames["graph"]))]
dangling_avPath = [[] for _ in range(len(topo.frames["graph"]))]

sol_nodes = [[] for _ in range(len(topo.frames["graph"]))]
sol_nnodes = [[] for _ in range(len(topo.frames["graph"]))]
sol_cyclo_number = [[] for _ in range(len(topo.frames["graph"]))]
sol_dia = [[] for _ in range(len(topo.frames["graph"]))]
sol_avPath = [[] for _ in range(len(topo.frames["graph"]))]

cycle_basis_cores = [[] for _ in range(len(topo.frames["graph"]))]
cycle_basis_connecting = [[] for _ in range(len(topo.frames["graph"]))]
cycle_basis_dangling = [[] for _ in range(len(topo.frames["graph"]))]
cycle_basis_sol = [[] for _ in range(len(topo.frames["graph"]))]
cc_cyclo = [[] for _ in range(len(topo.frames["graph"]))]

for idx,nxg in enumerate(topo.frames["graph"].iloc[:]):
    print(f"sorting idx {idx} into Giant Component (core and dangling) and Sol-Fraction")
    ##### ##### ##### ##### ##### ##### #####
    #  separate in core connecting and sol  #
    ##### ##### ##### ##### ##### ##### #####
    # sort by connection
    cc[idx] = sorted(list(nx.connected_components(nxg)),reverse=True,key=len)                                    # cc[idx][0] = giant_component
    # acyclic part of GC, dangling
    tmp_dang = set(nx.k_crust(nx.subgraph(nxg,list(cc[idx][0])), k=1).nodes())                                   # k=1 => dangling_chains
    dangling_nodes[idx] = [list(i) for i in sorted(nx.connected_components( nx.subgraph(nxg,tmp_dang)),reverse=True,key=len)]
    # if exist, find binconnected core
    tmp = list(nx.subgraph(nxg,list(cc[idx][0])).nodes() - tmp_dang)                                             # cores and connecting in GC
    tmp_gc_core = set([j for i in nx.biconnected_components(nx.subgraph(nxg,tmp)) for j in i if len(i)>2])       # biconnected cores in GC
    tmp = [list(i) for i in sorted(nx.connected_components(nx.subgraph(nxg,tmp_gc_core)),reverse=True,key=len)]
    gc_core_nodes[idx] = tmp if len(tmp)>0 else [[]]
    # acyclic part of GC, connecting
    tmp_conn = set(nx.subgraph(nxg, cc[idx][0]).nodes()) - tmp_gc_core - tmp_dang                                # giant_component - biconnected_cores - dangling_features
    tmp = [list(i) for i in sorted(nx.connected_components(nx.subgraph(nxg, tmp_conn)),reverse=True,key=len)]
    connecting_nodes[idx] = tmp if len(tmp)>0 else [[]]
    # active backbone of GC
    tmp_gc_active = set().union(tmp_gc_core, tmp_conn)
    tmp = [list(i) for i in sorted(nx.connected_components(nx.subgraph(nxg,tmp_gc_active)),reverse=True,key=len)]
    gc_active_nodes[idx] = tmp if len(tmp)>0 else [[]]
    # remaining part of the network.
    tmp = set(nxg.nodes()) - tmp_gc_core - tmp_conn - tmp_dang
    sol_nodes[idx] = [list(i) for i in sorted(nx.connected_components(nx.subgraph(nxg,tmp)),reverse=True,key=len)]
    ##### ##### #####
    # calc metrics  #
    ##### ##### #####
    # calc number of nodes
    gc_nnodes[idx] = len(nx.subgraph(nxg, cc[idx][0]).nodes())
    gc_active_nnodes[idx] = len(nx.subgraph(nxg, tmp_gc_active).nodes())
    gc_core_nnodes[idx] = [len(nx.subgraph(nxg,i).nodes()) for i in gc_core_nodes[idx]]
    connecting_nnodes[idx] = [len(nx.subgraph(nxg,i).nodes) for i in connecting_nodes[idx]]
    dangling_nnodes[idx] = [len(nx.subgraph(nxg,i).nodes) for i in dangling_nodes[idx]]
    sol_nnodes[idx] = [len(list(nx.subgraph(nxg,i).nodes())) for i in sol_nodes[idx]]
    # cyclomatic number
    gc_cyclo_number[idx] = len(nx.subgraph(nxg,cc[idx][0]).edges())-len(nx.subgraph(nxg,cc[idx][0]).nodes())+1
    gc_active_cyclo_number[idx] = len(nx.subgraph(nxg,tmp_gc_active).edges())-len(nx.subgraph(nxg,tmp_gc_active).nodes())+1
    gc_core_cyclo_number[idx] = [len(nx.subgraph(nxg,i).edges())-len(nx.subgraph(nxg,i).nodes())+1 for i in gc_core_nodes[idx] if len(i)>0]
    sol_cyclo_number[idx] = [len(nx.subgraph(nxg,i).edges())-len(nx.subgraph(nxg,i).nodes())+1 for i in sol_nodes[idx]]
    cc_cyclo[idx] = [len(nx.subgraph(nxg,i).edges())-len(nx.subgraph(nxg,i).nodes())+1 for i in cc[idx]]
    # diameter
    gc_dia[idx] = nx.diameter(nx.subgraph(nxg, cc[idx][0]))
    gc_active_dia[idx] = nx.diameter(nx.subgraph(nxg, tmp_gc_active)) if len(tmp_gc_active)>0 else 0
    gc_core_dia[idx] = [nx.diameter(nx.subgraph(nxg, i)) if len(i)>0 else [] for i in gc_core_nodes[idx]]
    connecting_dia[idx] = [nx.diameter(nx.subgraph(nxg,i)) if len(i)>0 else [] for i in connecting_nodes[idx]]
    dangling_dia[idx] = [nx.diameter(nx.subgraph(nxg,i)) for i in dangling_nodes[idx]]
    sol_dia[idx] = [nx.diameter(nx.subgraph(nxg,i)) for i in sol_nodes[idx]]
    # avg path
    gc_avPath[idx] = nx.average_shortest_path_length(nx.subgraph(nxg, cc[idx][0]))
    gc_active_avPath[idx] = nx.average_shortest_path_length(nx.subgraph(nxg, tmp_gc_active)) if len(tmp_gc_active)>0 else 0
    gc_core_avPath[idx] = [nx.average_shortest_path_length(nx.subgraph(nxg,i)) if len(i)>0 else []  for i in gc_core_nodes[idx]] if len(gc_core_nodes[idx][0]) > 0 else []
    connecting_avPath[idx] = [nx.average_shortest_path_length(nx.subgraph(nxg,i)) if len(i)>0 else [] for i in connecting_nodes[idx]]
    dangling_avPath[idx] = [nx.average_shortest_path_length(nx.subgraph(nxg,i)) for i in dangling_nodes[idx]]
    sol_avPath[idx] = [nx.average_shortest_path_length(nx.subgraph(nxg,i)) for i in sol_nodes[idx]]
    # avClus and Sigma (small world)
    gc_avClus[idx] = nx.average_clustering(nx.subgraph(nxg, cc[idx][0]))
    gc_sigma[idx] = (gc_avClus[idx] / (2.0 * len(nx.subgraph(nxg, cc[idx][0]).edges()) / gc_nnodes[idx] / gc_nnodes[idx])) / (gc_avPath[idx] / (np.log(gc_nnodes[idx]) / np.log(2.0 * len(nx.subgraph(nxg, cc[idx][0]).edges()) / gc_nnodes[idx])))
    gc_core_avClus[idx] = [nx.average_clustering(nx.subgraph(nxg,i)) for i in gc_core_nodes[idx]] if len(gc_core_nodes[idx][0]) > 0 else []
    gc_core_sigma[idx] = [ 
        (C / (k / N)) / (l / (np.log(N) / np.log(k)))
        for C, l, N, k in zip(gc_core_avClus[idx],gc_core_avPath[idx],gc_core_nnodes[idx],
                [2.0 * len(nx.subgraph(nxg, nodes).edges()) / len(nodes) for nodes in gc_core_nodes[idx]])] if len(gc_core_nodes[idx][0]) > 0 else []
    # minimum cycle basis: massively expensive
    #cycle_basis_active[idx] = [nx.minimum_cycle_basis(nx.subgraph(nxg,tmp_gc_active)]
    #cycle_basis_cores[idx] = sorted(nx.minimum_cycle_basis(nx.subgraph(nxg,[j for i in gc_core_nodes[idx] for j in i])),reverse=True,key=len)
    #cycle_basis_connecting_nodes[idx] = sorted(nx.minimum_cycle_basis(nx.subgraph(nxg,[j for i in connecting_nodes[idx] for j in i])),reverse=True,key=len)
    #cycle_basis_sol[idx] = sorted(nx.minimum_cycle_basis(nx.subgraph(nxg,list(sol_nodes[idx]))),reverse=True,key=len)
    # alternative if minimum cycle basis is not neccessary
    cycle_basis_active[idx] = [nx.cycle_basis(nx.subgraph(nxg,tmp_gc_active))] if len(tmp_gc_active)>0 else 0
    cycle_basis_cores[idx] = sorted(nx.cycle_basis(nx.subgraph(nxg,[j for i in gc_core_nodes[idx] for j in i])),reverse=True,key=len)
    cycle_basis_connecting[idx] = sorted(nx.cycle_basis(nx.subgraph(nxg,[j for i in connecting_nodes[idx] for j in i])),reverse=True,key=len)
    cycle_basis_dangling[idx] = sorted(nx.cycle_basis(nx.subgraph(nxg,[j for i in dangling_nodes[idx] for j in i])),reverse=True,key=len)
    cycle_basis_sol[idx] = sorted(nx.cycle_basis(nx.subgraph(nxg,[j for i in sol_nodes[idx] for j in i])),reverse=True,key=len)


#############################################
# node fractions of the network per frame   #
# giant component = cores + connecting + dangling
#############################################
nnodes = np.array([g.number_of_nodes() for g in topo.frames["graph"]], dtype=float)
fractions = pd.DataFrame({
    "timestep":   topo.frames["timestep"].to_numpy(),
    "cores":      100 * np.array([sum(n) for n in gc_core_nnodes]) / nnodes,
    "connecting": 100 * np.array([sum(n) for n in connecting_nnodes]) / nnodes,
    "dangling":   100 * np.array([sum(n) for n in dangling_nnodes]) / nnodes,
    "sol":        100 * np.array([sum(n) for n in sol_nnodes]) / nnodes,
})
fractions.to_csv("fractions.csv", index=False, float_format="%.2f")
print(fractions.round(1).to_string(index=False))

import matplotlib.pyplot as plt
parts = ["cores", "connecting", "dangling", "sol"]           # stacked bottom -> top
colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
x = fractions["timestep"] / 1e6
fig, ax = plt.subplots(figsize=(7, 4))
ax.stackplot(x, [fractions[p] for p in parts], labels=parts, colors=colors,
             edgecolor="white", linewidth=1.0)
# direct labels at the right edge for bands thick enough to hold them
bottom = 0.0
for p in parts:
    value = fractions[p].iloc[-1]
    if value > 6:
        ax.text(x.iloc[-1] * 0.98, bottom + value / 2, p, ha="right", va="center", fontsize=9, color="#1a1a19")
    bottom += value
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
ax.set(xlim=(x.min(), x.max()), ylim=(0, 100), xlabel="Timestep [10$^6$]", ylabel="Nodes [%]",
       title="Epoxy network fractions during curing")
ax.grid(axis="y", color="#e0e0e0", linewidth=0.8)
ax.set_axisbelow(True)
ax.spines[["top", "right"]].set_visible(False)
fig.savefig("fractions.png", dpi=200, bbox_inches="tight")

# all results, e.g. for further analysis (needs pickle)
#import pickle
#all_data = [DETDAdel, DETDApattern, EPONdel, EPONpattern, cc, cc_cyclo, connecting_avPath, connecting_dia, connecting_nnodes, connecting_nodes, cycle_basis_active, cycle_basis_connecting, cycle_basis_cores, cycle_basis_dangling, cycle_basis_sol, dangling_avPath, dangling_dia, dangling_nnodes, dangling_nodes, gc_active_avPath, gc_active_cyclo_number, gc_active_dia, gc_active_nnodes, gc_active_nodes, gc_avClus, gc_avPath, gc_core_avClus, gc_core_avPath, gc_core_cyclo_number, gc_core_dia, gc_core_nnodes, gc_core_nodes, gc_core_sigma, gc_cyclo_number, gc_dia, gc_nnodes, gc_sigma, sol_avPath, sol_cyclo_number, sol_dia, sol_nnodes, sol_nodes, topo]
#with open('all_data.pkl', 'wb') as file: pickle.dump(all_data,file)
