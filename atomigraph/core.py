# python atomigraph/core.py
#from fileinput import filename
import os, os.path

from networkx.algorithms.operators import union
from networkx.drawing import draw
os.environ.setdefault("MPLBACKEND", "Agg")
import sys, re
import warnings
import random

from typing import TextIO, Union, List

import matplotlib
# Use non-interactive backend to avoid Qt/X11 errors in headless environments
matplotlib.use("Agg")

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from .utils import k_nearest_neighs, convert_str2dict
from .utils import ON2ELEM, ON2HEX, ELEM2HEX, DEFAULT_COLOR
from .reader import read_bonds
from .logger import log, configure_log

configure_log(level="DEBUG", force=True)

#import pdb
#pdb.set_trace()

__all__ = ['DEFAULT_COLOR', 'ELEM2HEX', 'ON2ELEM', 'ON2HEX', 
           'AtomiGraph',
           'renumber_and_count_rxns', 'filter_transient_reactions', 'remove_atoms_by_type', 'remove_atoms_by_pattern', 
           'write_reactions', 'plot_reactions', 'plot_rxns', 'get_degrees', 'find_minimum_cycle_basis']


class AtomiGraph:
    ##############
    # initialize #
    ##############
    def __init__(self, 
                 infile: Union[str, List[str]] = "", informat: str = "reaxff", basename: str = "",
                 atom_type_map: str = "",
                 checkframe: int = 1, stepframe: int = 1, stabiframes: int = 0,
                 hash_by: str = "type", 
                 rxn_bond_cutoff: int = 1, plot_bonds_cutoff: int = 5, seed: int = 42,
                 ring_counter: bool = False, loop_limits:tuple[int,int]=None):
        """
        A class to extract changes in bond topology over time.

        infile : str or list[str]
            A bond information file, a list of files, or glob pattern(s). The
            reader will accept a single filename or multiple files. For
            backward compatibility `self.infile` is the first filename (or "").
        infile : str or list[str]
            A bond information file, a list of files, or glob pattern(s). The
            reader will accept a single filename or multiple files. For
            backward compatibility `self.infile` is the first filename (or "").
        informat : str
            file type of the file containing bond information. 
            reaxff, lammps_data
            Default "reaxff"
        basename : str
            base name for output. If not set the input file name is used as base name.
        startstep : int
            MD time step at start. Default: 0
        stopstep : int
            MD time step at end. Default: sys.maxsize (a very high umber)
        checkframe : int
            Number of frames difference to check for changed bonds. Default: 1
        stepframe : int
            Number of frames before the next evaluation is done. Default: 1
        hash_by : str
            hash for each reaction that is build upon 'element' or 'type'. 
            Hashes allow to identify if a similar reaction has already occured before or not.
        """
        
        log.info(f"Initializing AtomiGraph class object...")
        self.name: str = "AtomiGraph"

        # store input(s): `infiles` may be a single filename (str) or a list of filenames
        # reader.read_bonds accepts either a string, a list or glob pattern(s).
        self.infile: Union[str, List[str]] = infile
        self.informat: str = informat.lower()
        # derive basename based on only or first filename string
        if isinstance(self.infile, str) and len(self.infile) > 0:
            self.basename: str = re.sub(r'(\.(?:gz|txt|dat|data|dump))+$', '', os.path.basename(self.infile))
        elif isinstance(self.infile, list):
            self.basename: str = re.sub(r'(\.(?:gz|txt|dat|data|dump))+$', '', os.path.basename(self.infile[0]))
        
        self.checkframe:int = int(checkframe)   # necessary?
        self.stabiframe:int = int(stabiframes)
        self.stepframe:int = int(stepframe)     # necessary?
        self.maxframe_offset = 0
        
        # pandas DataFrame to store global bond topology
        self.frames:pd.DataFrame = pd.DataFrame(columns=("frame","timestep","graph")) # DataFrame with columns ['timestep','graph']
        
        # pandas DataFrame to store reactions found
        self.rxns:pd.DataFrame = pd.DataFrame(
            columns=("frame","timestep","rxnID","rxnCount",
                     "edges_before","edges_after",
                     "atoms_rxn","atoms_env","atoms_plot",
                     "Gbefore","Gafter",
                     "rxn_hash_before","rxn_hash_after"))                          # DataFrame with reactions found
        
        # build hash based on node attribute 'element' or 'type'
        self.hash_by:str = hash_by                      
        
        self.rxn_bond_cutoff:int = int(rxn_bond_cutoff)         # number of bonds distance that are considered one reaction, default: 0
        self.plot_bonds_cutoff:int = int(plot_bonds_cutoff)     # number of bonds distance to include for plots
        random.seed(a=int(seed))                                # reproducibility of pyplot plots

        # find and count loops in each global frame        
        self.loop_limits = loop_limits

        # expose ON2ELEM/ON2HEX as instance attributes so the rest of the class can use them
        self.on2elem:dict = ON2ELEM
        self.on2hex:dict = ON2HEX
        self.elem2hex:dict = ELEM2HEX
        self.default_color:str = DEFAULT_COLOR
        # atom mapping from type to ordinal number -> user input as str
        self.atom_type_map = atom_type_map
        self.type2on:dict = convert_str2dict(self.atom_type_map)
        self.elem2hex:dict = {self.on2elem[i]:self.on2hex[i] for i in self.on2elem.keys()}

    # read bonds #
    def read(self, infile: Union[str, List[str]] = None, informat: str = None):
        """
        Read bond file and populate instance state:
          - self.ts, self.nxg
          - self.startidx / self.stopidx (based on startstep/stopstep)
          - self.frames (pd.DataFrame) with columns ['timestep','graph']
        Returns the frames DataFrame.
        """
        # choose infile to pass to reader.read_bonds:
        # - if caller provided `infile` use that
        # - else prefer self.infiles (list) if present, otherwise self.infile (string)
        if infile is not None:
            infile_local = infile
        else:
            # use self.infile (string or list)
            infile_local = self.infile

        informat_local = informat or self.informat
        if not infile_local:
            raise ValueError("No infile supplied to read()")

        # read bond file(s)
        ts, nxg = read_bonds(infile_local, informat_local)

        # set element attribute for each node in each graph
        unmapped_types = set()
        for g in nxg:
            for n, data in g.nodes(data=True):
                atom_type = data.get("type", None)
                if atom_type is not None and atom_type in self.type2on:
                    data["element"] = self.type2on[atom_type]
                else:
                    data["element"] = "X"
                    unmapped_types.add(atom_type)
        if unmapped_types:
            log.warning(f"atom type(s) {sorted(unmapped_types, key=lambda t: (t is None, t or 0))} not in atom_type_map, element set to 'X'")
        
        frames_arr = list(range(len(ts)))
        if self.frames.empty:
            # first read, no existing frames, just set frames DataFrame
            self.frames = pd.DataFrame({"frame":frames_arr,"timestep": ts,"graph": nxg})
        else:
            # subsequent read
            # drop frames in case large files are read in multiple calls to avoid 
            # memory issues, but keep stabilize frames for reaction checking
            nframes = len(self.frames["frame"])
            if  nframes > 100 and nframes > self.stabiframe:
                maxframe_keep =  self.frames["frame"].iloc[-1] - self.stabiframe
                idx = np.where(self.frames["frame"].lt(maxframe_keep))[0].tolist()
                _ = self.frames.drop(index=idx, inplace=True)
                # renumber new frames to continue from last frame + 1
                maxframe_old = self.frames["frame"].max() if not self.frames.empty else -1
                frames_arr = [i + maxframe_old + 1 for i in frames_arr]
            # concatenate new frames to existing frames DataFrame
            self.frames = pd.concat([self.frames, 
                                     pd.DataFrame({"frame":frames_arr,"timestep": ts,"graph": nxg})],
                                    ignore_index=True)
        # check timesteps are monotonic increasing
        if not self.frames["timestep"].is_monotonic_increasing:
            log.warning(f"READER: Timesteps are not monotonic increasing!!!")

    # find reactions #
    def find_reactions(self):
        """
        Find reactions, i.e. changes in bond topology, between frames and store them in self.rxns.

        Frame idx is compared with frame idx - checkframe for idx = checkframe, checkframe + stepframe, ...
        up to the last frame minus stabiframes. If no frames have been read yet, self.read() is called.

        For each pair of frames:
        - changed bonds are the bonds present in only one of the two graphs (broken or formed);
        - atoms connected via changed bonds (and bonds between them) form the core of a reaction,
          so a broken and a newly formed bond sharing atoms (bond flip) are one reaction;
        - each core is expanded by rxn_bond_cutoff bonds in both frames; expanded sets that share
          atoms are merged (transitively) into one reaction;
        - each reaction is hashed (Weisfeiler-Lehman, node attribute hash_by) before and after.

        Results are appended to self.rxns (pandas.DataFrame, one row per reaction occurrence), then
        rxnID and rxnCount are renumbered over all stored reactions: reactions with the same
        before:after hash pair share a rxnID (in order of first appearance), rxnCount counts
        their occurrences. Columns:
            frame, timestep            frame index and MD timestep of the "after" frame
            rxnID, rxnCount            reaction type and running count of this type
            edges_before, edges_after  bonds broken / formed, list of {atom_i, atom_j}
            atoms_rxn                  atoms of the changed bonds
            atoms_env                  atoms of the reaction incl. rxn_bond_cutoff environment
            atoms_plot                 atoms_env plus plot_bonds_cutoff environment
            Gbefore, Gafter            subgraphs of atoms_plot before / after
            rxn_hash_before, rxn_hash_after  WL hashes of atoms_env before / after

        If no reaction is found at all, self.rxns is None.
        Use write_reactions(self.rxns) for a text summary and plot_reactions(self.rxns) for plots.

        Example:
            net = AtomiGraph(infile="bonds.reaxff.dump", atom_type_map="1:C,2:H,3:O")
            net.read()
            net.find_reactions()
            net.rxns[["timestep", "rxnID", "edges_before", "edges_after"]]
        """

        if self.frames.empty:
            log.info(f"No bond data read yet, reading now... {self.infile}")
            self.read()
        
        # check for ascending timesteps
        if not self.frames["timestep"].is_monotonic_increasing:
            log.warning(f"Timesteps are not monotonic increasing!!!")

        log.info("Searching reactions...")
        

        start = 0 
        stop =  len(self.frames["frame"]) - self.stabiframe 
        cf = self.checkframe
        fs = self.stepframe

        #self.rxn_id = []
        #self.rxn_count = []

        df_file = pd.DataFrame(columns=self.rxns.columns)
        
        for idx in range(start + cf, stop, fs):
            # dicts for conversion
            node2element = nx.get_node_attributes(self.frames["graph"].iloc[idx], name="element")
            node2type    = nx.get_node_attributes(self.frames["graph"].iloc[idx], name="type")

            before_idx = idx - cf
            after_idx = idx

            Gbefore = self.frames["graph"].iloc[idx-cf]
            Gafter = self.frames["graph"].iloc[idx]

            log.info(f"evaluating frame: {idx} ({self.frames['frame'].iloc[idx]}) timestep: {self.frames['timestep'].iloc[idx]}")
            # find reacting atoms
            [edges_sets_before,edges_sets_after, reaction_sets] = self._find_reacting_atoms_for_two_frames(Gbefore, Gafter)   # list(set(int,),set(int,))
            log.debug(f"reaction_sets: {reaction_sets} with edges_before: {edges_sets_before} and edges_after: {edges_sets_after}")
                
            if len(reaction_sets) > 0:
                df_frame = self._rsets_to_pd(before_idx, after_idx,
                                            reaction_sets, edges_sets_before, edges_sets_after)
                if df_frame is not None and not df_frame.empty:
                    df_file = pd.concat([df_file, df_frame], ignore_index=True)
                else:
                    log.warning("You should not be here. Maybe you have discovered a bug. Please consider reporting with a minimal example")

        
        # adding newly found reactions to self.rxns DataFrame
        # first search
        if self.rxns.empty:
            self.rxns = df_file.copy()
        else:
            log.info(f"appending new search to already stored reactions")
            df_file["frame"] = df_file["frame"] + self.maxframe_offset
            self.rxns = pd.concat([self.rxns,df_file],ignore_index=True)

        #self.df1 = df_file.copy()
        # renumber reactions and count unique reactions
        self.rxns = renumber_and_count_rxns(self.rxns)

    # find reacting atoms for two frames #
    def _find_reacting_atoms_for_two_frames(self,Gbefore:nx.Graph,Gafter:nx.Graph):
        # compare graphs, edges that are not in both graphs => reaction
        #reacting_edges = nx.symmetric_difference(Gbefore,Gafter).edges()
        all_edges_before = nx.difference(Gbefore,Gafter).edges()
        all_edges_after  = nx.difference(Gafter,Gbefore).edges()
        all_reacting_edges = set(all_edges_before).union(set(all_edges_after))
        log.debug(f"identified edges: {all_reacting_edges}")
        # atom IDs that are involved with changed bond connectivity
        all_reacting_atoms = set(i for j in set(all_reacting_edges) for i in j)       # set(int,)
        log.debug(f"identified atoms: {all_reacting_atoms}")

        # subgraph with reacting atoms and edges (before and after) 
        # to find connected components as individual reactions
        Gcombine = Gbefore.subgraph(all_reacting_atoms).copy()              # nx.Graph
        Gcombine.add_edges_from(all_reacting_edges)                         # nx.Graph, add reacting edges 
        reacting_atoms_core_sets = list(nx.connected_components(Gcombine))       # list(set(int,),)
        
        nsets = len(reacting_atoms_core_sets)
        # reactions found
        if nsets > 0:
            
            # include atoms rxn_bond_cutoff bonds away
            if self.rxn_bond_cutoff > 0:
                
                tmpsets = []
                reacting_atoms_sets = []

                # expand individual sets by rxn_bond_cutoff bonds
                log.debug(f"reacting_atoms_sets before expansion: {reacting_atoms_core_sets}")
                for i,iset in enumerate(reacting_atoms_core_sets):
                    reacting_atoms1 = k_nearest_neighs(Gbefore,iset,self.rxn_bond_cutoff)
                    reacting_atoms2 = k_nearest_neighs(Gafter ,iset,self.rxn_bond_cutoff)  
                    tmpsets.append(reacting_atoms1.union(reacting_atoms2))    # combine sets
                log.debug(f"reacting_atoms_sets after expansion: {tmpsets}")

                # more than one set, merge expanded sets that share atoms (transitively):
                # one node per expanded set, edge if two sets intersect, union per connected component
                if len(tmpsets) > 1:
                    log.debug(f"Merging: reacting_atoms_sets before merging: {tmpsets}")
                    Gmerge = nx.Graph()
                    Gmerge.add_nodes_from(range(len(tmpsets)))
                    for i in range(len(tmpsets)):
                        for j in range(i+1, len(tmpsets)):
                            if not tmpsets[i].isdisjoint(tmpsets[j]):
                                Gmerge.add_edge(i, j)
                    for component in nx.connected_components(Gmerge):
                        reacting_atoms_sets.append(set().union(*(tmpsets[i] for i in component)))
                    log.debug(f"Merged: reacting_atoms_sets after merging:   {reacting_atoms_sets}")

                else:
                    # only one set, just use the expanded set
                    reacting_atoms_sets = tmpsets
                    
            # no expansion, just use original sets of reacting atoms and edges
            else:
                reacting_atoms_sets = reacting_atoms_core_sets

            log.debug(f"reacting_atoms_sets: {reacting_atoms_sets}")

            # construct edge sets for each reaction set
            reacting_edges_before_sets = list(set())
            reacting_edges_after_sets = list(set())
            for i,rset in enumerate(reacting_atoms_sets):
                reacting_edges_before_sets.append([set(bond) for bond in all_edges_before if set(bond).issubset(rset)])
                reacting_edges_after_sets.append( [set(bond) for bond in all_edges_after  if set(bond).issubset(rset)])
            log.debug(f"reacting_edges_before: {reacting_edges_before_sets}")
            log.debug(f"reacting_edges_after: {reacting_edges_after_sets}")
        
        # no reactions found, return empty lists
        else:
            reacting_edges_before_sets = []
            reacting_edges_after_sets = []
            reacting_atoms_sets = []

        return reacting_edges_before_sets, reacting_edges_after_sets, reacting_atoms_sets

    # convert reaction sets to pandas DataFrame format for further analysis and plotting #
    def _rsets_to_pd(self,before:int,after:int,reaction_sets,edges_sets_before,edges_sets_after):
        """
        Convert reaction sets and their associated edge changes into 
        a pandas DataFrame format for further analysis and plotting.
        Parameters:
        before (int): Index of the "before" frame.
        after (int): Index of the "after" frame.
        reaction_sets (list of sets): List of sets of atom IDs involved in each reaction.
        edges_sets_before (list of sets): edges that vanish in this reaction.
        edges_sets_after (list of sets): edges that become created in this reaction.
        """
        tmp_list = []
        # before and after sets for each individual reaction
        for i,(rset,edges_before,edges_after) in enumerate(zip(reaction_sets,edges_sets_before,edges_sets_after)):
            log.debug(f"rset: {rset} with edges_before {edges_before} and edges_after {edges_after}")
            log.debug(f"rset: edge atoms {list(set().union(*edges_before).union(*edges_after))}")
            atoms_rxn = list(set().union(*edges_before).union(*edges_after))
            atoms_env = list(rset)
            plot_atoms1 = k_nearest_neighs(self.frames["graph"].iloc[before], atoms_env, self.plot_bonds_cutoff)
            plot_atoms2 = k_nearest_neighs(self.frames["graph"].iloc[after] , atoms_env, self.plot_bonds_cutoff)  
            atoms_plot = list(plot_atoms1.union(plot_atoms2))     # combine sets
            
            # before including plot atoms 
            Gbefore = self.frames["graph"].iloc[before].subgraph(atoms_plot).copy()
            hash_before = nx.weisfeiler_lehman_graph_hash(Gbefore.subgraph(atoms_env), node_attr=self.hash_by)
            # after including plot atoms 
            Gafter = self.frames["graph"].iloc[after].subgraph(atoms_plot).copy()
            hash_after = nx.weisfeiler_lehman_graph_hash(Gafter.subgraph(atoms_env), node_attr=self.hash_by)

            tmp_list.append({"frame":self.frames["frame"].iloc[after],
                             "timestep":self.frames["timestep"].iloc[after],
                             "rxnID":None,"rxnCount":None,
                             "edges_before":edges_before,"edges_after":edges_after,
                             "atoms_rxn":sorted(atoms_rxn),"atoms_env":sorted(atoms_env),"atoms_plot":sorted(atoms_plot),
                             "Gbefore":Gbefore,"Gafter":Gafter,
                             "rxn_hash_before":hash_before,"rxn_hash_after":hash_after})
        return pd.DataFrame(tmp_list)
# End of class



## work on reactions and topology ##
# renumber reactions and count unique reactions #
def renumber_and_count_rxns(df:pd.core.frame.DataFrame=None) -> pd.core.frame.DataFrame:
    """
    Renumber reactions and count unique reactions based on their hashes in pandas DataFrame. 
    This method updates two columns of the DataFrame: 'rxnID' and 'rxnCount'.
    If df is None, operates on self.rxns and updates it in place. Otherwise, operates on the 
    provided DataFrame and returns a new DataFrame with the updated columns.
    """
    # operate on a copy to avoid surprising in-place side effects for caller
    df_work = df.copy(deep=True)

    if df_work is None or df_work.empty:
        log.info("No reactions to renumber and count.")
        return

    # reset index to ensure consistent indexing for reaction ID assignment
    df_work.reset_index(drop=True, inplace=True)

    crxns = len(df_work["rxn_hash_before"]) # total count of reactions
    rxn_hashes = df_work["rxn_hash_before"] + [":"]*crxns + df_work["rxn_hash_after"]
    id2hash = dict(enumerate(pd.unique(rxn_hashes))).items()    # unique rxn_hashes in order of appearance
    hash2id = dict((v,k) for k,v in id2hash)                    # unique rxn_hashes in order of appearance
    nrxns = len(hash2id.keys())                                 # number of individual reactions

    rxn_id = np.array([None] * crxns)
    rxn_count = np.array([0] * crxns)
    counter_arr = np.array([0] * nrxns)

    for idx, h in enumerate(rxn_hashes):
        rxn_id[idx] = hash2id[h]
        counter_arr[hash2id[h]] += 1
        rxn_count[idx] = counter_arr[hash2id[h]]
    
    df_work["rxnID"] = rxn_id
    df_work["rxnCount"] = rxn_count

    return df_work


# remove reactions that reverse in stabilize frames #
def filter_transient_reactions(df:pd.core.frame.DataFrame=None, nframes:int=None)-> list[pd.core.frame.DataFrame,pd.core.frame.DataFrame]:
    """
    Remove reactions that reverse within a certain number of frames (stabiframe) after their 
    occurrence. Based on 'frame', 'rxn_hash_before' and 'rxn_hash_after' in pandas DataFrame.
    Returns a tuple of two DataFrames: (filtered_reactions, removed_reactions) where:
    - filtered_reactions: DataFrame with reactions that do not reverse within stabiframe frames.
    - removed_reactions: DataFrame with reactions that reverse within stabiframe frames.
    """
    # operate on a copy to avoid surprising in-place side effects for caller
    df_work = df.copy()
        
    # each reaction can cancel at most one reverse reaction; consumed rows are skipped,
    # so A->B, B->A, A->B removes the first pair and keeps the net reaction A->B
    rmv_idx = []
    consumed = set()
    atoms_env = df_work["atoms_env"].map(tuple)     # list comparison: convert to tuple
    for idx,row in df_work.iterrows():
        if idx in consumed:
            continue
        current_frame = row["frame"]
        max_frame = current_frame + nframes

        mask_frame = df_work["frame"].gt(current_frame) & df_work["frame"].le(max_frame)
        mask_hash = (df_work["rxn_hash_before"] == row["rxn_hash_after"]) & (df_work["rxn_hash_after"] == row["rxn_hash_before"])
        mask_atoms = atoms_env == tuple(row["atoms_env"])
        mask_free = ~df_work.index.isin(consumed)
        candidates = df_work.index[mask_frame & mask_hash & mask_atoms & mask_free]
        if len(candidates) > 0:
            # earliest reverse reaction, as index label
            rev_idx = df_work.loc[candidates, "frame"].idxmin()
            consumed.update((idx, rev_idx))
            rmv_idx.extend((idx, rev_idx))

    # removed reactions in original order
    df_rmv = df_work.loc[df_work.index.isin(rmv_idx)]
    df_work = df_work.drop(index=rmv_idx)
    if len(rmv_idx) > 0:
        log.info(f"{len(rmv_idx)} Reaction(s) found that reverse within {nframes} frames, removing reactions")

    return df_work, df_rmv


# remove_atoms #
def remove_atoms_by_type(df:pd.core.frame.DataFrame=None, target_atoms:tuple[int|str,...]=None) -> pd.core.frame.DataFrame:
    # copy DataFrame to avoid surprising in-place side effects for caller
    df_work = df.copy(deep=True)
            
    for idx, row in df_work.iterrows():
        # real copy: df.copy(deep=True) does not copy the nx.Graph objects
        nxg = row["graph"].copy()

        if target_atoms is None:
            nodes = list(nxg.nodes())
        else:
            nodes = [node for node, node_data in nxg.nodes(data=True)
                     if node_data.get('type') in target_atoms 
                     or node_data.get('element') in target_atoms]

        nxg.remove_nodes_from(nodes)
        
        # store the modified copy in the DataFrame
        df_work.at[idx, "graph"] = nxg
    
    # return independent DataFrame with modified graphs
    return df_work


# remove atoms by somorph search of subgraph #
def remove_atoms_by_pattern(df:pd.core.frame.DataFrame, template_node_ids:list|set, delete_node_ids:list|set, node_attr:str='type', pattern_from_frame:int=0) -> pd.core.frame.DataFrame:
    """
    Simplifies the graph by matching a template pattern and removing specific 
    nodes. Handles molecular symmetry by filtering unique node sets.
    """
    from networkx.algorithms import isomorphism

    # 1. Sanity Check
    template_set = set(template_node_ids)
    delete_set = set(delete_node_ids)
    if not delete_set.issubset(template_set):
        invalid_ids = list(delete_set - template_set)
        log.error(f"Nodes {invalid_ids} are not in template_node_ids!")
        raise ValueError("Invalid delete_node_ids provided.")
    
    df_work = df.copy(deep=True)

    # 2. Create the template graph
    template = df_work["graph"].iloc[pattern_from_frame].subgraph(template_node_ids).copy()
    log.info(f"Starting topology reduction")
    log.info(f"Template pattern nodes: {list(template.nodes())}")

    nm = nx.isomorphism.categorical_node_match(node_attr, None)

    for idx, (df_idx, frame) in enumerate(df_work.iterrows()):
        # real copy: df.copy(deep=True) does not copy the nx.Graph objects
        nxg = frame["graph"].copy()
        df_work.at[df_idx, "graph"] = nxg
        ncomp_before = nx.number_connected_components(nxg)

        # 3. Setup the GraphMatcher
        gm = nx.isomorphism.GraphMatcher(nxg, template, node_match=nm)

        # 4. Identify UNIQUE matches (Filtering Symmetries)
        nodes_to_remove = set()
        unique_molecule_footprints = set()

        for match in gm.subgraph_isomorphisms_iter():
            # match.keys() found pattern ID, match.values() corresponding template ID
            match_dict = dict(match)
            molecule_footprint = frozenset(match_dict.keys())
        
            if molecule_footprint not in unique_molecule_footprints:
                unique_molecule_footprints.add(molecule_footprint)
                # add key (pattern id) if value (template id) is in delete_node_ids
                #nodes_to_remove.add([k for k,v in match_dict.items() if v in delete_node_ids])
                # delete nodes only for the first isomorphism found per molecule
                inv_match = {v: k for k, v in match.items()}
                for d_id in delete_node_ids:
                    if d_id in inv_match:
                        nodes_to_remove.add(inv_match[d_id])
    
        match_count = len(unique_molecule_footprints)

        if match_count == 0:
            log.warning(f"Index {idx}: No matches found! Check atom types.")
            continue

        # 5. Global removal
        initial_count = nxg.number_of_nodes()
        nxg.remove_nodes_from(nodes_to_remove)
    
        # 6. Connectivity Check
        ncomp_after = nx.number_connected_components(nxg)
        log.info(f"Index {idx}: Pattern found {match_count} times. Reduced nodes from {initial_count} to {nxg.number_of_nodes()}.")

        if ncomp_after > ncomp_before:
            log.warning(f"Index {idx}: Network is no longer fully connected! Fragments: {ncomp_before} vs. {ncomp_after}")
    
    return df_work

# write reactions #
def write_reactions(df:pd.core.frame.DataFrame, filename:str="AtomiGraph_rxnIDs.dat") -> None:
    """
    Write a tab-separated summary with one line per reaction:
    timestep, rxnID, rxnCount, molecules before:after as atom IDs, atom types and elements,
    and the reaction hashes before:after. Molecules are the connected parts of the reaction
    environment (atoms_env) in the frame before and after the reaction.
    """
    header = "# Timestep\tRxnID\tRxnCount\tFromIDs:ToIDs\tFromType:ToType\tFromElem:ToElem\tRxn_hashes"
    with open(filename, "wt") as f:
        f.write(header + "\n")
        if df is None or df.empty:
            log.warning("No reactions found to write.")
            return
        for _, rxn in df.iterrows():
            ids, types, elems = [], [], []
            for G in (rxn["Gbefore"], rxn["Gafter"]):
                mols = sorted(sorted(c) for c in nx.connected_components(G.subgraph(rxn["atoms_env"])))
                ids.append(mols)
                types.append([[G.nodes[a].get("type") for a in m] for m in mols])
                elems.append([[G.nodes[a].get("element") for a in m] for m in mols])
            fields = [rxn["timestep"], rxn["rxnID"], rxn["rxnCount"],
                      f"{ids[0]}:{ids[1]}", f"{types[0]}:{types[1]}", f"{elems[0]}:{elems[1]}",
                      f"{rxn['rxn_hash_before']}:{rxn['rxn_hash_after']}"]
            f.write("\t".join(str(x) for x in fields) + "\n")
    log.info(f"{len(df)} reaction(s) written to {filename}")


# plot reactions #
def plot_reactions(df:pd.core.frame.DataFrame, basename:str="AtomiGraph", outformat:str="pdf") -> None:
    # check if DataFrame is empty
    if df.empty:
        log.warning("No reactions found to plot.")
        return

    outfolder = basename or "AtomiGraph_outdir"
    if not os.path.exists(outfolder):
        os.makedirs(outfolder, exist_ok=True)
    
    digitsIDX = np.max([4, len(str(len(df)))])
    digitsTS =  len(str(df["timestep"].max()))
    digitsID =  len(str(df["rxnID"].max()))
    digitsCount =  len(str(df["rxnCount"].max()))

    for idx, rxn in df.iterrows():
        frame = rxn["frame"]
        timestep = rxn["timestep"]
        rxnID = rxn["rxnID"]
        rxnCount = rxn["rxnCount"]
        hash_before = rxn["rxn_hash_before"]
        hash_after = rxn["rxn_hash_after"]

        Gbefore = rxn["Gbefore"]
        Gafter = rxn["Gafter"]
        #active_edges_before = nx.difference(Gbefore, Gafter).edges()
        #active_edges_after  = nx.difference(Gafter, Gbefore).edges()
        active_edges = nx.symmetric_difference(Gbefore, Gafter).edges()

        # initial positions
        Gprint = Gbefore.copy()
        Gprint.add_edges_from(active_edges)
        pos = nx.kamada_kawai_layout(Gprint)
        pos = nx.spring_layout(Gprint, iterations=200, pos=pos)
        del Gprint

        
        plt.clf()
        # before reaction
        # prepare data
        plt.figure(figsize=(14,6))
        plt.suptitle(f"Index: {idx} Timestep: {timestep} Type: {rxnID}\n {hash_before}:{hash_after}")
        plt.subplot(1, 2, 1)
        plt.title("Before")
        node2elem = nx.get_node_attributes(Gbefore, name="element")
        node2type = nx.get_node_attributes(Gbefore, name="type")
        color = [ELEM2HEX.get(v, DEFAULT_COLOR) for v in
                  nx.get_node_attributes(Gbefore, name="element").values()]
        node_labels = {k:  str(node2elem[k]) + ":" + str(node2type[k]) + "\n" + str(k) for k in Gbefore}
        if len(nx.get_edge_attributes(Gbefore,"bo")) > 0:
            bo = [tmp for tmp in list(nx.get_edge_attributes(Gbefore, "bo").values())]
        else:
            bo = [1.0 for i in Gbefore.edges()]
        pos_before = nx.spring_layout(Gbefore, iterations=75, pos=pos)
        # plot data
        #nx.draw_networkx_edges(Gbefore, edgelist=active_edges_before, alpha=0.6, width=5.0, edge_color="tab:red", pos=pos)
        nx.draw_networkx_edges(Gbefore, edgelist=active_edges, alpha=0.6, width=5.0, edge_color="tab:red", pos=pos_before)
        nx.draw(Gbefore, pos=pos_before, node_color=color, 
                with_labels=True, labels=node_labels, font_size=6,
                node_size=300, edge_color="black", width=bo)
        plt.tight_layout()
        
        # after reaction
        # prepare data
        plt.subplot(1, 2, 2)
        plt.title("After")
        color = [ELEM2HEX.get(v, DEFAULT_COLOR) for v in
                 nx.get_node_attributes(Gafter, name="element").values()]
        node_labels = {k:  str(node2elem[k]) + ":" + str(node2type[k]) + "\n" + str(k) for k in Gafter}
        if len(nx.get_edge_attributes(Gafter,"bo")) > 0:
            bo = [tmp for tmp in list(nx.get_edge_attributes(Gbefore, "bo").values())]
        else:
            bo = [1.0 for i in Gbefore.edges()]
        pos_after = nx.spring_layout(Gafter, iterations=75, pos=pos)
        # plot data
        #nx.draw_networkx_edges(Gafter, edgelist=active_edges_after, alpha=0.6, width=5.0, edge_color="tab:red", pos=pos)
        nx.draw_networkx_edges(Gafter, edgelist=active_edges, alpha=0.6, width=5.0, edge_color="tab:red", pos=pos_after)
        nx.draw(Gafter, pos=pos_after, node_color=color, 
                with_labels=True, labels=node_labels, font_size=6,
                node_size=300, edge_color="black", width=bo)
        plt.tight_layout()
        
        filename = f"Reaction{idx:0{digitsIDX}d}_timestep{timestep:0{digitsTS}d}_Type{rxnID:0{digitsID}d}_Count{rxnCount:0{digitsCount}d}"

        if outformat == "png":
            fileout = filename+".png"
            f_out = os.path.join(outfolder, fileout)
            plt.savefig(f_out,dpi=200)
        elif outformat == "pdf":
            fileout = filename+".pdf"
            f_out = os.path.join(outfolder, fileout)
            plt.savefig(f_out,format="pdf")
        else:
            log.warning(f"Unsupported output format {outformat}, defaulting to PDF")
            fileout = filename+".pdf"
            f_out = os.path.join(outfolder, fileout)
            plt.savefig(f_out,format="pdf")

        fig = plt.gcf()
        plt.close(fig)



# deprecated alias, kept for backward compatibility #
def plot_rxns(*args, **kwargs) -> None:
    """Deprecated: use plot_reactions()."""
    warnings.warn("plot_rxns() is deprecated, use plot_reactions() instead",
                  DeprecationWarning, stacklevel=2)
    return plot_reactions(*args, **kwargs)


## analyze topology ##
# get degrees #
def get_degrees(df:pd.core.frame.DataFrame=None, target_atoms:tuple[int|str,...]=None ) -> list[list[int]]:
    df_work = df.copy()
    degrees = [None]*len(df_work)
    for idx, (df_idx, frame) in enumerate(df_work.iterrows()):
        nxg = frame["graph"]

        if target_atoms is None:
            nodes = list(nxg.nodes())
        else:
            nodes = list(node for node, node_data in nxg.nodes(data=True)
                            if node_data.get('type') in target_atoms 
                            or node_data.get('element') in target_atoms)

        degrees[idx] = list(d for n, d in nxg.degree(nodes))
    return degrees

    
def find_minimum_cycle_basis(df: pd.DataFrame = None, min_size: int = 7, max_block_size: int = None) -> list[list[int]]:
    """
    Computes the Minimum Cycle Basis (MCB) for each frame in the trajectory.
    This method identifies the Smallest Set of Smallest Rings (SSSR) by 
    decomposing the graph into biconnected components. 

    Args:
        df (pd.DataFrame, optional): Input DataFrame containing 'graph' column. 
            Defaults to self.backbone or self.frames.
        min_size (int): Minimum number of nodes for a cycle to be included.
        max_block_size (int, optional): Safety threshold. Blocks with more nodes 
            than this will be skipped to avoid O(n^3) complexity stalls.

    Returns:
        list: A nested list [frames][cycles][node_ids].

    Complexity Note:
    The MCB algorithm is O(m^3 * n). For dense networks that 'gel' into a 
    single large block, max_block_size is highly recommended to avoid 
    computational stalls.
    """
    # Select data source
    df_in = df
    all_frame_cycles = [None] * len(df_in)

    for idx, (df_idx, frame) in enumerate(df_in.iterrows()):
        nxg = frame["graph"]
        frame_basis = []

        # Use biconnected components to isolate cyclic parts of the network
        for block_nodes in nx.biconnected_components(nxg):
            block_len = len(block_nodes)
        
            if block_len < min_size:
                continue
        
            # Check safety limit for computational cost
            if max_block_size and block_len > max_block_size:
                log.warning(
                    f"Skipping large block ({block_len} nodes) in frame {idx}. "
                    f"Increase max_block_size if this analysis is required."
                )
                continue
            
            subgraph = nxg.subgraph(block_nodes)
        
            # MCB calculation (the heavy lifting)
            block_basis = nx.minimum_cycle_basis(subgraph)
        
            # Filter and store results
            frame_basis.extend([c for c in block_basis if len(c) >= min_size])
    
        all_frame_cycles[idx] = frame_basis
    
    return all_frame_cycles

