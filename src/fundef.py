import numpy as np
from ase.io import read as aseread
from matscipy.neighbours import neighbour_list
import networkx as nx
from sys import argv
import itertools
import pandas as pd


def atoms_to_nxgraph(atoms, cutoff):
    ni, nj = neighbour_list('ij', atoms, cutoff)
    adjacency_matrix = np.zeros((len(atoms), len(atoms))).astype(np.int)
    for i, j in zip (ni, nj):
        adjacency_matrix[i,j] = 1
    graph = nx.from_numpy_matrix(np.array(adjacency_matrix))
    return graph


def minimal_cycles(graph, cutoff=9):
    all_cycles = []
    for node in graph.nodes():
        # Avoid A-B-A cycles with len(p) > 3. Avoid large non-minimal cycles with cutoff=9
        cycles = [set(p) for p in nx.algorithms.all_simple_paths(graph, node, node, cutoff=cutoff) if len(p) > 3]
        for cycle in cycles:
            if cycle not in all_cycles:
                all_cycles.append(cycle)
    # purge non minimal cycles and non-cycles
    for c0 in [cycle for cycle in all_cycles if len(cycle) > 6]:
        sub = nx.Graph(graph.subgraph(list(c0)))
        if sub.number_of_edges() != sub.number_of_nodes():
            all_cycles.remove(c0)
    return all_cycles


def cycle_dual_graph(all_cycles):
    # create the network of connected cycles
    cycle_adjacency = np.zeros((len(all_cycles), len(all_cycles))).astype(np.int)
    for i, ci in enumerate(all_cycles):
        for j, cj in enumerate(all_cycles):
            if j > i:
                if len(ci.intersection(cj)) > 0:
                    cycle_adjacency[i,j] = 1
    cycle_adjacency += cycle_adjacency.T
    # create the dual network: a node is a minimal ring, an edge is a shared edge between two rings (e.g. an O atom in the real system)
    graph_dual = nx.from_numpy_matrix(cycle_adjacency)
    return graph_dual


def get_angle_distribution(at, cutoff=1.98):

    from matscipy.neighbours import neighbour_list
    from itertools import combinations

    ii, jj, DD = neighbour_list('ijD', at, cutoff) 

    angles = []
    for atom in at:
        if atom.number == 8:
            neighs = np.where(ii == atom.index)[0]
            # Si_1, Si_2 = jj[neighs]
            # D_1, D_2 = DD[neighs]
            # print(np.linalg.norm(D_2 - D_1))
            # print(Si_1, Si_2)
            for Si_1, Si_2 in combinations(jj[neighs], 2):
                angle = at.get_angle([Si_1, atom.index, Si_2])
                angles.append(angle)
    return np.array(angles)