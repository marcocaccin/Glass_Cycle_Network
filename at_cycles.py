import numpy as np
from quippy.atoms import Atoms
import networkx as nx
from quippy.ringstat import distance_map

at = Atoms('../reduced_1ayer.xyz')
positions = at.positions[:,:2]
at.set_cutoff(3.8) # 2 * 1.6 + some extra for elongation. visual inspection first!
at.calc_connect()
adjacency_matrix = np.array(distance_map(at, at.n, at.n) == 1).astype(int)
graph = nx.from_numpy_matrix(np.array(adjacency_matrix))

all_cycles = []
for node in graph.nodes():
    # Avoid A-B-A cycles with len(p) > 3. Avoid large non-minimal cycles with cutoff=9
    cycles = [set(p) for p in nx.algorithms.all_simple_paths(graph, node, node, cutoff=9) if len(p) > 3]
    for cycle in cycles:
        if cycle not in all_cycles:
            all_cycles.append(cycle)

# purge non minimal cycles
for c0 in [cycle for cycle in all_cycles if len(cycle) > 6]:
    sub = nx.Graph(graph.subgraph(list(c0)))
    if sub.number_of_edges() > sub.number_of_nodes():
        all_cycles.remove(c0)


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
positions_dual = np.array([positions[np.array(list(cycle))].mean(axis=0) for cycle in all_cycles])

cycle_n_nodes = {}
for i, c in enumerate(all_cycles):
    cycle_n_nodes[i] = len(c)
nx.set_node_attributes(graph_dual, 'cycle_size', cycle_n_nodes)


# print out the graphs corresponding to the glass network and its dual.
from bokeh import palettes
import matplotlib.pyplot as plt
lengths = np.array(cycle_n_nodes.values())
lmin = lengths.min()
colours = [palettes.RdYlBu6[i - lmin] for i in lengths]
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.figsize'] = [24., 18.]
plt.rcParams['savefig.transparent'] = True
plt.clf()
nx.draw(graph_dual, positions_dual, node_color=colours, node_size=plt.rcParams['figure.figsize'][0]*lengths**2, alpha=0.7, width=2)
plt.savefig('onlydual.eps')
plt.clf()
nx.draw(graph_dual, positions_dual, node_color=colours, node_size=plt.rcParams['figure.figsize'][0]*lengths**2, alpha=0.7, width=2, style='dotted', linewidth=0)
nx.draw(graph, positions, node_color='#000000', node_size=plt.rcParams['figure.figsize'][0]*3**2, width=3)
plt.savefig('superimposed.eps')
plt.clf()
