import numpy as np
from ase.io import read as aseread
import networkx as nx
import itertools
import pandas as pd
from bokeh import palettes
import matplotlib.pyplot as plt

from fundef import atoms_to_nxgraph, minimal_cycles, cycle_dual_graph

########################################################################
########################################################################

at = aseread('../data/reduced_1ayer.xyz', format='extxyz')
do_plots = True

top = at[np.where(at.positions[:,2] > -.5)[0]]
at = top[np.where(top.get_atomic_numbers() == 14)[0]]
cutoff = 3.8 # 2 * 1.6 + some extra for elongation. visual inspection first!

graph = atoms_to_nxgraph(at, cutoff)

all_cycles = minimal_cycles(graph, cutoff=9)

graph_dual = cycle_dual_graph(all_cycles)

cycle_n_nodes = {}
for i, c in enumerate(all_cycles):
    cycle_n_nodes[i] = len(c)
nx.set_node_attributes(graph_dual, 'cycle_size', cycle_n_nodes)


if do_plots:
    # print out the graphs corresponding to the glass network and its dual.
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['figure.figsize'] = [24., 18.]
    plt.rcParams['savefig.transparent'] = True
    
    positions = at.get_positions()[:,:2]
    positions_dual = np.array([positions[np.array(list(cycle))].mean(axis=0) for cycle in all_cycles])
    
    lengths = np.array(cycle_n_nodes.values())
    lmin = lengths.min()
    colours = [palettes.RdYlBu6[i - lmin] for i in lengths]
    plt.clf()
    nx.draw(graph_dual, positions_dual, node_color=colours, node_size=plt.rcParams['figure.figsize'][0]*lengths**2, alpha=0.7, width=2)
    plt.savefig('../figs/onlydual.eps')

    plt.clf()
    nx.draw(graph_dual, positions_dual, node_color=colours, node_size=plt.rcParams['figure.figsize'][0]*lengths**2, alpha=0.7, width=2, style='dotted', linewidth=0)
    nx.draw(graph, positions, node_color='#000000', node_size=plt.rcParams['figure.figsize'][0]*3**2, width=3)
    plt.savefig('../figs/superimposed.eps')
    plt.clf()


lengths = np.array(cycle_n_nodes.values())
smallest, largest = lengths.min(), lengths.max()
allsizes = np.arange(smallest, largest + 1)
neighbours = [[lengths[u] for u in graph_dual.neighbors(i)] for i in range(len(lengths))]

nneighs = [len(n) for n in neighbours]
inner_indices = []
for i, (length, nn) in enumerate(zip(lengths, nneighs)):
    if length == nn:
        inner_indices.append(i)

# frequency of ring of given size
freqs = {}
for i in allsizes:
    freqs[i] = (lengths == i).sum()
n_rings = np.float(np.sum(freqs.values()))
for i in allsizes:
    freqs[i] /= n_rings

# indices of neighbours of each ring
n_indices = {}
for i in allsizes:
    n_indices[i] = np.where(lengths == i)[0]
    
import string
table = []
for size in allsizes:
    all_neigh_size = [neighbours[i] for i in n_indices[size]]
    all_neigh_size = np.array(list(itertools.chain.from_iterable(all_neigh_size)))
    result = np.array([(all_neigh_size == i).sum() for i in allsizes.astype('float')]) 
    result = result / float(result.sum())
    table.append(result)
    # print out probability of neighbour size
    print(string.join(['%.2f' % res for res in result], sep=' & '))

import matplotlib.pyplot as plt
import seaborn as sns
from pylab import rcParams
rcParams['figure.figsize'] = 8, 6
plt.clf()
sns.set_context('paper', font_scale=2, rc={"lines.linewidth": 2})
[plt.plot(allsizes, v, 'o-', label=k) for k, v in zip(allsizes, table)]
sns.despine()
plt.xlabel("Neighbouring ring size")
plt.ylabel("Frequency")
plt.legend()
plt.savefig("../figs/neighbour_sizes.pdf")

plt.clf()
df = pd.DataFrame(data=np.array(table), index=[4,5,6,7,8,9], columns=[4,5,6,7,8,9])
sns.set_context('paper', font_scale=2, rc={"lines.linewidth": 2, 'figure.figsize': [8, 6]})
cbar_kws = { 'label': 'Frequency'}
ax = sns.heatmap(df, cmap="YlGnBu", cbar_kws=cbar_kws)
plt.ylabel(r"Ring size $N$")
plt.xlabel(r"Neighbouring ring size $M$")

rcParams['figure.figsize'] = 8, 6
plt.clf()
plt.hist(lengths, bins=[4,5,6,7,8,9,10], align='left')
plt.xlabel("Ring size")
plt.ylabel("Count")
plt.savefig("../figs/size_count.pdf")

