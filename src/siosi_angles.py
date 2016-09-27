import numpy as np
from ase.io import read as aseread
from matscipy.neighbours import neighbour_list
import networkx as nx
from sys import argv

from fundef import minimal_cycles
from itertools import combinations
import _matscipy

########################################################################
########################################################################


at = aseread('../data/bilayer_TSmin.xyz', format='extxyz')

top = at[np.where(at.positions[:,2] > -.5)[0]]
top_si = top[np.where(top.get_atomic_numbers() == 14)[0]]
top_o = top[np.where(top.get_atomic_numbers() == 8)[0]]

cutoff = 1.97
ii, jj = neighbour_list('ij', top, cutoff)
neighbour_si = []
neighbour_o = []
angles = []
graph = nx.Graph()
for atom in top:
    if atom.number == 8:
        neighs = np.where(ii == atom.index)[0]
        for Si_1, Si_2 in combinations(jj[neighs], 2):
            graph.add_edge(Si_1, Si_2)
            graph[Si_1][Si_2]['oxygen'] = atom.index
            neighbour_o.append(atom.index)
            neighbour_si.append([Si_1, Si_2])
            angle = at.get_angle([Si_1, atom.index, Si_2])
            angles.append(angle)


dists = _matscipy.distances_on_graph(ii,jj)
adj_si =  np.zeros((len(top_si), len(top_si))).astype(np.int)

si_indices = np.where(top.get_atomic_numbers() == 14)[0]
for i_si, dists_row in enumerate(dists):
    if i_si in si_indices:
        adj_si[i_si] = (dists_row[si_indices] == 2)

all_cycles = minimal_cycles(graph, cutoff=9)
angles_ext = {}
all_angles = {}
for cycle in all_cycles:
    sub = graph.subgraph(cycle)
    n_ring = sub.number_of_nodes()
    for si_si_o in sub.edges(data='oxygen'):
        angle = at.get_angle([si_si_o[0], si_si_o[2], si_si_o[1]])
        if n_ring not in all_angles.keys():
            all_angles[len(sub)] = [angle]
        else:
            all_angles[len(sub)].append(angle)
        angle_name = '%d-%d' % (si_si_o[0], si_si_o[1])
        if angle_name not in angles_ext.keys():
            angles_ext[angle_name] = [[n_ring, angle]]
        else:
            angles_ext[angle_name].append([n_ring, angle])
angles = angles_ext.values()
angles_internal = [a for a in angles if len(angles) == 2]
angles = {}
for (n1, angle), (n2, _) in angles_internal:
    n1, n2 = np.sort([n1, n2])
    name = '%d-%d' % (n1, n2)
    angle = angle * 180 / np.pi
    if name not in angles.keys():
        angles[name] = [angle]
    else:
        angles[name].append(angle)

meanstd = {}
for k,v in angles.iteritems():
    meanstd[k] = [np.mean(v), np.std(v), len(v)]
data = []
for k,v in meanstd.iteritems():
    n1, n2 = k.split('-')
    if n1 == '4' or n2 =='9':
        pass
    else:
        data.append([int(n1), int(n2), v[0], v[1]])
data = np.array(data)

dat = np.zeros((4,4))
means = dat.copy()
stds = dat.copy()
for i,j,mm, ss in data:
    means[i-5, j-5] = mm
    stds[i-5, j-5] = ss
    means[j-5, i-5] = mm
    stds[j-5, i-5] = ss
for mm, ss in zip(means, stds):
    line = ['%.1f \pm %.1f' % (m,s) for m, s in zip(mm,ss)]
    line = ' & '.join(line)
    print(line)

import matplotlib.pyplot as plt
import seaborn as sns
plt.close()
from pylab import rcParams
rcParams['figure.figsize'] = 16, 12
sns.set_style('white')
sns.set_context('paper', font_scale=2, rc={"lines.linewidth": 2})

f, axarr = plt.subplots(2, 2)

for i, jj in zip([5,6,7,8], [[0,0], [0,1], [1,0], [1,1]]):
    angles = np.array(all_angles[i]) *180 / np.pi
    axarr[jj[0], jj[1]].hist(angles, bins=30, label='%d' % i, range=(110,170), color='grey', normed=False)
    axarr[jj[0], jj[1]].locator_params(axis='y',nbins=7)
    # axarr[jj[0], jj[1]].locator_params(axis='x',nbins=3)
    axarr[jj[0], jj[1]].set_xticks([110,130,150,170])
axarr[0,0].set_ylabel(r'Count')
axarr[1,0].set_ylabel(r'Count')
axarr[1,1].set_xlabel(r'Si-O-Si angle')
axarr[1,0].set_xlabel(r'Si-O-Si angle')
plt.show()
plt.legend(markerscale=0)
plt.tight_layout()
plt.savefig('../figs/siosi_angles.pdf')

