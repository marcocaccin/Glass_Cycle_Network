# Glass_Cycle_Network

Small utility for visualising 2D networks from atomic environments of amorphous silica.

- Read in an Atoms configuration of a SiO2 network, already without the bridging oxygens between Si atoms (node=Si atom, edge=bridging O)
- Create a graph from it
- Identify all cycles in the network up to a given size
- Purge all cycles that are not minimal i.e., cycles that are composed of two or more smaller cycles.
- Identify neighbouring cycles
- Generate a dual network in which each node is a Si ring (or cycle), and each edge is a shared bridging oxygen.
- Print it with pretty colours


### Requirements:
numpy
ase
matscipy
networkx
