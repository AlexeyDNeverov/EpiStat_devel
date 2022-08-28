#!/usr/bin/env python
import louvain
import igraph as ig
G = ig.Graph.Erdos_Renyi(100, p=5./100);
partition = louvain.find_partition(G, louvain.ModularityVertexPartition)
#optimiser = louvain.Optimiser()
##optimiser.consider_comms=louvain.ALL_COMMS
#
#G = ig.Graph.Erdos_Renyi(100, p=5./100);
#partition = louvain.ModularityVertexPartition(G);
#improv = 1;
#while improv > 0:
#	improv = optimiser.optimise_partition(partition);