#!/usr/bin/env python
#http://louvain-igraph.readthedocs.io/en/latest/index.html
import louvain
import igraph as ig
from tabulate import tabulate
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-n", "--neg", dest="neg",
                  help="specifies negative graph layer in GraphML format")
parser.add_option("-p", "--pos", dest="pos",
                  help="specifies positive graph layer in GraphML format")
parser.add_option("-f", "--funct", dest="funct", default="modularity",
                  help="specifies optimizing functional")
parser.add_option("-r", "--repeats", dest="n_iters", type="int", default=0,
                  help="specifies number of repeats of the optimal partition searching")
parser.set_defaults(verbose=True)
(options, args) = parser.parse_args()
if not(options.pos):
    parser.error("Graph layer with positive weights on branches is required!")
if not((options.funct == "modularity") or 
	(options.funct == "surprise")):
	parser.error("Unsupported functional accounted (-f)!")
print("Optimizing functional= ",options.funct)
n_iters=options.n_iters

G_pos=ig.Graph.Read_GraphML(options.pos)
G_neg=ig.Graph()
if options.neg:
	G_negative=ig.Graph.Read_GraphML(options.neg)
	G_neg = G_negative.subgraph_edges(G_negative.es.select(weight_lt = 0), delete_vertices=False)
	G_neg.es['weight'] = [-w for w in G_neg.es['weight']]
optimiser = louvain.Optimiser()
optimiser.consider_comms=louvain.ALL_COMMS
if options.funct == "modularity":
	part_pos = louvain.ModularityVertexPartition(G_pos, weights='weight')
elif options.funct == "surprise":
	part_pos = louvain.SurpriseVertexPartition(G_pos, weights='weight')
diff=1
if options.neg:
	print("Graph with negative layer:")
	if options.funct == "modularity":
		part_neg = louvain.ModularityVertexPartition(G_neg, weights='weight')
	elif options.funct == "surprise":
		part_neg = louvain.SurpriseVertexPartition(G_neg, weights='weight')
	diff = optimiser.optimise_partition_multiplex([part_pos, part_neg],layer_weights=[1,-1])
else:
	print("Graph without negative layer:")
	diff = optimiser.optimise_partition(part_pos)
quality=part_pos.quality()
if options.neg:
	quality=quality-part_neg.quality()
membership=part_pos.membership
significance=louvain.SignificanceVertexPartition.FromPartition(part_pos).quality()
print("Q_0=",quality,"Significance_0=",significance)
while n_iters>0:
	if options.neg:
		diff = optimiser.optimise_partition_multiplex([part_pos, part_neg],layer_weights=[1,-1])
	else:
		diff = optimiser.optimise_partition(part_pos)
	n_iters=n_iters-1
quality=part_pos.quality()
if options.neg:
	quality=quality-part_neg.quality()
membership=part_pos.membership
print("Q=",quality,"Significance=",louvain.SignificanceVertexPartition.FromPartition(part_pos).quality())
print(tabulate(zip(G_pos.vs["id"],[x+1 for x in membership]),headers=['site','group'],tablefmt='plain'))