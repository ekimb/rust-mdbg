import networkx

g = networkx.read_graphml('reggraph-k10-p0.001-l12.u.graphml')

def save_cc(i,sg):
    networkx.write_graphml(sg, "reggraph-k10-p0.001-l12.u_wcc/wcc_{}.graphml".format(i)) 

[save_cc(list(wcc_nodes)[0], networkx.subgraph(g,wcc_nodes)) for wcc_nodes in networkx.weakly_connected_components(g)] 
