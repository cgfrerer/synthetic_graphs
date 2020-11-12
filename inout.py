#!/usr/bin/python

import networkx as nx

"""input and output procedures"""


def write_communities(filename, communities):
	communitiesfile = open(filename, 'w')
	for community in communities:
		if len(community) > 2:
			communitiesfile.write(' '.join(map(str, community)))
			communitiesfile.write('\n')
	communitiesfile.close()

def write_graph_to_adjlist(G, filename):
	print('write output')
	nx.write_adjlist(G, filename + '.adjlist')

def read_configuration(filename):
	configurationfile = open(filename, 'r')
	configuration = {}
	initial_blocks = {}
	extended_blocks = {}
	for l in configurationfile:
		l = l.strip('\n')
		if len(l) > 1:
			if l.startswith('#'):
				#ignore comments
				continue
			elif l.startswith('blocks'):
				#blocks type c n count more
				data = l.split(' ')
				c = int(data[2])
				n = int(data[3])
				count = int(data[4])
				if data[1].startswith('extend'):
					extended_blocks[(c, n)] = count
				else:
					initial_blocks[(c, n)] = count

			elif '=' in l:
				key, value = l.split('=')
				if 'number' in key:
					configuration[key] = int(value)
				elif 'size' in key:
					configuration[key] = int(value)
				elif 'portion' in key:
					configuration[key] = float(value)
				elif 'file' in key:
					configuration[key] = value

	configurationfile.close()
	return configuration, initial_blocks, extended_blocks

def read_metis_to_undirected_graph(graphfilepath, startindex=1):
    #indices in metis graph file start at 1 !!
    G = nx.Graph()
    i = startindex-1
    graphfile = open(graphfilepath, 'r')
    for line in graphfile:
        ls = line.strip(' \n').split(' ')
        if i >= startindex:
            ls = [int(x) for x in ls]
            for j in ls:
                G.add_edge(i, j)
        elif i == startindex-1:
            no_nodes = ls[0]
            no_edges = ls[1]
            for j in range(startindex, startindex+int(no_nodes)):
                G.add_node(j)
        i += 1
    print('read graph done')
    print('vertices', G.number_of_nodes())
    print('edges', G.number_of_edges())
    graphfile.close()
    return G

def read_adjlist_to_undirected_graph(graphfilepath, reduce=False):
	G = nx.Graph()
	i = 0
	graphfile = open(graphfilepath, 'r')
	for line in graphfile:
		ls = line.strip(' \n').split(' ')
		if len(ls) > 1:
			ls = [int(x) for x in ls]
			startnode = ls[0]
			for j in ls: #j[0] is connected to all others
				if startnode != j:
					G.add_edge(startnode, j)
		i += 1
		#reduce read data for local test
		if reduce == True and i == 50000:
			break
	graphfile.close()
	return G

def read_adjlist_to_directed_graph(graphfilepath):
	G = nx.DiGraph()
	graphfile = open(graphfilepath, 'r')
	for line in graphfile:
		ls = line.strip(' \n').split(' ')
		if len(ls) > 1:
			ls = [int(x) for x in ls]
			startnode = ls[0]
			for j in ls: #j[0] is connected to all others
				if startnode != j:
					G.add_edge(startnode, j)

	print('read graph done')
	print('vertices', G.number_of_nodes())
	print('edges', G.number_of_edges())
	graphfile.close()
	return G

def read_overlaps_per_node(metafilepath):
	overlapfactors = {}
	i=0
	metafile = open(metafilepath, 'r')
	#node; graphdegree; communities; sum_communitydegrees; overlap(sum/graphdegree)
	for line in metafile:
		if i == 0:
			i += 1
			continue
		ls = line.strip(' \n').split(';')
		n = int(ls[0])
		o = int(ls[4])
		if o < 0:
			o = 0
		overlapfactors[n] = o+1 #plus 1 to avoid zeros, of are used for probabilities
	metafile.close()
	return overlapfactors
