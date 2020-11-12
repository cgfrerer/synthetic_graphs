#!/usr/bin/python

#produces graph in terms of overlapping communities

import networkx as nx
import numpy as np
import random
from random import sample
import time
import sys
import getopt
#from this project
import inout


def create_graph(communityfilename, number_of_vertices, number_of_edges, degree_pareto, degree_min, startindex, minimum_community_size):

    #produce graph with required number of vertices
    G = nx.Graph()
    rg = range(startindex, startindex+number_of_vertices)
    G.add_nodes_from(rg)

    i = 0
    communityfile = open(communityfilename + '.txt', 'r')
    for c in communityfile:

        community = list(map(int, c.strip(' \n').split(' ')))
        if len(community) < minimum_community_size:
            continue

        # e.g. community = [1 107 128 144 152 170 174 177 180 194 195]

        #distribution of vertex degrees
        community_degrees = []
        while len(community_degrees) < len(community):
           d = int(round(degree_min * random.paretovariate(degree_pareto)))
           if d < len(community):
               community_degrees.append(d)
       
        #correction if sum of vertex degrees is odd
        if sum(community_degrees) % 2 == 1:
            dx = community_degrees.index(min(community_degrees))
            community_degrees[dx] += 1


        #count the existing vertex degrees in the community
        allocated_degrees = {}
        for j in community:
            allocated_degrees[j] = 0
        for j in range(0, len(community)-1):
            u = community[j]
            for k in range(j+1, len(community)):
                v = community[k]
                if G.has_edge(u, v):
                    allocated_degrees[u] += 1
                    allocated_degrees[v] += 1

        #distribution of edges: flip a coin
        community_degrees.sort()
        flip = 0
        while len(community_degrees) > 1:
            if flip == 0:
                #smallest remaining degree and remaining node with smallest allocated degree
                degree_choice = community_degrees[0]
                community_degrees = community_degrees[1:]
                node_zero = min(allocated_degrees, key=allocated_degrees.get)
                flip = 1

            elif flip == 1:
                #highest remaining degree and acording node
                degree_choice = community_degrees.pop()
                node_zero = max(allocated_degrees, key=allocated_degrees.get)
                flip = 0

            assign_degree = degree_choice - allocated_degrees[node_zero]
            allocated_degrees.pop(node_zero)
            community.remove(node_zero)

            if assign_degree > 0:

                #if only one other vertex is left, the situation is clear
                #also the case if the left-over community is smaller equal the degree
                if len(community) == 1 or len(community) <= assign_degree:
                     for ix in community:
                         G.add_edge(node_zero, ix)
                         allocated_degrees[ix] += 1

                else:
                     community_sub = []
                     for ix in community:
                         if G.has_edge(node_zero, ix) == False:
                             community_sub.append(ix)

                     to_be_allocated = min(assign_degree, len(community_sub))
                     choice = sample(community_sub, to_be_allocated)
                     for ix in choice:
                         G.add_edge(node_zero, ix)
                         allocated_degrees[ix] += 1

        if i % 100000 == 0:
            number_of_current_edges = G.number_of_edges() #expensive operation
            print('round, edges', i, number_of_current_edges)


        i += 1
    communityfile.close()

    #count unconnected vertices and integrate them
    unconnected = 0
    for v in G:
        if G.degree(v) == 0:
            unconnected += 1
    print('number of unconnected vertices', unconnected)
    for v in G:
        if G.degree(v) == 0:
            u = np.random.randint(startindex, startindex+number_of_vertices)
            if u != v:
                G.add_edge(u, v)

    #result of graph generation
    print('number of edges ', G.number_of_edges())
    return G



def complement_edges(G, number_of_edges, triangle_portion, startindex):
    #create connections throughout communities

    #create a dedicated number of triangles
    e = int((number_of_edges - G.number_of_edges()) * triangle_portion)
    #e = int(G.number_of_edges()*triangle_portion)
    while e > 0:
        u = np.random.randint(startindex, startindex+number_of_vertices)
        neighbors = list(G.neighbors(u))
        if len(neighbors) >= 2:
            v = np.random.choice(neighbors, 2, replace=False)
            if G.has_edge(v[0], v[1]) == False:
                G.add_edge(v[0], v[1])
                e -= 1
    print('number of edges after triangle creation', G.number_of_edges())

    #create random connections to complete edges
    e = number_of_edges - G.number_of_edges()
    while e > 0:
        u = np.random.randint(startindex, startindex+number_of_vertices)
        v = np.random.randint(startindex, startindex+number_of_vertices)
        if u != v and G.has_edge(u, v) == False:
            G.add_edge(u, v)
            e -= 1
    print('number of edges after explorartory integration', G.number_of_edges())
    return G





if __name__ == "__main__":
    start_time = time.time()

    configurationfile = ''

    ##default parameters
    communityfilename = 'communities-pokec-NCe1n0'
    number_of_vertices = 990908
    number_of_edges = 19684191
    startindex = 0
    minimum_community_size = 6
    triangle_portion = 0.5   #[0.0, 1.0]
    #derrived in another script from communityfile
    degree_pareto = 1.33735194561
    degree_min = 2

    if len(sys.argv) >= 2: #configuration is provided
        configurationfile = sys.argv[1]
       
        configuration, initial_blocks, extended_blocks = inout.read_configuration(configurationfile)
        if 'communities_filename' in configuration.keys():
            communityfilename = configuration['communities_filename']
        if 'number_of_vertices' in configuration.keys():
            number_of_vertices = configuration['number_of_vertices']
        if 'number_of_edges' in configuration.keys():
            number_of_edges = configuration['number_of_edges']
        if 'triangle_portion' in configuration.keys():
            triangle_portion = configuration['triangle_portion']

    print('communities from', communityfilename)
    print('number of vertices', number_of_vertices)
    print('number of edges', number_of_edges)
    print('triangle portion', triangle_portion)
    print('minimum community size', minimum_community_size)


    G = create_graph(communityfilename, number_of_vertices, number_of_edges, degree_pareto, degree_min, startindex, minimum_community_size)
    complement_edges(G, number_of_edges, triangle_portion, startindex)
    inout.write_graph_to_adjlist(G, 'syngraph-' + communityfilename)

    duration = time.time() - start_time
    print('program run took ' + str(duration) + ' seconds')
    
