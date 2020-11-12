#!/usr/bin/python

"""random distributions for community lengths and nodes in communities"""
"""another distribution - blocks, different c*n size blocks"""

import numpy as np
import random
from random import sample
import time
import sys
import getopt
#from this project
import inout


def get_exponential_distributed_community_lengths(exp_scale_c, number_of_communities, min_community_length):
    community_lengths = np.random.exponential(scale=exp_scale_c, size=number_of_communities)+min_community_length
    for i in range(0, number_of_communities):
        community_lengths[i] = int(round(community_lengths[i]))
    community_lengths = community_lengths.astype(np.int)
    return community_lengths

def get_powerlaw_distributed_community_lengths(alpha, xmin, xmax, size):
    cumulative = []
    x = xmin
    #last_sum = pow(xmin, -1*alpha)
    last_sum = (alpha-1)/xmin * pow(x/xmin, -1*alpha)
    cumulative.append(last_sum)
    while len(cumulative) < (xmax-xmin) or cumulative[-1] < 1.:
        x += 1
        #last_sum += pow(x, -1*alpha)
        last_sum += (alpha-1)/xmin * pow(x/xmin, -1*alpha)
        cumulative.append(last_sum)
    cumulative.append(1.)

    randoms = np.random.rand(size)
    community_lengths = []
    for r in randoms:
        t = 0
        while cumulative[t] < r:
            t += 1
        community_lengths.append(xmin + t)
    return community_lengths

def get_exponential_distributed_nodes_occurrences(exp_scale_n, number_of_vertices):
    nodes_occurrences = np.random.exponential(scale=exp_scale_n, size=number_of_vertices)
    for i in range(0, number_of_vertices):
        nodes_occurrences[i] = int(round(nodes_occurrences[i]))
    nodes_occurrences = nodes_occurrences.astype(np.int)
    return nodes_occurrences

def get_lognormal_distributed_nodes_occurrences(mu, sigma, size):
    rng = np.random.default_rng()
    s = rng.lognormal(mu, sigma, size)
    s=[int(round(x)) for x in s]
    return s

def print_communities_creation_status(communities, goal):
    c_sum = 0
    for i in range(0, len(communities)):
        c_sum += len(communities[i])
    print('sum new communities', c_sum)
    if goal != 0:
        print('portion achieved', round(c_sum * 100 / goal))
        print('')


def initiate_blocks(block_tuple, number_of_blocks, new_communities, nodes_occurrences, community_lengths, blocks_c, blocks_n):

    block_size_c = block_tuple[0]
    block_size_n = block_tuple[1]
    sequence_c = []
    i = 0
    while i < len(community_lengths):
        j = int(community_lengths[i] / block_size_n)
        sequence_c.extend([i] * j)
        i += 1
    np.random.shuffle(sequence_c)

    sequence_n = []
    i = 0
    while i < len(nodes_occurrences):
        j = int(nodes_occurrences[i] / block_size_c)
        sequence_n.extend([i] * j)
        i += 1
    np.random.shuffle(sequence_n)
    #print('sequences initiated')

    pos_c = 0
    pos_n = 0
    i = 0
    while i < number_of_blocks:
        selected_communities = sequence_c[pos_c:pos_c+block_size_c]
        pos_c += block_size_c
        selected_nodes = sequence_n[pos_n:pos_n+block_size_n]
        pos_n += block_size_n
        for c in selected_communities:
            community_lengths[c] -= block_size_n
            new_communities[c].update(selected_nodes)

        for n in selected_nodes:
            nodes_occurrences[n] -= block_size_c

        blocks_c.append(selected_communities)
        blocks_n.append(selected_nodes)
        i += 1


def extend_blocks(block_tuple, number_of_blocks, new_communities, nodes_occurrences, community_lengths, blocks_c, blocks_n):

    block_size_c = block_tuple[0]
    block_size_n = block_tuple[1]
    #set sequence to proper selection
    sequence_n = []
    i = 0
    while i < len(nodes_occurrences):
        j = int(nodes_occurrences[i] / block_size_c)
        sequence_n.extend([i] * j)
        i += 1
    np.random.shuffle(sequence_n)

    pos_n = 0
    i = 0
    while i < number_of_blocks:
        b = np.random.randint(0, len(blocks_c))

        #selected_nodes = set(sample(sequence_n, block_size_n))
        selected_nodes = set(sequence_n[pos_n:pos_n+block_size_n])
        pos_n += block_size_n

        #exclude nodes that already appear in the block/overlap
        selected_nodes.difference_update(blocks_n[b])

        block_subselection = blocks_c[b].copy()
        for c in block_subselection:
            if community_lengths[c] < len(selected_nodes):
                block_subselection.remove(c)

        if len(block_subselection) >= block_size_c:
            selected_communities = sample(block_subselection, block_size_c)

            for c in selected_communities:
                community_lengths[c] -= len(selected_nodes)
                new_communities[c].update(selected_nodes)
            for n in selected_nodes:
                nodes_occurrences[n] -= block_size_c

            blocks_c.append(selected_communities)
            blocks_n.append(selected_nodes)
            i += 1

def distribute_remaining_nodes(new_communities, nodes_occurrences, community_lengths):
    sequence_of_remaining_nodes = []
    for i in range(0, len(nodes_occurrences)):
        n = nodes_occurrences[i]
        if n > 0:
            sequence_of_remaining_nodes.extend([i] * n)
        if n < 0:
            print(i, n)
    #print('len remaining sequence', len(sequence_of_remaining_nodes))
    np.random.shuffle(sequence_of_remaining_nodes)

    i = 0
    for c in range(0, len(community_lengths)):
        r = community_lengths[c]
        #only extend existing communities
        if len(new_communities[c]) > 0 and r > 0:
            selected_nodes = sequence_of_remaining_nodes[i:i+r]
            new_communities[c].update(selected_nodes)
            i += r

if __name__ == '__main__':
    start_time = time.time()

    number_of_vertices = 990908
    number_of_communities = 1580364
    outputfilename = 'test_output'

    #exp_scale_c = 5.2432
    exp_scale_c = 11.22185
    #min_community_length = 3
    min_community_length = 5
    #exp_scale_n = 24.75
    exp_scale_n = 22.68193

    initial_blocks = {}
    extended_blocks = {}

    configurationfile = ''

    if len(sys.argv) >= 2: #configuration is provided
        configurationfile = sys.argv[1]
        
        configuration, initial_blocks, extended_blocks = inout.read_configuration(configurationfile)
        if 'number_of_vertices' in configuration.keys():
            number_of_vertices = configuration['number_of_vertices']
        if 'number_of_communities' in configuration.keys():
            number_of_communities = configuration['number_of_communities']
        if 'communities_filename' in configuration.keys():
            outputfilename = configuration['communities_filename']

    print('vertices', number_of_vertices)
    print('communities', number_of_communities)
    print('outputfile', outputfilename)
    print('')

    community_lengths = get_exponential_distributed_community_lengths(exp_scale_c, number_of_communities, min_community_length)
    #community_lengths = get_powerlaw_distributed_community_lengths(alpha=1.55, xmin=min_community_length, xmax=1800, size=number_of_communities)
    print('sum community lengths', sum(community_lengths))

    nodes_occurrences = get_exponential_distributed_nodes_occurrences(exp_scale_n, number_of_vertices)
    #nodes_occurrences = get_lognormal_distributed_nodes_occurrences(size=number_of_vertices, mu=4.11041141929245, sigma=1.47061422357945)
    print('sum nodes occurrences', sum(nodes_occurrences))

    #ensure sum(nodes_occurrences) == sum(community_lengths)
    if sum(nodes_occurrences) < sum(community_lengths):
        diff = sum(community_lengths) - sum(nodes_occurrences)
        selected_nodes = np.random.randint(0, number_of_vertices, size=diff)
        for s in selected_nodes:
            nodes_occurrences[s] += 1
    elif sum(community_lengths) < sum(nodes_occurrences):
        diff = sum(nodes_occurrences) - sum(community_lengths)
        selected_communities = np.random.randint(0, number_of_communities, size=diff)
        for s in selected_communities:
            community_lengths[s] += 1

    print('')
    goal = sum(nodes_occurrences)

    new_communities = [None] * number_of_communities
    for i in range(0, number_of_communities):
        new_communities[i] = set()

    blocks_c = []
    blocks_n = []

    #produce blocks randomly
    for k, v in initial_blocks.items():
        print('blocks initial (c, n) count:', k, v)
        initiate_blocks(k, v, new_communities, nodes_occurrences, community_lengths, blocks_c, blocks_n)
        print_communities_creation_status(new_communities, goal)

    #extend the blocks by smaller overlaps
    for k, v in extended_blocks.items():
        print('blocks extended (c, n) count:', k, v)
        extend_blocks(k, v, new_communities, nodes_occurrences, community_lengths, blocks_c, blocks_n)
        print_communities_creation_status(new_communities, goal)

    distribute_remaining_nodes(new_communities, nodes_occurrences, community_lengths)
    print_communities_creation_status(new_communities, goal)

    inout.write_communities(outputfilename + '.txt', new_communities)

    duration = time.time() - start_time
    print('program run took ' + str(duration) + ' seconds')


