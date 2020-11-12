#!/usr/bin/python

"""find overlaps of communities"""

import time
import numpy as np
import matplotlib.pyplot as plt

start_time = time.time()

communitystring = 'NCe1n0'

threshold = 3
min_cluster_length = 3

min_community_length = 3
max_node_number = 1 #tbd
min_node_number = 1 #tbd

#read communities and put them into a list of lists
communities = []
entries = 0
communityfile = open(communitystring + '.txt', 'r')
for line in communityfile:
#	if len(line) > 1:
	community = list(map(int, line.strip(' \n').split(' ')))
	# e.g. community = [125, 1, 143, 165, 175, 8757, 6273, 8761]
	if len(community) >= min_community_length:
		communities.append(community)
		entries += len(community)
	for n in community:
		if n > max_node_number:
			max_node_number = n
		if n < min_node_number:
			min_node_number = n
communityfile.close()
print('communities', len(communities))
print('entries', entries)
print('vertex entries from to:', min_node_number, max_node_number)



#create a list of lists with the communities for each node
nodes_in_communities = [[] for i in range(max_node_number+1)]
for i in range(0, len(communities)):
	community = communities[i]
	for n in community:
		nodes_in_communities[n].append(i)


#find all communities where the node occurres
#find block of size at least a*b
overlapsfile = open(communitystring + '_' + str(min_cluster_length) + 'v' + \
               str(threshold) + 't_overlaps_all.txt', 'w')

o_index = 0
last_o_index = -1

while last_o_index < o_index*0.99:
	last_o_index = o_index

	seq = np.arange(0, max_node_number+1)
	np.random.shuffle(seq)

	for x in seq:

		if len(nodes_in_communities[x]) >= threshold:
			potential = []
			count = {}
			for c in nodes_in_communities[x]:
				if len(communities[c]) >= min_cluster_length:
					potential.append(c)
					for n in communities[c]:
						if n in count.keys():
							count[n] += 1
						else:
							count[n] = 1

			#reduce count by elements below threshold
			count = { k : v for k, v in count.items() if v >= threshold}

			if len(count) >= min_cluster_length:
				#find set of elements ('cluster') with min_cluster_length,
				#which occur more often than the threshold
				cluster = set()
				extend = []
				for k, v in sorted(count.items(), key=lambda item: item[1], reverse=True):
					if len(cluster) < min_cluster_length:
						cluster.add(k)
					else:
						extend.append(k)

				if len(cluster) == min_cluster_length:
					rounds = 1 + len(extend)
					last_cluster_count = 0
					while rounds > 0:
						cluster_count = 0
						for p in potential:
							b = cluster.intersection(set(communities[p]))
							if cluster == b:
								cluster_count += 1

						if cluster_count >= threshold:
							last_cluster_count = cluster_count
							rounds -= 1
							if len(extend) > 0:
								last_extension = extend.pop(0)
								cluster.add(last_extension)
						else:
							if len(cluster) > min_cluster_length:
								cluster.remove(last_extension)
							else:
								cluster.clear()
							rounds = 0

					if len(cluster) >= min_cluster_length:
						for p in potential:
							b = cluster.intersection(set(communities[p]))
							if cluster == b:
								for n in cluster:
									communities[p].remove(n)
									nodes_in_communities[n].remove(p)

						overlapsfile.write(str(o_index))
						overlapsfile.write(';')
						overlapsfile.write(str(last_cluster_count))
						overlapsfile.write(';')
						overlapsfile.write(','.join(map(str, cluster)))
						overlapsfile.write('\n')
						overlapsfile.flush()

						o_index += 1
						if o_index % 10000 == 0:
							print('overlaps', o_index)

overlapsfile.close()
print('')
print('total overlaps', o_index)


duration = time.time() - start_time
print('program run took ' + str(round(duration, 2)) + ' seconds')

