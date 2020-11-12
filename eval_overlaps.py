#!/usr/bin/python

"""summarizes size and frequency of overlaps"""
"""prepares values for the configuration"""


import numpy as np
import matplotlib.pyplot as plt


clusterstring = 'NCe1n0_3v3t_overlaps_all'
clusterfile = open(clusterstring + '.txt', 'r')
cluster_count = []
cluster_lengths = []
cluster_tuples = {}
number_of_clusters = 0
block_sum = 0
for line in clusterfile:
    #index;count;node_1,...,node_x
    ix, count, nodes = line.strip(' \n').split(';')
    cluster = list(map(int, nodes.split(',')))
    cluster_lengths.append(len(cluster))
    cluster_count.append(int(count))
    tuple = (int(count), len(cluster))
    if tuple in cluster_tuples:
        cluster_tuples[tuple] += 1
    else:
        cluster_tuples[tuple] = 1

    number_of_clusters += 1
    block_sum += len(cluster)*int(count)
clusterfile.close()

print('block sum', block_sum)


print('cluster lengths results')
print('min', min(cluster_lengths))
print('max', max(cluster_lengths))
print('mean', np.mean(cluster_lengths))
print('std', np.std(cluster_lengths))
print('var', np.var(cluster_lengths))
print('sum', sum(cluster_lengths))

len_dict = {}
for x in range(min(cluster_lengths), max(cluster_lengths)+1):
    len_dict[x] = 0
for x in cluster_lengths:
    len_dict[x] += 1
for k, v in sorted(len_dict.items(), key=lambda item: item[0]):
    per_cent = round(v/len(cluster_lengths), 2)
    if per_cent > 0:
        print(k, per_cent)


print('cluster count results')
print('min', min(cluster_count))
print('max', max(cluster_count))
print('mean', np.mean(cluster_count))
print('std', np.std(cluster_count))
print('var', np.var(cluster_count))
print('sum', sum(cluster_count))

cnt_dict = {}
for x in range(min(cluster_count), max(cluster_count)+1):
    cnt_dict[x] = 0
for x in cluster_count:
    cnt_dict[x] += 1
for k, v in sorted(cnt_dict.items(), key=lambda item: item[0]):
    per_cent = round(v/len(cluster_count), 2)
    if per_cent > 0:
        print(k, per_cent)


print('overlaps c*n')
configfile = open(clusterstring + '_config.txt', 'w')
print('count', number_of_clusters)
for k, v in sorted(cluster_tuples.items(), key=lambda item: item[1], reverse=False):
    per_cent = round(v/number_of_clusters, 2)
    if per_cent > 0:
        sz = str(k).strip('()').replace(',','')
        print('blocks initial', sz, v, per_cent)
        configfile.write('blocks initial ')
        configfile.write(sz + ' ') 
        configfile.write(str(v) + ' ')
        configfile.write(str(per_cent) + '\n')

configfile.close()

#bin_min = min(cluster_lengths)
#bin_max = 15 #max(cluster_lengths)
#plt.hist(cluster_lengths, bins=bin_max+1-bin_min, range=(bin_min,bin_max))
#plt.show()

#bin_min = min(cluster_count)
#bin_max = 15 #max(cluster_count)
#plt.hist(cluster_count, bins=bin_max+1-bin_min, range=(bin_min,bin_max))
#plt.show()

#log log plot
# bin_min = min(cluster_count)
# bin_max = max(cluster_count)
# bins = 10**(np.linspace(np.log10(bin_min), np.log10(bin_max), 100))
# counts, edges = np.histogram(cluster_count, bins, density=True)
# centers = (edges[1:] + edges[:-1]/2)
# plt.plot(centers, counts, '.')
# plt.xscale("log")
# plt.yscale("log")
# plt.show()
