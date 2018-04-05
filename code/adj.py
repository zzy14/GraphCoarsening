# coding:utf-8
# 不更改节点的编号

import sys
from collections import defaultdict

if len(sys.argv) < 3:
    print "need input and output!"
    exit(1)

inputfile = sys.argv[1]
outputfile = sys.argv[2]

# read graph
graph = defaultdict(list)
count = 0
with open(inputfile, 'r') as edgeFile:
    while 1:
        l = edgeFile.readline()
        if l == '':
            break
        vec = l.split()
        node1, node2 = [int(x) for x in vec[:2]]
        if node2 == node1:
            continue
        graph[node1].append(node2)
        # graph[node2].append(node1)
    for n in graph.keys():
        clean = set(graph[n])
        graph[n] = list(clean)
        count += len(graph[n])

nsize = max(graph.keys())+1
fout = open(outputfile, 'w')
fout.write('{} {}\n'.format(nsize, count))
for i in range(nsize):
    if i in graph:
        fout.write('{} {}\n'.format(len(graph[i]), ' '.join([str(x) for x in graph[i]])))
    else:
        fout.write('0\n')
fout.close()
