import os
from collections import defaultdict
import numpy as np
DATA_PATH = './'
SUPPORT_NUM = 5

used_group = {}
idx = 0
with open(os.path.join(DATA_PATH, 'look_up_raw.txt'), 'r') as f:
    with open(os.path.join(DATA_PATH, 'look_up.txt'), 'w') as fout:
        node_size, _ = f.readline().split()
        node_size = int(node_size)
        labeled_node_size = 0
        count = 0.
        while 1:
            l = f.readline()
            if l == '':
                break
            vec = l.split()
            for x in vec[1:]:
                if x not in used_group:
                    used_group[x] = str(idx)
                    idx += 1
            fout.write('{} {}\n'.format(vec[0], ' '.join([used_group[x] for x in vec[1:]])))
        print 'used group size', len(used_group)

with open(os.path.join(DATA_PATH, 'look_up.txt'), "r") as f:
    old = f.read()

with open(os.path.join(DATA_PATH, 'look_up.txt'), "w") as f:
    f.write('{} {}\n'.format(node_size, len(used_group)))
    f.write(old)

with open(os.path.join(DATA_PATH, 'coarse.txt'), 'r') as f:
    node_size, _ = f.readline().split()
    node_size = int(node_size)
    d = defaultdict(int)
    while 1:
        l = f.readline()
        if l == '':
            break
        vec = l.split()[1:]
        for x in vec:
            d[x] += 1
    print 'Num of groups:', len(d.keys())