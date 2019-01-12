import numpy as np
import time
import itertools
from sklearn.externals import joblib

k = 6
t1 = time.time()
file_name = '/mnt/HA/groups/rosenGrp/sw1/embed_samples/data/kegg/seqs.fasta'


err = 0
ids = []
read_idx = 0

kmer_table = np.zeros((16399, 4**k))
kmer2idx = {}
idx2kmer = []

for idx, kmer in enumerate(itertools.product('ATGC', repeat=k)):
    kmer2idx[''.join(kmer)] = idx
    idx2kmer.append(''.join(kmer))

print('kmer num:', len(kmer2idx), flush = True)

with open(file_name) as f:
    for line in f:
        l = line.strip('\n')
        if l[0] == '>':
            ids.append(l.strip('>'))
        else:
            M = len(l) - k + 1
            for n in range(M):
                nts = l[n:n+k]
                try:
                    kmer_table[read_idx, kmer2idx[nts]] += 1
                except:
                    err += 1
            read_idx += 1
            if read_idx % 163 == 0:
                t_diff = str(round((time.time() - t1)/60,1)) + ' min.'
                print('KEGG:\tProcessed ' + str(read_idx) + ' reads in ' + t_diff, flush = True)
                print('ERR:\t', err, flush = True)
print(np.sum(kmer_table, axis = 0), flush = True)

joblib.dump(kmer2idx, 'kegg_6_mer_kmer2idx.pkl')
joblib.dump(ids, 'kegg_6_mer_ids.pkl')
np.save('kegg_6_mer_kmerfreq.npy', kmer_table)
