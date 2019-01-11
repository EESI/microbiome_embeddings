import numpy as np
import time
import itertools
from sklearn.externals import joblib
import glob

k = 10
t1 = time.time()
file_list = glob.glob('/mnt/HA/groups/rosenGrp/embed_samples/data/ag/fasta/*.fasta')


err = 0
ids = []

kmer_table = np.zeros((len(file_list), 4**k))
kmer2idx = {}

for idx, kmer in enumerate(itertools.product('ATGC', repeat=k)):
    kmer2idx[''.join(kmer)] = idx
print('kmer num:', len(kmer2idx), flush = True)


for file_idx, file_name in enumerate(file_list):
    ids.append(file_name.strip('.fasta'))
    with open(file_name) as f:
        for line in f:
            l = line.strip('\n')
            if l[0] != '>':
                M = len(l) - k + 1
                for n in range(M):
                    nts = l[n:n+k]
                    try:
                        kmer_table[file_idx, kmer2idx[nts]] += 1
                    except:
                        err += 1
    if file_idx % 150 == 0:
        t_diff = str(round((time.time() - t1)/60,1)) + ' min.'
        print(' AG:\tProcessed ' + str(file_idx) + ' files in ' + t_diff, flush = True)
        print('ERR:\t', err, flush = True)
            
joblib.dump(kmer2idx, 'ag_10_mer_kmer2idx.pkl')   
joblib.dump(ids, 'ag_10_mer_ids.pkl')
np.save('ag_10_mer_kmerfreq.npy', kmer_table)            
print(np.sum(kmer_table[:,2084]), flush = True)
print('err kmer:', err, flush = True)
