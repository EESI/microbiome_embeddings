#!/usr/bin/env python

import sys
from sys import argv
import os
import zipfile
from os.path import splitext, isfile
import gzip
import csv
import six.moves.cPickle
from collections import Counter
from shutil import copyfile

import numpy as np
import pandas as pd
from itertools import product, islice
from operator import itemgetter
from sklearn.manifold import TSNE
import random
import math
from scipy import sparse
import time
import multiprocessing as mp

from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence

import embed_functions as emb
import embed_params as p


def get_nn(line,r_ids,r_kmers,r_counts):
    
    max_count = 0
    tie_breaker = 0
    
    line = line.decode('utf-8').split()
    #q_kmer = set(line)
    q_count = Counter(line)
    q_kmer = set(q_count.keys())
    
    nn = []

    for r,r_kmer in enumerate(r_kmers):

            i_total = len(q_kmer & r_kmer)

            if i_total > max_count:

                tie_breaker = 0
                max_count = i_total
                nn = [r_ids[r]]

            elif i_total == max_count and i_total > 50:

                counts = q_count & r_counts[r]
                i_counts = sum(counts.values())

                if i_counts > tie_breaker:
                    tie_breaker = i_counts
                    nn = [r_ids[r]]
                elif i_counts == tie_breaker:
                    nn.append(r_ids[r])
            else:
                pass

    return nn

fn_row = int(argv[1]) - 1
nc = int(argv[2])
nl = int(argv[3])

k = ['6','10','15'][fn_row]

ref_fn_in = 'kegg_' + k + '_kmers.csv.gz'
query_fn_in = 'query_' + k + '_kmers.csv.gz'
hits_fn = 'ngrams' + k + '_kmers.csv.gz'

ngram_fn = 'ngrams' + k + '_kmers.pkl'
ngram_dir = os.path.expanduser('~/embedding/ngrams3')

ref_fn = 'tmp' + str(fn_row) + '_' + ref_fn_in
query_fn = 'tmp' + str(fn_row) + '_' + query_fn_in

r_ids_fn = 'kegg_' + k + '_ids' + '.pkl'
q_ids_fn = 'query_' + k + '_ids' + '.pkl'

if not os.path.exists(os.path.join(ngram_dir,ngram_fn)):

    while True:
        print('Copying ' + ref_fn_in)
        copyfile(ref_fn_in,ref_fn)
        if os.path.getsize(ref_fn_in) == os.path.getsize(ref_fn):
            break

    while True:
        print('Copying ' + query_fn_in)
        copyfile(query_fn_in,query_fn)
        if os.path.getsize(query_fn_in) == os.path.getsize(query_fn):
            break

    r_ids = six.moves.cPickle.load(open(r_ids_fn,'rb'))['ids']
    q_ids = six.moves.cPickle.load(open(q_ids_fn,'rb'))['ids']

    file_open = emb.open_file_method(ref_fn)
    r_file = file_open(ref_fn)
    
    print('Calculating reference counts.')
    r_lines = [line.decode('utf-8').split() for line in r_file]
    r_kmers = [set(line) for line in r_lines]
    r_counts = [Counter(line) for line in r_lines]

    file_open = emb.open_file_method(query_fn)
    q_file = file_open(query_fn)
    
    print('Finding nearest neighbors')

    def worker(lines):
 
        result = {id:get_nn(line,r_ids,r_kmers,r_counts) for id,line in lines}
        
        return result

    if __name__ == '__main__':

        sys.stdout = open('ngram_log_' + str(fn_row) + '.txt','w')

        print('Preparing nodes.')

        n_lines = nl
        n_cores = nc

        nns = {}
        q_ids_iter = iter(q_ids)

        pool = mp.Pool(processes=n_cores)
        
        b = 0
        while True:
            
            batch_lines = list(islice(q_file, n_lines * n_cores))
            batch_ids = list(islice(q_ids_iter, n_lines * n_cores))

            if not batch_lines or not batch_ids:
                break
            else:
                batch = list(zip(batch_ids,batch_lines))

            tmp = pool.map(worker,(batch[line:line + n_lines] for line in range(0,len(batch),n_lines)))

            for d in tmp:
                nns.update(d)


            if (b + 1) % 10 == 0:
                print('Running batch ' + str(b + 1) + '.')
                sys.stdout.flush()

            if (b + 1) % 100 == 0:
                print('Writing temp. output.')
                six.moves.cPickle.dump(nns,open('tmp_nns_' + str(fn_row) + '.pkl','wb'))
            
            b += 1

        pool.close()
        sys.stdout.close()

    print('Complete.')

    z = [(q,r[0]) for q,r in nns.items()]
    df = pd.DataFrame(list(z), columns=['query','organism:16s'])
    df.to_csv(os.path.join(ngram_dir, hits_fn), compression='gzip', index=False)

    print('Saving ngrams.')
    six.moves.cPickle.dump({'nns':nns},
            open(os.path.join(ngram_dir,ngram_fn),'wb'),protocol=4)

else:

    print(ngram_fn + ' exists')
