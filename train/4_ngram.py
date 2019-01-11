#!/usr/bin/env python

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
from itertools import product
from operator import itemgetter
from sklearn.manifold import TSNE
import random
import math
from scipy import sparse
import time

from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence

import embed_functions as emb
import embed_params as p

fn_row = int(argv[1]) - 1
kegg_fn_in = open('fns.dat','r').readlines()[fn_row].rstrip()
query_fn_in = kegg_fn_in.replace('kegg','query')
k = int(kegg_fn_in.split('_')[1])

ngram_fn = kegg_fn_in.replace('kegg','ngrams').replace('.csv.gz','.pkl')
ngram_dir = os.path.expanduser('~/embedding/ngrams')

kegg_fn = 'tmp' + str(fn_row) + '_' + kegg_fn_in
query_fn = 'tmp' + str(fn_row) + '_' + query_fn_in

kegg_ids_fn = 'kegg_' + str(k) + '_ids' + '.pkl'
query_ids_fn = 'query_' + str(k) + '_ids' + '.pkl'

if not os.path.exists(os.path.join(ngram_dir,ngram_fn)):

    while True:
        print('Copying ' + kegg_fn_in)
        copyfile(kegg_fn_in,kegg_fn)
        if os.path.getsize(kegg_fn_in) == os.path.getsize(kegg_fn):
            break

    while True:
        print('Copying ' + query_fn_in)
        copyfile(query_fn_in,query_fn)
        if os.path.getsize(query_fn_in) == os.path.getsize(query_fn):
            break

    kegg_ids = six.moves.cPickle.load(open(kegg_ids_fn,'rb'))['ids']
    query_ids = six.moves.cPickle.load(open(query_ids_fn,'rb'))['ids']

    file_open = emb.open_file_method(kegg_fn)
    kegg_file = file_open(kegg_fn)
    
    print('Calculating reference counts.')
    r_lines = [line.decode('utf-8').split() for line in kegg_file]
    r_kmers = [set(line) for line in r_lines]
    r_counts = [Counter(line) for line in r_lines]

    file_open = emb.open_file_method(query_fn)
    query_file = file_open(query_fn)

    nn_count = dict()
    nn = dict()
    
    print('Finding nearest neighbors')
    t1 = time.time()
    for q,q_line in enumerate(query_file):

        max_count = 0
        tie_breaker = 0

        q_line = q_line.decode('utf-8').split()
        q_kmer = set(q_line)
        q_count = Counter(q_line)
        
        for r,r_kmer in enumerate(r_kmers):

            i_total = len(q_kmer & r_kmer)

            if i_total > max_count:

                tie_breaker = 0
                max_count = i_total
                nn[query_ids[q]] = [kegg_ids[r]]

            elif i_total == max_count and i_total > 50:

                counts = q_count & r_counts[r]
                i_counts = sum(counts.values())

                if i_counts > tie_breaker:
                    tie_breaker = i_counts
                    nn[query_ids[q]] = [kegg_ids[r]]
                elif i_counts == tie_breaker:
                    nn[query_ids[q]].append(kegg_ids[r])
            else:
                pass

            nn_count[query_ids[q]] = max_count

        print('Read ' + str(q) + ' (len=' + str(len(q_kmer)) + '): ' + str(len(nn[query_ids[q]])) + \
                ' matches with score ' + str(max_count) + ', tie-breaker ' + str(tie_breaker))

    t2 = time.time()
    print('Time elapsed: ' + str(t2 - t1))

    print('Saving ngrams.')
    six.moves.cPickle.dump({'nn':nn,'nn_count':nn_count},
            open(os.path.join(ngram_dir,ngram_fn),'wb'),protocol=4)

else:

    print(ngram_fn + ' exists')
