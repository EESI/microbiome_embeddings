#!/usr/bin/env python

import sys
from sys import argv
import os
import zipfile
from os.path import splitext
import gzip
import csv
import six.moves.cPickle
import collections

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

import numpy as np
import pandas as pd
from itertools import product, islice
from operator import itemgetter
from sklearn.manifold import TSNE
from scipy import spatial
import random
import math

from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence

import embed_functions as emb
import embed_params as p

fn_row = int(argv[1]) - 1
query_fn = open('fns.dat', 'r').readlines()[fn_row].rstrip()
kegg_fn = query_fn.replace('query', 'kegg')

if not os.path.exists(query_fn):
    print('Aborting; query kmer file does not exist.')
    sys.exit(0)

if not os.path.exists(kegg_fn):
    print('Aborting; kegg kmer file does not exist.')
    sys.exit(0)

hits_fn = query_fn.replace('query', 'hits').replace('pkl', 'csv.gz')
sim_dir = os.path.expanduser('~/embedding/cosin')

if not os.path.exists(os.path.join(sim_dir, hits_fn)):
    print("Loading query and reference.")
    dat_kegg = six.moves.cPickle.load(open(kegg_fn, 'rb'))
    dat_query = six.moves.cPickle.load(open(query_fn, 'rb'))
    ids_kegg, embed_kegg = dat_kegg['ids'], dat_kegg['embeddings']
    ids_query, embed_query = dat_query['ids'], iter(dat_query['embeddings'])

    print("Calculating cosine similarity.")
    nn_min_idx = list()
    counter = 1

    while True:
        
        print('Processing batch ' + str(counter) + '.')
        embed_query_batch = list(islice(embed_query,10000))
        
        if len(embed_query_batch) == 0:
            break
        
        nn = spatial.distance.cdist(embed_query_batch, embed_kegg, metric='cosine')
        nn_min_idx.extend(np.argmin(nn, axis=1))
        counter += 1

    nn_org_id = [ids_kegg[i] for i in nn_min_idx]

    print("Exporting as dataframe.")
    z = zip(ids_query, nn_min_idx, nn_org_id)
    df = pd.DataFrame(list(z), columns=['query', 'reference', 'organism:16s'])
    df.to_csv(os.path.join(sim_dir, hits_fn), compression='gzip', index=False)

else:
    print(hits_fn + ' exists')
