#!/usr/bin/env python

from sys import argv
import os
import zipfile
from os.path import splitext, isfile
import gzip
import csv
import six.moves.cPickle
import collections
from shutil import copyfile

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

import numpy as np
import pandas as pd
from itertools import product
from operator import itemgetter
from sklearn.manifold import TSNE
import random
import math
from scipy import sparse

from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence

import embed_functions as emb
import embed_params as p

if argv[1] == 'query' or argv[1] == 'kegg':
  name = argv[1]
else:
  raise NameError('Must provide query or kegg as input.')

fn_row = int(argv[2])
model_fn = open('fns.dat','r').readlines()[fn_row].rstrip()
embed_fn = name + '_' + model_fn
k = int(model_fn.split('_')[1])

embed_dir = os.path.expanduser('~/embedding/embeddings')

kmers_fn_in = name + '_' + str(k) + '_kmers.csv.gz'
kmers_fn = name + '_' + str(k) + '_kmers_' + str(fn_row) + '.csv.gz'

ids_fn = name + '_' + str(k) + '_ids' + '.pkl'

if not os.path.exists(os.path.join(embed_dir,embed_fn)):

    while True:
        print('Copying ' + kmers_fn_in)
        copyfile(kmers_fn_in,kmers_fn)
        if os.path.getsize(kmers_fn_in) == os.path.getsize(kmers_fn):
            break

    ids = six.moves.cPickle.load(open(ids_fn,'rb'))['ids']

    model = Word2Vec.load(model_fn)
    
    embeddings = emb.read2vec(kmers_fn,model,ids,k=k)

    print('Saving ids and embeddings.')
    six.moves.cPickle.dump({'ids':ids,'embeddings':embeddings},
            open(os.path.join(embed_dir,embed_fn),'wb'),protocol=4)

else:

    print(embed_fn + ' exists')
