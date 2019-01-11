#!/usr/bin/env python

import sys
from sys import argv
import os
import zipfile
from os.path import splitext, isfile
import gzip
import csv
import six.moves.cPickle
import collections
from shutil import copyfile

import numpy as np
import pandas as pd
from itertools import product
from operator import itemgetter
from sklearn.manifold import TSNE
import random
import math
import time
from glob import glob

from collections import defaultdict

import embed_functions as emb
import r2v_functions as r2v

name = 'kegg'
params = argv[1]

dir_split = '/mnt/HA/groups/rosenGrp/embed_cons/out/%s/total_kmers_split' % (name)
path_total_kmers = glob(dir_split + '/' + params + '/*total_kmers_split.pkl')

dir_out = '/mnt/HA/groups/rosenGrp/embed_cons/out/%s' % (name)
path_out = os.path.join(dir_out,'%s_%s_clust_total_kmers.pkl' % (name,params))


kmer_counts = defaultdict(int)
for i,f in enumerate(path_total_kmers):
    if i % 50 == 0:
        print('Merging file %s' % (i))
    d = six.moves.cPickle.load(open(f,'rb'))
    for kmer in d:
        kmer_counts[kmer] += d[kmer]

print('Normalizing counts to frequencies.')
total_kmers = np.sum([count for count in kmer_counts.values()])
kmer_counts = {kmer:count/total_kmers for kmer,count in kmer_counts.items()}
kmer_counts['total_kmers'] = total_kmers

print('Saving dictionary.')
six.moves.cPickle.dump(kmer_counts,open(path_out,'wb'),protocol=4)
