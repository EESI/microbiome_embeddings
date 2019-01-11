#!/usr/bin/env python

import sys
from sys import argv
import os
import zipfile
from os.path import splitext, isfile
import gzip
import csv
import six.moves.cPickle

import numpy as np
import pandas as pd
import random
import math
import time
from glob import glob

#path = '/mnt/HA/groups/rosenGrp/embed_ag_samples/out/site_embedding_split'
path = '/mnt/HA/groups/rosenGrp/embed_ag_samples/out/site_embedding_split_trainset'
files = glob(path + '/*csv')

samples = six.moves.cPickle.load(open('/mnt/HA/groups/rosenGrp/embed_ag_samples/sample_site.pkl','rb'))

out = {}
n = {}

for i,f in enumerate(files):
    
    site = f.split('/')[-1].split('_')[-2]
    emb = np.loadtxt(f,delimiter=',')

    print('Merged file %s/%s' % (i,len(files)))

    try:
        out[site] += emb[1:]
        n[site] += emb[0]
    except KeyError:
        out[site] = emb[1:]
        n[site] = emb[0]

embs = {}
for site in out:
    print('Saving %s.' % (site))
    embs[site] = out[site]/n[site]

embs = pd.DataFrame(embs)
#embs.to_csv('/mnt/HA/groups/rosenGrp/embed_ag_samples/out/seqs_10_256_5_50_10_1e-06_100_embedding_sites.csv')
embs.to_csv('/mnt/HA/groups/rosenGrp/embed_ag_samples/out/seqs_10_256_5_50_10_1e-06_100_embedding_sites_trainset.csv')
