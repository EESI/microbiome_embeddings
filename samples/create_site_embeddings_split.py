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

row = int(argv[1]) - 1
files = open('/mnt/HA/groups/rosenGrp/embed_ag_samples/ag_sample_file_split.txt').readlines()[row].rstrip().split(',')

samples = six.moves.cPickle.load(open('/mnt/HA/groups/rosenGrp/embed_ag_samples/sample_site.pkl','rb'))
site_embedding = {site:np.zeros((257,)) for site in set(samples.values())} # num reads + num nodes

for i,f in enumerate(files):
    
    print('Embedding sample %s/%s.' % (i,len(files)))
    sys.stdout.flush()

    sample = f.split('/')[-1].split('_')[0]
    
    if sample in samples:
        site = samples[sample]
        emb = np.loadtxt(files[0],delimiter=',',skiprows=1,usecols=range(1,257))
        emb_n = emb.shape[0]
        emb_sum = emb.sum(axis=0)
        site_embedding[site][0] += emb_n
        site_embedding[site][1:] += emb_sum

for site in site_embedding:
    s = site.split(':')[1]
    print('Saving %s.' % (s))
    np.savetxt('/mnt/HA/groups/rosenGrp/embed_ag_samples/out/site_embedding_split/seqs_10_256_5_50_10_1e-06_100_embedding_np_%s_%s.csv' % (s,row),
            site_embedding[site],delimiter=',')
