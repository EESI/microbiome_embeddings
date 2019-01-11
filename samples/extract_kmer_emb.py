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

from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence

import embed_functions as emb

name = 'ag'
a = 1e-05

path_model = argv[1]

fn_model_base = path_model.split('/')[-1]
k = int(fn_model_base.split('_')[1])
d = int(fn_model_base.split('_')[2])
fn_model_base = '_'.join(fn_model_base.split('_')[1:-1])
fn_out = '%s_%s_total_kmers_split.pkl' % (name,fn_model_base)

dir_totalkmers = '/mnt/HA/groups/rosenGrp/embed_samples/data/ag/total_kmers'
path_totalkmers = os.path.join(dir_totalkmers,fn_out)
total_kmers = six.moves.cPickle.load(open(path_totalkmers,'rb'))

fn_out = 'model_kmeremb_%s_%s_raw.csv.gz' % (fn_model_base,a)
dir_out = '/mnt/HA/groups/rosenGrp/embed_ag_samples/out'
path_out = os.path.join(dir_out,fn_out)

model = Word2Vec.load(path_model)
d = model.layer1_size
model = model.wv

with gzip.open(path_out,'wb') as f:
    header = 'kmer' + ',' + ','.join(str(i+1) for i in range(d)) + '\n'
    f.write(header.encode('utf-8'))
    for kmer in total_kmers:
        embedding = model[kmer] * a/(a + total_kmers[kmer])
        line = kmer + ',' + ','.join(str(e) for e in embedding) + '\n'
        f.write(line.encode('utf-8'))
