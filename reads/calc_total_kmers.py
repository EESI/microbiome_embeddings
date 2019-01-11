#!/usr/bin/env python

import sys
sys.path.insert(0, '/data/sw1/embeddings/code/')

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
import r2v_functions as r2v

name_reads = 'kegg'
path_reads = '/home/sw424/embed_samples/data/%s/seqs.fasta' % (name_reads)
path_params = '/home/sw424/embed_samples/code/read_params.csv'
dir_out = '/home/sw424/embed_samples/data/%s/' % (name_reads)

model_idx = int(argv[1]) - 1
path_model = glob('/home/sw424/embed_samples/data/models/*model.pkl')
path_model.sort()
path_model = path_model[model_idx]

fn_model_base = path_model.split('/')[-1]
fn_model_base = '_'.join(fn_model_base.split('_')[1:-1])
fn_out = '%s_%s_total_kmers.pkl' % (name_reads,fn_model_base)
path_out = os.path.join(dir_out,fn_out)
model_fn = path_model.split('/')[-1]
k = int(model_fn.split('_')[1])

if not os.path.exists(dir_out):
    os.makedirs(dir_out)

print('Calculating kmer totals for %s using model %s.' % (path_reads,path_model))

total_kmers = r2v.calc_total_kmers(path_reads,path_model,k,verbose=True,v=10000)

six.moves.cPickle.dump(total_kmers,open(path_out,'wb'),protocol=4)
