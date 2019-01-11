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

name_reads = 'ag'

path_reads = argv[1]
path_model = argv[2]

name_sample = path_reads.split('/')[-1].replace('.fasta','')
fn_model_base = path_model.split('/')[-1]
fn_model_base = '_'.join(fn_model_base.split('_')[1:-1])
fn_out = '%s_total_kmers_split.pkl' % (name_sample)
model_fn = path_model.split('/')[-1]
k = int(model_fn.split('_')[1])

dir_out = '/mnt/HA/groups/rosenGrp/embed_samples/data/%s/total_kmers_split/%s' % (name_reads,fn_model_base)
path_out = os.path.join(dir_out,fn_out)
if not os.path.exists(dir_out):
    os.makedirs(dir_out)

print('Calculating kmer totals for %s using model %s.' % (path_reads,path_model))

total_kmers = r2v.calc_total_kmers_split(path_reads,path_model,k,verbose=True,v=10000)
six.moves.cPickle.dump(total_kmers,open(path_out,'wb'),protocol=4)
