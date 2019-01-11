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
import r2v_functions as r2v

param_row = int(argv[1])

name_sample = 'kegg'
path_params = '/home/sw424/embed_samples/code/read_params.csv'
path_totalkmers = '/home/sw424/embed_samples/out/total_kmers.pkl'
path_out = '/home/sw424/embed_samples/out/%s_reads_embeddings/' % (name_sample)
path_sample = '/home/sw424/embed_samples/data/%s/seqs.fasta' % (name_sample)

params = open(path_params,'r').readlines()[param_row].rstrip().split(',')
a = float(params[0])
path_model = params[1]

model_fn = path_model.split('/')[-1]
k = int(model_fn.split('_')[1])

dir_totalkmers = '/home/sw424/embed_samples/data/%s/' % (name_sample)
fn_model_base = path_model.split('/')[-1]
fn_model_base = '_'.join(fn_model_base.split('_')[1:-1])
fn_out = '%s_%s_total_kmers.pkl' % (name_sample,fn_model_base)
path_totalkmers = os.path.join(dir_totalkmers,fn_out)

#path_out = '/home/sw424/embed_samples/out/group_reads_embeddings/'
if not os.path.exists(path_out):
    os.makedirs(path_out)

#path_model = '/home/sw424/embed_samples/data/models/gg_%s_%s_5_%s_%s_%s_100_model.pkl' \
#        % (k,d,50,10,1e-06)
#path_samples = '/home/sw424/embed_samples/data/groups'
#sample_ids = ['group_1']
#path_samples = glob('/home/sw424/embed_samples/data/kegg/*.fasta')
#for path_sample in path_samples:
    #fn = '%s.fasta' % (sample)
    #path_sample = os.path.join(path_samples,fn)
r2v.embed_reads(path_sample,path_totalkmers,path_model,path_out,k=k,a=a,delim=None,verbose=True,v=1000)
