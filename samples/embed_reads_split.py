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

name = 'ag'
a = 1e-05

path_sample = argv[1]
path_model = argv[2]

fn_model_base = path_model.split('/')[-1]
k = int(fn_model_base.split('_')[1])
fn_model_base = '_'.join(fn_model_base.split('_')[1:-1])
fn_out = '%s_%s_total_kmers_split.pkl' % (name,fn_model_base)

dir_totalkmers = '/mnt/HA/groups/rosenGrp/embed_samples/data/ag/total_kmers'
path_totalkmers = os.path.join(dir_totalkmers,fn_out)

path_out = '/mnt/HA/groups/rosenGrp/embed_ag_samples/out/embeddings_split'
path_out = os.path.join(path_out,fn_model_base)
if not os.path.exists(path_out):
    os.makedirs(path_out)

r2v.embed_reads(path_sample,path_totalkmers,path_model,path_out,k=k,a=a,delim='',svm=False,
        normread=False,to_sample=True,verbose=True,v=1000)
