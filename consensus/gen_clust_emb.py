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
from collections import defaultdict
import pandas as pd
from itertools import product
from operator import itemgetter
from sklearn.manifold import TSNE
import random
import math
import time

from sklearn.decomposition import TruncatedSVD
from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence

import embed_functions as emb
from glob import glob

path_work = '/mnt/HA/groups/rosenGrp/embed_cons/out/kegg'
path_params = glob(os.path.join(path_work,'embeddings_split') + '/*')

for p in path_params:
    params = p.split('/')[-1]
    path_out_raw = os.path.join(path_work,'seqs_cluster_' + params + '_1e-05_remb_raw.csv.gz')
    path_out = os.path.join(path_work,'seqs_cluster_' + params + '_1e-05_remb.csv.gz')
    path_embs = glob(p + '/*.csv.gz')
    out = np.zeros((256,len(path_embs)))
    
    for e in path_embs:
        cl = int(e.split('/')[-1].split('_')[1])
        m = np.genfromtxt(e,delimiter=',',skip_header=True,usecols=range(1,257,1))
        m = np.sum(m,axis=0)
        out[:,cl] = m

    cl_ids = ['cl_' + str(i) for i in range(len(path_embs))]
    
    df = pd.DataFrame(out.T,index=cl_ids)
    df.to_csv(path_out_raw,compression='gzip')

    svd = TruncatedSVD(n_components=1, n_iter=7, random_state=0)
    svd.fit(out)
    pc = svd.components_
    out -= out.dot(pc.T) * pc
    df = pd.DataFrame(out.T,index=cl_ids)
    df.to_csv(path_out,compression='gzip')
