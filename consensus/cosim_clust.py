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
from sklearn.metrics.pairwise import cosine_similarity as cosim
from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence

import embed_functions as emb
from glob import glob

def cosim_search(e,q=None,qw=None):
    if (q is None) or (qw is None):
        return cosim(e,e)
    else:
        idx = [i for i,v in enumerate(qw) if v != 0]
        e = e.as_matrix()
        e_idx = e[:,idx]
        q_idx = q[idx].reshape(1,-1)
        return cosim(q_idx,e_idx).flatten()

path_work = '/mnt/HA/groups/rosenGrp/embed_cons/out/kegg'
params = ['6_256_5_50_10_1e-06_100','10_256_5_50_10_1e-06_100']
suffix = ['remb.csv.gz','remb_raw.csv.gz']

for p in params:
    path_emb = glob(path_work + '/*' + p + '*.csv.gz')

    for suf in suffix:
        mods = [path for path in path_emb if suf in path]
        if 'vcons' not in mods[0]:
            mods = [mods[1],mods[0]]
        remb_cons = pd.read_csv(mods[0],index_col=0)
        remb_clust = pd.read_csv(mods[1],index_col=0)
        sim = pd.DataFrame(cosim(remb_cons.as_matrix(),remb_clust.as_matrix()),
                index=remb_cons.index,columns=remb_clust.index)
        sim.to_csv(os.path.join(path_work,'cosim_cons_' + p + '_1e-05_' + suf))
