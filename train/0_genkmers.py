#!/usr/bin/env python

from sys import argv
import os
import time
import zipfile
from os.path import splitext, isfile
import gzip
import csv
import six.moves.cPickle
import collections
from smart_open import smart_open

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

import numpy as np
import pandas as pd
from itertools import product, islice
from operator import itemgetter
from sklearn.manifold import TSNE
import random
import math
import multiprocessing as mp

import logging
import gensim
from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence

import embed_functions as emb
import embed_params as p

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

v = int(argv[3])
k = int(argv[2])
name = argv[1]  #['gg','kegg','query']

print(name + ':\tGenerating kmers.')

if name == 'gg':
  reads_fn = 'gg_13_5.fasta.gz'
if name == 'kegg':
  reads_fn = 'reference_seqs.fna.gz'
if name == 'query': # query_hmp
  reads_fn = 'seqs.fna.gz'
if name == 'query_oral':
  reads_fn = 'seqs_oral.fna.gz'

ids_fn = name + '_' + str(k) + '_ids.pkl'
kmers_fn = name + '_' + str(k) + '_kmers.csv.gz'

alphabet = {'A':'A','C':'C','G':'G','T':'T'}


print(name + ':\tLoading reads.')

open_file = emb.open_file_method(reads_fn)
in_file = open_file(reads_fn)

out_kmers = gzip.open(kmers_fn,'w')

ids = []
read_idx = 0
t1 = time.time()

for line in in_file:
  l = line.decode("utf-8").strip('\n')
  if l[0] == '>':
    ids.append(l.strip('>'))
    read_idx += 1
    if read_idx % v == 0:
      t_diff = str(round((time.time() - t1)/60,1)) + ' min.'
      print(name + ':\tProcessed ' + str(read_idx) + ' reads in ' + t_diff)
  else:
    read = ''
    M = len(l) - k + 1
    for n in range(M):
      l = list(l)
      nts = l[n:n+k]
      kmer = ''
      for nt in nts:
        try:
          kmer += alphabet[nt]
        except:
          continue
      if len(kmer) == k:
        read = read + kmer + ' '
    read += '\n'
    out_kmers.write(read.encode())
out_kmers.close()

print(name + ':\tDumping ids.')
six.moves.cPickle.dump({'ids':ids},open(ids_fn,'wb'))
