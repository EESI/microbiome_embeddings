#!/usr/bin/env python

from sys import argv
import os
import zipfile
from os.path import splitext
import gzip
import csv
import six.moves.cPickle
import collections
from smart_open import smart_open

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

import numpy as np
import pandas as pd
from itertools import product
from operator import itemgetter
from sklearn.manifold import TSNE
import random
import math
from scipy import sparse

import logging
import gensim
from gensim.models import Word2Vec, Doc2Vec
from gensim.models.word2vec import LineSentence
from gensim.models.doc2vec import TaggedDocument

import embed_functions as emb
import embed_params as p

def open_file_method(path):
    ext = splitext(path)[1]
    if (ext == '.gz'):
        def open_file(path):
            return gzip.open(path)
    else:
        def open_file(path):
            return open(path,'r')
    return open_file

def generate_kmers(k,nts='ACGT'):
    dictionary = {''.join(x):-1 for x in product(nts, repeat=k)}
    return dictionary

def extract_kmers(path_in,path_out,k,nts='ACGT',v=10000):

    assert k <= 12

    dictionary = generate_kmers(k,nts)

    unk = dict()
    ids = []
    read_idx = -1

    open_file = open_file_method(path_in)
    in_file = open_file(path_in)
    out_file = gzip.open(path_out,'w')
    for line in in_file:
        l = line.decode("utf-8").strip('\n')
        if l[0] == '>':
            ids.append(l.strip('>'))
            read_idx += 1
            if read_idx % v == 0:
                print('Processing read: ' + str(read_idx))
        else:
            read = ''
            unk_kmers = []
            M = len(l) - k + 1
            for n in range(M):
                kmer = l[n:n+k]
                try:
                    check = dictionary[kmer]
                    read += kmer + ' '
                except:
                    unk_kmers.append(kmer)
            read += '\n'
            out_file.write(read.encode())
            if len(unk_kmers) > 0:
                unk[read_idx] = unk_kmers

    out_file.close()

    return ids, unk

def extract_taxonomy(path):
    open_tax = open_file_method(path)
    tax = dict()
    file = open_tax(path)
    line = file.readline().decode("utf-8")
    id_end = line.find('\t')
    id = line[0:id_end]
    tax[id] = line[id_end+1:].strip('\n')
    for line in file:
        line = line.decode("utf-8")
        id_end = line.find('\t')
        id = line[0:id_end]
        tax[id] = line[id_end+1:].strip('\n')
    return tax

def read2vec(path,model,ids,k=None,nts='ACGT',v=2500):

    read_embeddings = np.empty((len(ids),model.layer1_size),dtype='float64')
    failures = {}
    file_open = open_file_method(path)
    file = file_open(path)

    for i,line in enumerate(file):

        kmers = line.decode("utf-8").split()
        if len(kmers) > 0:
            try:
                read_embeddings[i,] = np.mean(model.wv[kmers],axis=0)

            except KeyError:
                kmers = [k for k in kmers if k in model.wv.vocab]
                read_embeddings[i,] = np.mean(model.wv[kmers],axis=0)

        else:
            failures[str(i)] = ids[i]

        if i % v == 0:
            print('Processing read: ' + str(i) + '/' + str(len(ids)))

    return read_embeddings

def ngrams(path,ids,k=None,nts='ACGT',v=10000):

    file_open = open_file_method(path)
    file = file_open(path)

    kmer_counts = dict()
    for i, line in enumerate(file):

        kmers = line.decode('utf-8').split()
        if len(kmers) > 0:
            kmer_counts[ids[i]] = collections.Counter(kmers)

        if i % v == 0:
            print('Processing read: ' + str(i) + '/' + str(len(ids)))

    return kmer_counts

def doc2vec(path,model,ids,kmer_freq=False,k=None,nts='ACGT',v=100):

    if kmer_freq:
        kmer_i = {''.join(kmer):i for i,kmer in enumerate(product(nts, repeat=k))}
        kmer_f = np.zeros((len(ids), len(kmer_i)))

    read_embeddings = np.empty((len(ids),model.layer1_size),dtype='float64')
    failures = {}
    file_open = open_file_method(path)
    file = file_open(path)

    for i,line in enumerate(file):

        kmers = line.decode("utf-8").split()
        if len(kmers) > 0:
            try:
                read_embeddings[i,] = model.infer_vector(kmers)

                if kmer_freq:
                    kmer_idx,kmer_count = np.unique([kmer_i[kmer] for kmer in kmers], return_counts=True)
                    kmer_f[i,kmer_idx] += kmer_count/len(kmers)

            except IndexError:
                print(i)
                pass

        else:
            failures[str(i)] = ids[i]

        if i % v == 0:
            print('Processing read: ' + str(i) + '/' + str(len(ids)))

    if kmer_freq:
        kmer_f = sparse.csr_matrix(kmer_f)
        return read_embeddings, kmer_f
    else:
        return read_embeddings

def plot_read_embeddings(path,read_embeddings,ids,tax,plot_only=750,tax_level=1,col_num=8):
    tsne = TSNE(perplexity=30, n_components=2, init='pca', n_iter=5000)
    low_dim_embs = tsne.fit_transform(read_embeddings[:plot_only, :])
    label = [tax[i].split('; ')[tax_level] for i in ids][:plot_only]
    tax_counts = collections.Counter(label).most_common()
    colors = plt.cm.Set1(np.linspace(0,1,col_num))
    colors = np.concatenate((colors,np.tile([0.,0.,0.,1.],(len(tax_counts)-col_num,1))),axis=0)
    col_dict= {t[0]:colors[i] for i,t in enumerate(tax_counts)}
    sdx = 0.0001*(np.amax(low_dim_embs[:,0]) - np.amin(low_dim_embs[:,0]))
    sdy = 0.0001*(np.amax(low_dim_embs[:,1]) - np.amin(low_dim_embs[:,1]))
    plt.figure(figsize=(28, 28))
    for i, lab in enumerate(label):
        x, y = low_dim_embs[i,:]
        x += np.random.randn(1) * sdx
        y += np.random.randn(1) * sdy
        plt.scatter(x, y, c=col_dict[lab])
    plt.savefig(path)
    plt.close()
