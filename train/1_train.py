#!/usr/bin/env python

from sys import argv, stdout
import os
from os.path import splitext, isfile
from shutil import copyfile

import logging
import gensim
from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence

assert gensim.models.doc2vec.FAST_VERSION > -1

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

param_row = int(argv[1])
n_cores = int(argv[2])
seed = int(argv[3])

params = open('params.csv','r').readlines()[param_row].rstrip().split(',')

k = int(params[0])
d = int(params[1])
w = int(params[2])
neg_samps = int(params[3])
samp_freq = float(params[4])
n_min = int(params[5])

epochs = 5
#seed = 564 #original seed

ids_fn =  'gg_' + str(k) + '_ids.pkl'
kmers_fn_in = 'gg_' + str(k) + '_kmers.csv.gz'
kmers_fn = 'gg_' + str(k) + '_kmers_' + str(param_row) + '.csv.gz'
model_fn = 'gg_' + str(k) + '_' + str(d) + \
        '_' + str(epochs) + '_' + str(w) + '_' + \
        str(neg_samps).replace('0.','') + '_' + \
        str(samp_freq) + '_' + str(n_min) + '_model.pkl'

final_dir = os.path.expanduser('~/embedding/models')

if not os.path.exists(os.path.join(final_dir,model_fn)):

    while True:
        print('Copying ' + kmers_fn_in)
        copyfile(kmers_fn_in,kmers_fn)
        if os.path.getsize(kmers_fn_in) == os.path.getsize(kmers_fn):
            break

    kmers_init = LineSentence(kmers_fn,max_sentence_length=100000)

    model = Word2Vec(kmers_init,sg=1,size=d,window=w,min_count=n_min,negative=neg_samps,
            sample=samp_freq,iter=epochs,workers=n_cores,seed=seed)

    os.remove(kmers_fn)

    model.save(os.path.join(final_dir,model_fn))

else:

    print(model_fn + ' exists')
