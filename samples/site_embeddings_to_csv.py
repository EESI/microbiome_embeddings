#!/usr/bin/env python

import sys
from sys import argv
import os
import zipfile
from os.path import splitext, isfile
import gzip
import csv
import six.moves.cPickle

import numpy as np
import pandas as pd
import random
import math
import time
from glob import glob

import dask.array as da
import dask.dataframe as dd


dat = six.moves.cPickle.load(open('/mnt/HA/groups/rosenGrp/embed_ag_samples/out/seqs_10_256_5_50_10_1e-06_100_site_embeddings.pkl','rb'))

for site in dat:
    s = site.split(':')[1]
    df = dat[site]
    df.to_csv('seqs_10_256_5_50_10_1e-06_100_%s_embedding.csv' % (s))

