#!/usr/bin/env python

from itertools import product
from csv import writer

k = [12] #[6,8,10,12,15,20]
d = [64,128,256]
w = [20,50] #[10,20,50]
neg_samps = [10,20] #[5,10,20]
samp_freq = [.0001,.000001] #[.001,.00001,.0000001]
n_min = [100] #[50,100]

params = [k,d,w,neg_samps,samp_freq,n_min]
params_grid = list(product(*params))

with open('params.csv','w') as f:
    out = writer(f)
    out.writerows(params_grid)

