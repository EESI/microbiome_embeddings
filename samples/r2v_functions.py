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

def generate_kmers(line,k,alphabet={'A':'A','C':'C','G':'G','T':'T'}):

    read = ''
    l = line.strip('\n')
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
            yield kmer

def embed_reads(path_sample,path_totalkmers,path_model,path_out,k=None,
        normread=True,to_sample=False,a=1e-5,n_components=1,
        delim=None,svm=True,verbose=True,v=1000):
    
    sample_id = os.path.basename(os.path.normpath(path_sample)).split('.')[0]

    fn_sample_base = path_model.split('/')[-1]
    fn_sample_base = '_'.join(fn_sample_base.split('_')[1:-1]) + '_' + str(a) 

    if verbose:
        print('Loading total kmers.')
    total_kmers = six.moves.cPickle.load(open(path_totalkmers,'rb'))
    
    if verbose:
        print('Loading model.')
    model = Word2Vec.load(path_model)
    d = model.layer1_size
    model = model.wv

    total_reads = sum([1 for i in open(path_sample,'r')])//2
    if verbose:
        print('Total reads in sample %s: %s.' % (sample_id,total_reads))

    if normread:
        print('Normalizing each read by total number of kmers in that read.')
    else:
        print('Normalizing each read by total number of kmers in the sample.')

    file = open(path_sample,'r')

    i = 0
    read_ids = []
    n_kmer = 0
    reads = np.zeros((d,total_reads),dtype='float64')
    for line in file:
        if line[0] == '>':
            if delim is None:
                read_id = line[1:-1]
            else:
                read_id = line[1:line.find(delim)]
            read_ids.append(read_id)
            if verbose:
                if i % v == 0:
                    print('Processing %s: %s/%s.' % (read_id,i,total_reads))
        else:
            r = np.zeros(d,dtype='float64')
            kmers = generate_kmers(line,k)
            for kmer in kmers:
                try:
                    r += model[kmer] * a/(a + total_kmers[kmer])
                    n_kmer += 1
                except KeyError:
                    continue
            reads[:,i] = r
            if normread:
                reads[:,i] /= n_kmer
                n_kmer = 0
            i += 1

    if not normread:
        reads /= n_kmer

    fn_sample = '%s_%s_remb_raw.csv.gz' % (sample_id,fn_sample_base) 
    path_out2 = os.path.join(path_out,fn_sample)
    if verbose:
        print('Saving reads to %s.' % (path_out2))
    
    if to_sample:
        df = pd.DataFrame(np.sum(reads.T,axis=0).reshape(1,-1),index=[sample_id])
        df.to_csv(path_out2,compression='gzip',header=False)
    else:
        df = pd.DataFrame(reads.T,index=read_ids)
        df.to_csv(path_out2,compression='gzip')

    if svm:
        if verbose:
            print('Performing SVD: (%s,%s).' % (d,total_reads))
        svd = TruncatedSVD(n_components=n_components, n_iter=7, random_state=0)
        svd.fit(reads)
        pc = svd.components_
        reads -= reads.dot(pc.T) * pc

        fn_sample = '%s_%s_remb.csv.gz' % (sample_id,fn_sample_base) 
        path_out2 = os.path.join(path_out,fn_sample)
        if verbose:
            print('Saving reads to %s.' % (path_out2))
        if to_sample:
            df = pd.DataFrame(np.sum(reads.T,axis=0).reshape(1,-1),index=[sample_id])
            df.to_csv(path_out2,compression='gzip',header=False)
        else:
            df = pd.DataFrame(reads.T,index=read_ids)
            df.to_csv(path_out2,compression='gzip')

def plot_embeddings(embedding, ids, taxa, taxon_level = 'genus', taxon_parent_rank = 0,
                    n = 8, min_taxa = 10, n_iter = 1000, p = None, seed = None,
                    return_df = False, verbose = True):

    if seed is not None:
        np.random.seed(seed)

    if verbose:
        print('Preparing data.')

    taxon_level = taxon_level.lower()

    taxon_levels = {'phylum':1,'class':2,'order':3,'family':4,'genus':5}
    taxon_levels_rev = {j:i for i,j in taxon_levels.items()}
    taxon_idx = taxon_levels[taxon_level]

    taxon_upper = [t[taxon_idx - 1] for i,t in taxa.items() if t[taxon_idx - 1] != 'NA']
    top_upper_taxon = collections.Counter(taxon_upper).most_common()[taxon_parent_rank][0]

    read_data = [(i,j,taxa[j][taxon_idx - 1],taxa[j][taxon_idx]) for i,j in enumerate(ids)
        if (taxa[j][taxon_idx - 1] == top_upper_taxon) & (taxa[j][taxon_idx - 1] != 'NA') & (taxa[j][taxon_idx] != 'NA')]

    top_target_taxon = collections.Counter([t2 for i,j,t1,t2 in read_data]).most_common()
    top_target_taxon = [t for t,c in top_target_taxon if c >= min_taxa]
    top_target_taxon = np.random.choice(top_target_taxon,n)

    read_data = [(i,j,t1,t2) for i,j,t1,t2 in read_data if t2 in top_target_taxon]

    embedding = np.array([embedding[:,i] for i,_,_,_ in read_data])

    if p is None:
        p = min_taxa

    if verbose:
        print('Running t-SNE with perplexity %s for %s iterations.' % (p,n_iter))

    tsne = TSNE(perplexity=p, n_components=2, init='pca', n_iter=n_iter)
    low_dim_embs = tsne.fit_transform(embedding)

    read_data = [data + tuple(data[1].split(':')) + tuple(low_dim_embs[i]) for i,data in enumerate(read_data)]

    df = pd.DataFrame(data=read_data,columns=['i', 'id', taxon_levels_rev[taxon_idx - 1], taxon_level, 'org', 'gene' ,'x','y'])

    if return_df:
        return(df)
    else:
        p = ggplot(aes(x='x',y='y',color='genus'),data=df) + geom_point() + labs(x='Axis 1',y='Axis 2')
        print(p)

def dump_vsearch_results(path_in, path_out = None, verbose = True, v=1000):

    with open(path_in, 'r') as f:

        res = {}
        current_q = ''
        rank = 1
        i = -1

        for line in f:

            items = line.replace('\n','').split('\t')

            if items[0] != current_q:

                i += 1

                res[str(i)] = {}
                current_q = items[0]
                rank = 1

                if verbose:
                    if (i+1) % v == 0:
                        print('Extracting results for read ' + str(i+1) + '.')

            if items[0] == items[1]:
                continue

            res[str(i)][str(rank)] = {'query':items[0],
                                       'match':items[1],
                                       'id':items[2],
                                       'alnlen':items[3],
                                       'mism':items[4],
                                       'opens':items[5],
                                       'qlo':items[6],
                                       'qhi':items[7],
                                       'tlo':items[8],
                                       'thi':items[9]}

            rank += 1

        if path_out is not None:
            if verbose:
                print('\nSaving results.')
            six.moves.cPickle.dump(res,open(path_out,'wb'),protocol=4)

        return res

def extract_cosin_results(idx, ids, cossim, max_rank = 50):

    res = []
    rank = 1

    for i in np.argsort(-1 * cossim[idx,]):

        if idx == i:
            continue

        res.append((ids[i], rank, cossim[idx,i]))

        if rank >= max_rank:
            break

        rank += 1

    return res

def extract_vsearch_results(idx, results, max_rank = 50):
    return [(v['match'],int(k),v['id']) for k,v in results[str(idx)].items()
               if int(k) <= max_rank]

def embed_samples(path,model,samp_ids,k=None,a=1e-4,n_components=1,
                  path_match=None,checks=False,verbose=True,v=2500):

    if path_match is not None:
      file_match = gzip.open(path_matches,'r')

    samp_ids_counts = collections.Counter(samp_ids.values())
    samp_idx = {samp:i for i,samp in enumerate(set(samp_ids.values()))}
    n_samps = len(samp_idx)

    f_kmers = [(w, model.wv.vocab[w].count) for w in model.wv.vocab]
    n_corpus = np.sum([count for _,count in f_kmers])
    f_kmers = {kmer:count/n_corpus for kmer,count in f_kmers}

    file_open = emb.open_file_method(path)
    file = file_open(path)

    wemb_samp = np.zeros((model.layer1_size,n_samps),dtype='float64')

    if verbose:
        t1 = time.time()
        t_total = 0
        print('Beginning file sweep.\n')

    wemb_read_check = 0
    wemb_samp_check = 0
    for i,line in enumerate(file):

        kmers = line.decode('utf-8').split(' ')

        if len(kmers) > 0:

            kmers = [kmer for kmer in kmers[:-1] if kmer in model.wv]

            samp_id = samp_ids[i]

            wemb_read = np.zeros(model.layer1_size,dtype='float64')
            wemb_read_count = 0

            for kmer in kmers:
              wemb_read += model.wv[kmer] * a/(a + f_kmers[kmer])
              wemb_read_count += 1

            if path_match is not None:
              match = str(wemb_read_count) + ','
              file_match.write(match.encode)

            wemb_read /= wemb_read_count
            wemb_samp[:,samp_idx[samp_id]] += wemb_read/samp_ids_counts[samp_id]

            if checks:
              wemb_read_tmp = np.max(np.abs(wemb_read))
              if wemb_read_tmp > wemb_read_check:
                  wemb_read_check = wemb_read_tmp
              wemb_samp_tmp = np.max(np.abs(wemb_samp[:,samp_idx[samp_id]]))
              if wemb_samp_tmp > wemb_samp_check:
                  wemb_samp_check = wemb_samp_tmp

            if verbose:
                if i % v == 0:
                    t2 = time.time()
                    t_diff = (t2-t1)/60
                    t_total += t_diff/60
                    print('Processed read %s/%s in %.2f minutes (Total: %.2f hours).' % (str(i),str(len(samp_ids)),t_diff,t_total))
                    if checks:
                      print('Max w-read = %.3f\nMax w-samp = %.3f\n' % (wemb_read_check, wemb_samp_check))
                    t1 = time.time()

    file_match.close()

    if verbose:
        print('Performing SVD on weighted embedding matrix (%s, %s)' % (wemb_samp.shape[0],wemb_samp.shape[1]))

    svd = TruncatedSVD(n_components=n_components, n_iter=7, random_state=0)
    svd.fit(wemb_samp)
    pc = svd.components_
    wemb_samp -= wemb_samp.dot(pc.T) * pc

    return wemb_samp

def gen_kmer_profile(path_reads,path_model,path_kmers,k,verbose=True,v=10000):
  model = Word2Vec.load(path_model)
  model_wv = model.wv
  del model

  alphabet = {'A':'A','C':'C','G':'G','T':'T'}

  open_file = emb.open_file_method(path_reads)
  in_file = open_file(path_reads)

  kmer_counter = {kmer:0 for kmer in model_wv.vocab}
  read_counter = {}
  sample_profile = {}

  v_counter = 0
  t_total = 0
  t1 = time.time()
  for line in in_file:

      l = line.decode('utf-8').strip('\n')
      if l[0] == '>':
          l = l.strip('>')
          sample = l[:l.find('_')]
      else:

          if verbose:
              if v_counter % v == 0:
                  t2 = time.time()
                  t_diff = (t2-t1)/60
                  t_total += t_diff
                  print('Processing read %s. Last batch: %.3f minutes. Total time: %.3f hours.' % (v_counter, t_diff, t_total/60))
                  t1 = time.time()
              v_counter += 1

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
                  if kmer in model_wv:
                      kmer_counter[kmer] += 1

                      try:
                          read_counter[sample] += 1
                      except:
                          read_counter[sample] = 1
                          sample_profile[sample] = {}
                          sample_profile[sample][kmer] = 1

                          continue

                      try:
                          sample_profile[sample][kmer] += 1
                      except:
                          sample_profile[sample][kmer] = 1

  total_kmers = np.sum([count for count in kmer_counter.values()])
  kmer_counter = {kmer:count/total_kmers for kmer,count in kmer_counter.items()}

  print('\nSaving results.')
  six.moves.cPickle.dump({'sample_profile':sample_profile,
                          'kmer_counter':kmer_counter,
                          'read_counter':read_counter},
                         open(path_kmers,'wb'))

def gen_kmer_profile_matrix(path_reads,path_model,path_samples,path_profile,k,verbose=True,v=10000):

    model = Word2Vec.load(path_model)
    model_wv = model.wv
    del model

    kmer_dict = {k:i for i,k in enumerate(model_wv.vocab)}

    ids_dict = six.moves.cPickle.load(open(path_samples,'rb'))
    ids_dict = {s:i for i,s in enumerate(ids_dict)}

    profile = np.zeros((len(ids_dict),len(kmer_dict)),dtype='int')

    alphabet = {'A':'A','C':'C','G':'G','T':'T'}

    in_file = open(path_reads,'r')

    v_counter = 0
    t_total = 0
    t1 = time.time()
    for line in in_file:

        l = line.strip('\n')
        if l[0] == '>':
            sample = l[1:l.find('_')]
        else:
            if verbose:
                if v_counter % v == 0:
                    t2 = time.time()
                    t_diff = (t2-t1)/60
                    t_total += t_diff
                    print('Processing read %s. Last batch: %.3f minutes. Total time: %.3f hours.' % (v_counter, t_diff, t_total/60))
                    t1 = time.time()
                v_counter += 1

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
                    if kmer in model_wv:
                        profile[ids_dict[sample],kmer_dict[kmer]] += 1

    print('\nSaving results.')
    np.save(path_profile,profile)


def gen_kmer_profile_dict(path_reads,path_model,path_out,k):

    model = Word2Vec.load(path_model)
    model_vocab = model.wv.vocab
    del model
  
    alphabet = {'A':'A','C':'C','G':'G','T':'T'}
    
    in_file = open(path_reads,'r')
  
    sample_profile = {}
  
    for line in in_file:
        l = line.strip('\n')
        if l[0] == '>':
            sample = l[1:l.find('_')]
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
                    if kmer in model_vocab:
                        try:
                            sample_profile[sample][kmer] += 1
                        except:
                            try:
                                sample_profile[sample][kmer] = 1
                            except:
                                sample_profile[sample] = collections.Counter({})
                                sample_profile[sample][kmer] = 1

    six.moves.cPickle.dump(sample_profile,open(path_out,'wb'),protocol=4)


def calc_total_kmers(path_reads,path_model,k,verbose=True,v=10000):
    
    model = Word2Vec.load(path_model)
    model_wv = model.wv
    del model
    
    kmer_counter = {kmer:0 for kmer in model_wv.vocab}

    file = open(path_reads,'r')
    v_counter = 0
    t_total = 0
    t1 = time.time()
    for line in file:
        if line[0] == '>':
            continue
        else:
            if verbose:
                if v_counter % v == 0:
                    t2 = time.time()
                    t_diff = (t2-t1)/60
                    t_total += t_diff
                    print('Processing read %s. Last batch: %.3f minutes. Total time: %.3f hours.' % (v_counter, t_diff, t_total/60))
                    t1 = time.time()
                v_counter += 1
            kmers = generate_kmers(line,k)
            for kmer in kmers:
                if kmer in model_wv:
                    kmer_counter[kmer] += 1

    total_kmers = np.sum([count for count in kmer_counter.values()])
    kmer_counter = {kmer:count/total_kmers for kmer,count in kmer_counter.items()}

    return kmer_counter


def calc_total_kmers_split(path_reads,path_model,k,verbose=True,v=10000):
    
    model = Word2Vec.load(path_model)
    model_wv = model.wv
    del model
    
    kmer_counter = defaultdict(int) #{kmer:0 for kmer in model_wv.vocab}

    file = open(path_reads,'r')
    v_counter = 0
    t_total = 0
    t1 = time.time()
    for line in file:
        if line[0] == '>':
            continue
        else:
            if verbose:
                if v_counter % v == 0:
                    t2 = time.time()
                    t_diff = (t2-t1)/60
                    t_total += t_diff
                    print('Processing read %s. Last batch: %.3f minutes. Total time: %.3f hours.' % (v_counter, t_diff, t_total/60))
                    t1 = time.time()
                v_counter += 1
            kmers = generate_kmers(line,k)
            for kmer in kmers:
                if kmer in model_wv:
                    kmer_counter[kmer] += 1

    return kmer_counter


def dictsum(*dicts):
    ret = defaultdict(int)
    for d in dicts:
        for k, v in d.items():
            ret[k] += v
    return dict(ret)
