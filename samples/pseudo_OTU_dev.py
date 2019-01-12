import pickle
import gzip
import numpy as np
from glob import glob
def get_file_list(model_name, pattern):
    model_folder = '/mnt/HA/groups/rosenGrp/embed_samples/data/ag/embeddings_split/'
    model_folder = '' # need to be commented out
    model_name = ''
    file_list = glob(model_folder + model_name + pattern)
    return file_list
def load_emd(filename):
    f=gzip.open(filename,'rb')
    # skip the header line
    line = f.readline()
    c = 0
    for line in f:
        c += 1
    f.close()
    return c
def read_emb(c, filename):
    f=gzip.open(filename,'rb')
    X = np.zeros((c,256))
    # skip the header line
    line = f.readline()
    idx = 0
    for line in f:
        med = line.decode('utf-8').split(',')
        vec = np.array([float(item) for item in med[1:]])
        X[idx,:] = vec
        idx += 1
    f.close()
    return X
    
    
    
import numpy as np
import pickle

X = pickle.load(open('sampes_for_seed_gen.p','rb'))

import time
from sklearn.cluster import MiniBatchKMeans
mbk = MiniBatchKMeans(init='k-means++', n_clusters=1000, batch_size=10000, max_iter=10,
                      n_init=10, random_state = 0, max_no_improvement=1170, verbose=1, compute_labels = False)
t0 = time.time()
mbk.fit(X)
t_mini_batch = time.time() - t0

joblib.dump(mbk, 'miniKM_model.pkl') 


import pickle
samples_info = pickle.load(open('ag_sample_label.p','rb'))

from sklearn.externals import joblib
mini_KM_otu_1= joblib.load('miniKM_otu_dict_1.pkl') 
mini_KM_otu_2= joblib.load('miniKM_otu_dict_2.pkl') 
mini_KM_otu_3= joblib.load('miniKM_otu_dict_3.pkl') 

c = 0
merged_dict = {}
for key in mini_KM_otu_1:
    if key in samples_info:
        if key in merged_dict:
            c += 1
        else:
            merged_dict[key] = mini_KM_otu_1[key]
for key in mini_KM_otu_2:
    if key in samples_info:
        if key in merged_dict:
            c += 1
        else:
            merged_dict[key] = mini_KM_otu_2[key]
for key in mini_KM_otu_3:
    if key in samples_info:
        if key in merged_dict:
            c += 1
        else:
            merged_dict[key] = mini_KM_otu_3[key]
print(len(merged_dict),c)

# raw dict without check the existence of label info
c = 0
merged_dict = {}
for key in mini_KM_otu_1:
    if key in merged_dict:
        c += 1
    else:
        merged_dict[key] = mini_KM_otu_1[key]
for key in mini_KM_otu_2:
    if key in merged_dict:
        c += 1
    else:
        merged_dict[key] = mini_KM_otu_2[key]
for key in mini_KM_otu_3:
    if key in merged_dict:
        c += 1
    else:
        merged_dict[key] = mini_KM_otu_3[key]
print(len(merged_dict),c)

header = ['x' + str(i) for i in range(1001)]
header[0] = 'Sample_ID'
header = ','.join(header)
header += '\n'
f = open('pseudo_otu.csv', 'wb')
f.write(header.encode())
for key in merged_dict:
    tmp = [key] + [str(i) for i in merged_dict[key]]
    med = ','.join(tmp)
    med += '\n'
    f.write(med.encode())
f.close()
print('Done')

import numpy as np
X = np.zeros((len(merged_dict), 1000))
Y = np.zeros(len(merged_dict))
for i, key in enumerate(merged_dict):
    X[i, :] = merged_dict[key]
    Y[i] = samples_info[key]
    
np.save('ag_emb_otu_X.npy', X)
np.save('ag_emb_otu_Y.npy',Y)

import numpy as np
X = np.load('ag_emb_otu_X.npy')
Y = np.load('ag_emb_otu_Y.npy')
X = X/np.sum(X, axis = 1).reshape(7142,1)

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.33, random_state=42)

X_test_bal = np.concatenate((X_test[y_test==0][:116],X_test[y_test==1][:116],X_test[y_test==2][:116]))
y_test_bal = np.concatenate((y_test[y_test==0][:116],y_test[y_test==1][:116],y_test[y_test==2][:116]))

from sklearn.ensemble import RandomForestClassifier
clf = RandomForestClassifier(max_depth=100, random_state=130)
clf.fit(X_train, y_train)
pred = clf.predict(X_test)
from sklearn.metrics import accuracy_score
acc = accuracy_score(y_test, pred)
acc

from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, pred)
print(cm)




