import gzip
import numpy as np
import pickle
import pandas as pd
from sklearn.cluster import KMeans
from sklearn import metrics
[tax_dict, tax2label] = pickle.load(open('tax_label.p','rb'))
def load_emd(filename):
    f=gzip.open(filename,'rb')
    # skip the header line
    line = f.readline()
    c = 0
    d = 0
    emd_dict = {}
    for line in f:
        med = line.decode('utf-8').split(',')
        vec = np.array([float(item) for item in med[1:]])
        if med[0] in emd_dict:
            c += 1
        emd_dict[med[0]] = vec
    f.close()
    print('sample size:', len(emd_dict))
    print('abnormal:', c)
    return emd_dict
def get_tax_given_level(tax_dict, tax2label, l):
    label_dict = {}
    for key in tax_dict:
        flag = 0
        for item in tax_dict[key]:
            level_tax = item.split('__')
            if len(level_tax) == 2:
                level = level_tax[0]
                tax = level_tax[1]
                if len(tax) != 0:
                    if level == l:
                        label_dict[key] = tax2label[level][tax]
                        flag = 1
                        break
        if flag == 0:
            label_dict[key] = -1
    return label_dict
def convert_to_mat(emd_dict, label_dict, l, filename):
    X = np.zeros((len(label_dict),4096))
    Y = np.zeros(len(label_dict))
    idx = 0
    for key in emd_dict:
        X[idx,:] = emd_dict[key]
        Y[idx] = label_dict[key]
        idx += 1
    file_name_X = filename + '_X_' + l + '6_mer.npy'
    file_name_Y = filename + '_Y_' + l + '.npy'
    np.save(file_name_X, X)
    np.save(file_name_Y, X)
    return X, Y
def KM_clustering(X, Y):
    labeled_X = X[Y!=-1,:]
    labeled_Y = Y[Y!=-1]
    n_C = np.unique(labeled_Y).shape[0]
    model = KMeans(n_clusters=n_C, random_state=0).fit(labeled_X)
    pred = model.labels_
    n_C_KM = np.unique(pred).shape[0]
    HS = metrics.homogeneity_score(labeled_Y,pred)
    CS = metrics.completeness_score(labeled_Y, pred)
    ARI = metrics.adjusted_rand_score(labeled_Y, pred)
    AMI = metrics.adjusted_mutual_info_score(labeled_Y, pred)
    NMI = metrics.normalized_mutual_info_score(labeled_Y, pred)
    return pd.Series([n_C, n_C_KM, HS, CS, ARI, AMI, NMI],index=['n_C_truth', 'n_C_KM', 
                                                                 'homogeneity', 'completeness', 
                                                                 'ARI','AMI','NMI'])
filename = 'kegg_6_mer_cleaned.csv.gz'
filename_trimed = filename[:-7]
emd_dict = load_emd(filename)
result_dict = {}
for l in ['k', 'p', 'c', 'o', 'f', 'g', 's']:
    label_dict = get_tax_given_level(tax_dict, tax2label, l)
    X, Y = convert_to_mat(emd_dict, label_dict, l, filename_trimed)
    result = KM_clustering(X, Y)
    result_dict[l] = result
df = pd.DataFrame(result_dict)
output_name = filename + '_cluster_perf.p'
pickle.dump(df, open(output_name,'wb'))
