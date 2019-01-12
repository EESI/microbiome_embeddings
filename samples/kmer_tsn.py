import numpy as np
from sklearn.externals import joblib
from sklearn.manifold import TSNE

k = 6

[train_kmer_table, train_kmer_label] = joblib.load('ag' + str(k)  + '_mer_train.pkl')
[test_kmer_table, test_kmer_label] = joblib.load('ag' + str(k)  + '_mer_test.pkl')

X = np.concatenate((train_kmer_table, test_kmer_table), axis=0)
y = np.concatenate((train_kmer_label, test_kmer_label), axis=0)

X = X/np.sum(X, axis = 1).reshape(-1,1)
del train_kmer_table
del test_kmer_table

X_embedded = TSNE(n_components=2, n_iter = 500, n_iter_without_progress = 60).fit_transform(X)
print(X_embedded.shape)

joblib.dump([X_embedded, y], 'ag' + str(k)  + '_mer_reduced_tsn.pkl')
